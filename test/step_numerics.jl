using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Random
using TOML
using Statistics

@testset "Bleaching fit degeneracy guard" begin
    Random.seed!(42)

    # Flat data: constant ~30 locs/frame with Poisson-like noise, 19k frames.
    # Reproduces a real-data case that produced a=-42873, k≈0 with unbounded NelderMead.
    # Expected: reject as degenerate, return nothing.
    flat = [max(0, round(Int, 30 + randn() * sqrt(30))) for _ in 1:19_000]
    @test SMLMAnalysis._estimate_bleaching_rate(flat) === nothing

    # Pure exponential decay: should recover parameters accurately (non-regression).
    t_exp = 1:5000
    exp_data = [max(0, round(Int, 10 + 50 * exp(-0.0005 * i) + randn() * 2)) for i in t_exp]
    res_exp = SMLMAnalysis._estimate_bleaching_rate(exp_data)
    @test res_exp !== nothing
    @test isapprox(res_exp.k_bleach, 5.0e-4, rtol = 0.1)
    @test isapprox(res_exp.offset, 10.0, atol = 2.0)
    @test isapprox(res_exp.N_0, 50.0, atol = 3.0)
    @test res_exp.r_squared > 0.9

    # Bleach-then-flat (realistic DNA-PAINT): should still fit with physical params.
    t_bf = 1:10_000
    bf_data = [max(0, round(Int, 20 + 30 * exp(-0.001 * i) + randn() * 1.5)) for i in t_bf]
    res_bf = SMLMAnalysis._estimate_bleaching_rate(bf_data)
    @test res_bf !== nothing
    @test res_bf.offset >= 0     # physical bound held
    @test res_bf.N_0 >= 0        # physical bound held
    @test res_bf.k_bleach > 0
end

@testset "intensity-filter p₂ mixture estimator" begin
    # Regression for the unbiased double-emitter fraction estimator.
    # The mixture model must (a) recover a small KNOWN p₂ and, critically,
    # (b) return ≈0 when there are NO doubles — NOT the ~(1 - rate_percentile)
    # ≈ 5% floor that the legacy tail-ratio estimate is pinned at.
    rng = MersenneTwister(20260717)
    # Peaked, identifiable single-emitter model: Gamma(k=3, θ=1) drawn as the sum
    # of three exponentials. Field-normalized by the single-emitter p95 (→ ~1),
    # mirroring the pipeline's photons/λ(x,y).
    gamma3(r) = -log(rand(r)) - log(rand(r)) - log(rand(r))
    function make_normalized(p_true, N, r)
        n_dbl = round(Int, N * p_true)
        n_sgl = N - n_dbl
        raw = [gamma3(r) for _ in 1:(n_sgl + 2 * n_dbl)]
        scale = quantile(raw, 0.95)
        scale <= 0 && (scale = 1.0)
        singles = raw[1:n_sgl] ./ scale
        da = raw[(n_sgl + 1):(n_sgl + n_dbl)] ./ scale
        db = raw[(n_sgl + n_dbl + 1):(n_sgl + 2 * n_dbl)] ./ scale
        vcat(singles, da .+ db)   # f_double = f_single ⊛ f_single (sum of two singles)
    end
    cfg = IntensityFilterConfig(p2_n_bins = 200)
    @test cfg.p2_method === :mixture   # new behavior is the default

    # (a) p_true = 1% → recovered within a few × and clearly off the 5% floor.
    p_hat = SMLMAnalysis._estimate_p2_mixture(make_normalized(0.01, 60_000, rng), cfg)
    @test p_hat !== nothing
    @test 0.003 <= p_hat <= 0.03

    # (b) p_true = 0% → near zero, NOT pinned near 1 - rate_percentile (~0.05).
    p_hat0 = SMLMAnalysis._estimate_p2_mixture(make_normalized(0.0, 60_000, rng), cfg)
    @test p_hat0 !== nothing
    @test p_hat0 < 0.02
end

@testset "intensity-filter spatial binning matches the reference scan" begin
    # Reference implementation: origin/main's _spatial_bin_rates, an O(n_bins²·n)
    # per-bin scan against the exact half-open edges (see git show
    # origin/main:src/steps/intensityfilter.jl). SMLMAnalysis._spatial_bin_rates
    # computes each emitter's bin directly via _bin_index instead of scanning, so
    # this pins the single-pass version to the scan's exact edge comparisons —
    # including the float-rounding cases a floor-division guess gets wrong.
    function ref_spatial_bin_rates(xs, ys, photons, cfg)
        n_bins = cfg.n_bins
        min_count = cfg.min_bin_count
        pct = cfg.rate_percentile

        x_min, x_max = extrema(xs)
        y_min, y_max = extrema(ys)

        dx = (x_max - x_min) / n_bins
        dy = (y_max - y_min) / n_bins
        dx == 0 && (dx = 1.0)
        dy == 0 && (dy = 1.0)

        centers_x = Float64[]
        centers_y = Float64[]
        rates = Float64[]
        counts = Int[]
        rate_grid = fill(NaN, n_bins, n_bins)

        for ix in 1:n_bins, iy in 1:n_bins
            bx_lo = x_min + (ix - 1) * dx
            bx_hi = ix == n_bins ? x_max + eps(x_max) : x_min + ix * dx
            by_lo = y_min + (iy - 1) * dy
            by_hi = iy == n_bins ? y_max + eps(y_max) : y_min + iy * dy

            bin_photons = Float64[]
            for i in eachindex(xs)
                if bx_lo <= xs[i] < bx_hi && by_lo <= ys[i] < by_hi
                    push!(bin_photons, photons[i])
                end
            end

            if length(bin_photons) >= min_count
                rate = quantile(bin_photons, pct)
                push!(centers_x, (bx_lo + min(bx_hi, x_max)) / 2)
                push!(centers_y, (by_lo + min(by_hi, y_max)) / 2)
                push!(rates, rate)
                push!(counts, length(bin_photons))
                rate_grid[ix, iy] = rate
            end
        end

        x_edges = range(x_min, x_max, length = n_bins + 1)
        y_edges = range(y_min, y_max, length = n_bins + 1)

        return (
            centers_x = centers_x, centers_y = centers_y, rates = rates, counts = counts,
            rate_grid = rate_grid, x_edges = x_edges, y_edges = y_edges,
        )
    end

    # rate_grid carries NaN in unfilled bins; isequal (not ==) is the correct
    # "identical field-by-field" check since NaN == NaN is false but
    # isequal(NaN, NaN) is true — a bare == would fail even when the two
    # implementations agree exactly.
    function check_matches_reference(xs, ys, photons, cfg)
        ref = ref_spatial_bin_rates(xs, ys, photons, cfg)
        got = SMLMAnalysis._spatial_bin_rates(xs, ys, photons, cfg)
        @test isequal(ref, got)
    end

    # (a) Interior-edge float rounding: extrema (1.0, 2.0), 10 bins → dx=0.1;
    # x=1.2 lands in bin 3 under exact edge comparison, but a naive
    # floor((x-x_min)/dx)+1 guess can drift to bin 2 on this input.
    cfg_a = IntensityFilterConfig(n_bins = 10, min_bin_count = 1)
    xs_a = [1.0, 1.2, 2.0, 1.5, 1.3, 1.7, 1.9, 1.05, 1.85, 1.45]
    ys_a = collect(range(0.0, 1.0, length = 10))
    photons_a = collect(100.0:100.0:1000.0)
    check_matches_reference(xs_a, ys_a, photons_a, cfg_a)
    @test SMLMAnalysis._bin_index(1.2, 1.0, 2.0, 0.1, 10) == 3

    # (b) Constant x: x_min == x_max forces dx = 1.0 (the /0 guard), which must
    # still clip bin centers to min(bx_hi, x_max).
    cfg_b = IntensityFilterConfig(n_bins = 4, min_bin_count = 1)
    xs_b = fill(5.0, 12)
    ys_b = collect(range(0.0, 3.0, length = 12))
    photons_b = collect(10.0:10.0:120.0)
    check_matches_reference(xs_b, ys_b, photons_b, cfg_b)

    # (c) Constant y: same, on the other axis.
    xs_c = collect(range(0.0, 3.0, length = 12))
    ys_c = fill(-2.0, 12)
    check_matches_reference(xs_c, ys_c, photons_b, cfg_b)

    # (d) Points exactly on interior edges x_min + k*dx (k = 0..n_bins), plus a
    # repeated x_max point ((e) the max point) to confirm the widened last-bin
    # edge still claims it rather than dropping it as out-of-range.
    cfg_d = IntensityFilterConfig(n_bins = 10, min_bin_count = 1)
    xs_d = vcat(collect(0.0:1.0:10.0), 10.0)   # 0,1,...,10, and a second 10.0
    ys_d = collect(range(0.0, 5.0, length = length(xs_d)))
    photons_d = collect(1.0:length(xs_d))
    check_matches_reference(xs_d, ys_d, photons_d, cfg_d)

    # (f) ~20 random seeded clouds.
    rng_bins = MersenneTwister(20260926)
    for _ in 1:20
        n = rand(rng_bins, 40:200)
        xs_r = rand(rng_bins, n) .* 10 .- 3
        ys_r = rand(rng_bins, n) .* 6 .+ 1
        photons_r = rand(rng_bins, n) .* 900 .+ 100
        cfg_r = IntensityFilterConfig(n_bins = rand(rng_bins, (4, 5, 8, 10)), min_bin_count = 1)
        check_matches_reference(xs_r, ys_r, photons_r, cfg_r)
    end
end
