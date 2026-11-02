using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Random
using TOML
using Statistics
using Logging

@testset "Density auto threshold preserves a sparse uniform population" begin
    Random.seed!(42)
    cam = IdealCamera(100, 100, 0.1)
    # Same pair radius as DensityFilterConfig(): 2 * hypot(hypot(0.1, 0.1), ...).
    radius = 0.4
    n = round(Int, 3.0 * 100 / (π * radius^2))
    emitters = [
        SMLMAnalysis.Emitter2DFit{Float64}(
            10 * rand(), 10 * rand(), 1000.0, 10.0, 0.1, 0.1, 30.0, 1.0
        ) for _ in 1:n
    ]
    counts = [
        count(j -> j != i && hypot(emitters[i].x - emitters[j].x, emitters[i].y - emitters[j].y) < radius, 1:n)
            for i in 1:n
    ]
    mode_count = argmax([count(==(k), counts) for k in 0:maximum(counts)]) - 1
    @test 2 <= mean(counts) <= 4
    smld = SMLMAnalysis.BasicSMLD(emitters, cam, 1, 1, Dict{String, Any}())
    filtered, info = analyze(smld, DensityFilterConfig(); verbose = Verbosity.SILENT)
    @test length(filtered.emitters) / n > 0.9
    @test info.info.threshold <= mode_count
end

@testset "Filter diagnostics without p-values" begin
    Random.seed!(53)
    cam = IdealCamera(16, 16, 0.1)
    for E in (SMLMAnalysis.Emitter2DFit{Float64}, SMLMAnalysis.Emitter3DFit{Float64})
        emitters = [
            if E <: SMLMAnalysis.Emitter3DFit
                E(
                    0.8, 0.8, 0.1 * randn(), 500.0 + 1000 * rand(), 10.0 + rand(),
                    0.005 + 0.01 * rand(), 0.005 + 0.01 * rand(), 0.02 + 0.01 * rand(),
                    30.0, 1.0; frame = i
                )
            else
                E(
                    0.8, 0.8, 500.0 + 1000 * rand(), 10.0 + rand(),
                    0.005 + 0.01 * rand(), 0.005 + 0.01 * rand(), 30.0, 1.0; frame = i
                )
            end for i in 1:12
        ]
        smld = SMLMAnalysis.BasicSMLD(emitters, cam, 12, 1, Dict{String, Any}())
        for verbose in (Verbosity.STANDARD, Verbosity.DEBUG),
                cfg in (
                    FilterConfig(photons = (1000.0, Inf)),
                    FilterConfig(precision = (0.0, 0.01)), FilterConfig(pvalue = (0.001, 1.0)),
                )
            mktempdir() do outdir
                # A real cache exercises the overlay path, including box colors and labels.
                batch = SMLMAnalysis.ROIBatch(fill(20.0f0, 7, 7, 1), [5], [5], [1], cam)
                SMLMAnalysis._save_detectfit_sample_cache(
                    outdir, smld, fill(20.0f0, 16, 16, 1), batch, [1]
                )
                filtered, info = if cfg.pvalue !== nothing
                    @test_logs (:warn, r"FilterConfig.pvalue has no effect:.*no p-value") min_level = Logging.Warn analyze(
                        smld, cfg; outdir = outdir, verbose = verbose, checkpoint = Checkpoint.NONE
                    )
                else
                    analyze(smld, cfg; outdir = outdir, verbose = verbose, checkpoint = Checkpoint.NONE)
                end
                expected = if cfg.photons !== nothing
                    filter(e -> e.photons >= 1000.0, emitters)
                elseif cfg.precision !== nothing
                    filter(e -> max(e.σ_x, e.σ_y) <= 0.01, emitters)
                else
                    emitters
                end
                @test filtered.emitters == expected
                @test info.info.n_after == length(expected)
                dir = joinpath(outdir, "00_filter")
                @test isfile(joinpath(dir, "fit_quality.png"))
                @test isfile(joinpath(dir, "fit_overlay.png"))
                @test isfile(joinpath(dir, "stats.md"))
                @test !occursin("- pvalue:", read(joinpath(dir, "stats.md"), String))
                if verbose == Verbosity.DEBUG
                    detailed = read(joinpath(dir, "detailed_stats.md"), String)
                    @test !occursin("| P-value", detailed)
                end
                @test SMLMAnalysis._fit_box_color(emitters[1], FilterConfig(pvalue = (0.001, 1.0)), nothing) == :green
            end
        end
        @test !occursin("pval", SMLMAnalysis._fit_overlay_title(FilterConfig(pvalue = (0.001, 1.0)), emitters))
        # Empty vectors must also pass every writer without inspecting emitters[1].
        empty_smld = SMLMAnalysis.BasicSMLD(E[], cam, 12, 1, Dict{String, Any}())
        mktempdir() do outdir
            filtered, _ = analyze(
                empty_smld, FilterConfig(pvalue = (0.001, 1.0));
                outdir = outdir, verbose = Verbosity.DEBUG, checkpoint = Checkpoint.NONE
            )
            @test isempty(filtered.emitters)
        end
    end
end

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

@testset "Driftcorrect max drift: reference frame by dataset_mode" begin
    # (a) Dataset 2 sits 5 um off (inter shift) with a P1 intra term from -0.1 to +0.1 um.
    # :registered measures from dataset 2's own first frame, so the 5 um shift cancels and
    # the excursion is 200 nm. _calc_inter_shifts reports the shift separately.
    m = SMLMDriftCorrection.LegendrePolynomial(2, 2, 100; degree = 1)
    m.inter[2].dm[1] = 5.0
    m.intra[2].dm[1].coefficients = [0.1]
    @test SMLMAnalysis._calc_max_drift(m, :registered) ≈ 200.0 atol = 1.0e-6
    @test SMLMAnalysis._calc_inter_shifts(m)[2] ≈ 5000.0 atol = 1.0e-6

    # (b) Two chunks chained as one acquisition: chunk 1 runs 0 to 0.2 um and chunk 2 runs
    # 0.2 to 0.4 um. :continuous measures from the acquisition's first frame, so the whole
    # 400 nm is drift. :registered measures each chunk from its own first frame, so it
    # sees two 200 nm excursions and the 0.2 um step between them is excluded.
    mc = SMLMDriftCorrection.LegendrePolynomial(2, 2, 100; degree = 1)
    mc.inter[1].dm[1] = 0.1
    mc.intra[1].dm[1].coefficients = [0.1]
    mc.inter[2].dm[1] = 0.3
    mc.intra[2].dm[1].coefficients = [0.1]
    @test SMLMAnalysis._calc_max_drift(mc, :continuous) ≈ 400.0 atol = 1.0e-6
    @test SMLMAnalysis._calc_max_drift(mc, :registered) ≈ 200.0 atol = 1.0e-6

    # (c) One dataset whose trajectory starts at zero: both modes give the old
    # drift_trajectory excursion (max |x| over the trajectory).
    m1 = SMLMDriftCorrection.LegendrePolynomial(2, 1, 100; degree = 1)
    m1.inter[1].dm[1] = 0.1
    m1.intra[1].dm[1].coefficients = [0.1]
    old = maximum(abs.(SMLMDriftCorrection.drift_trajectory(m1; dataset = 1).x)) * 1000
    @test old ≈ 200.0 atol = 1.0e-6
    @test SMLMAnalysis._calc_max_drift(m1, :continuous) ≈ old atol = 1.0e-6
    @test SMLMAnalysis._calc_max_drift(m1, :registered) ≈ old atol = 1.0e-6
end
