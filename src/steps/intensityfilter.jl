"""
Intensity filter step - rejects multi-emitter events via Poisson upper-tail test.

When two fluorophores are simultaneously active within one diffraction-limited spot,
the fitter produces a single localization with corrupted position and abnormally high
photon count. This step estimates the spatially-varying expected single-emitter rate
using a high percentile of the per-bin photon distribution, then fits a smooth field
model and applies a Poisson upper-tail test to reject statistical outliers.

The per-bin percentile captures the upper bound of single-emitter emission at each
spatial position, accounting for both the excitation field profile and on-time variance.
"""

"""
    IntensityFilterConfig <: AbstractSMLMConfig

Intensity-based multi-emitter rejection. Estimates the spatially-varying expected
single-emitter photon rate and rejects localizations whose photon count is statistically
inconsistent with single-emitter emission (Poisson upper-tail test).

The rate estimate uses the `rate_percentile` of the per-bin photon distribution
(default: 95th percentile). This captures the upper bound of single-emitter emission
at each spatial position, accounting for both excitation field variation and stochastic
on-time variance. The Gaussian field fit then smooths across bins.

# Keywords
- `cutoff::Float64`: P-value cutoff; emitters with p < cutoff are rejected (default: 0.01)
- `field_mode::Symbol`: Excitation field model — `:uniform` (single global rate) or
  `:gaussian` (spatially-varying Gaussian beam fit) (default: `:gaussian`)
- `n_bins::Int`: Spatial grid bins per axis for field estimation (default: 10)
- `min_bin_count::Int`: Minimum emitters per bin for rate estimation (default: 30)
- `rate_percentile::Float64`: Percentile of per-bin photon distribution used as the
  expected single-emitter rate (default: 0.95). Higher = more permissive.

# p₂ (double-emitter fraction, info-only)
`p₂` is the estimated fraction of localizations that are actually two overlapping
emitters. It is reported in `IntensityFilterInfo` and does **not** affect which
localizations are filtered.

- `estimate_p2::Bool`: Whether to estimate p₂ (default: `true`).
- `p2_method::Symbol`: Estimator for p₂ (default: `:mixture`).
    - `:mixture` — threshold-free two-component mixture fit of the field-normalized
      photon distribution: `h(x) = (1-p)·f_single(x) + p·f_double(x)` with
      `f_double = f_single ⊛ f_single`. `f_single` is a Gamma model fit jointly with
      the mixing weight `p`; `p` is returned directly as the double fraction. Uses the
      whole distribution (no tail threshold), so it is unbiased and does not floor near
      `1 - rate_percentile`.
    - `:tail` — legacy tail-ratio estimate `Σh(>τ) / Σf_double(>τ)` at `p2_tail_threshold`.
      Biased high (floored near `1 - rate_percentile`) because it attributes the entire
      single-emitter tail to doubles. Kept for comparison/back-compat only.
- `p2_tail_threshold::Float64`: Tail threshold τ (in units of λ) for the `:tail` method
  and for the tail-mass diagnostics reported in `IntensityFilterInfo` (default: 1.0).
- `p2_n_bins::Int`: Number of histogram bins used for p₂ estimation (default: 200).
"""
@kwdef struct IntensityFilterConfig <: SMLMData.AbstractSMLMConfig
    cutoff::Float64 = 0.01
    field_mode::Symbol = :gaussian
    n_bins::Int = 10
    min_bin_count::Int = 30
    rate_percentile::Float64 = 0.95
    estimate_p2::Bool = true
    p2_method::Symbol = :mixture
    p2_tail_threshold::Float64 = 1.0
    p2_n_bins::Int = 200
end

"""
    intensityfilter_step(smld, cfg; outdir=nothing, step_number=0, verbose=Verbosity.STANDARD)

Filter localizations by intensity-based multi-emitter rejection. Returns `(filtered_smld, IntensityFilterInfo)`.
"""
function intensityfilter_step(smld::BasicSMLD, cfg::IntensityFilterConfig;
                               outdir::Union{String,Nothing}=nothing,
                               step_number::Int=0,
                               verbose::Int=Verbosity.STANDARD)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)

    v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg))" cutoff=cfg.cutoff field_mode=cfg.field_mode rate_percentile=cfg.rate_percentile

    emitters = smld.emitters
    n_before = length(emitters)

    # Edge case: too few emitters to estimate field
    if n_before < 100
        v >= Verbosity.PROGRESS && @warn "  → Too few emitters ($n_before < 100), skipping intensity filter"
        info = IntensityFilterInfo(n_before, n_before, cfg.field_mode, 0.0, nothing, 0.0, 0.0,
                                    nothing, cfg.p2_tail_threshold, nothing, nothing)
        return (smld, info)
    end

    t_total = @elapsed begin
        xs = [e.x for e in emitters]
        ys = [e.y for e in emitters]
        photons = [e.photons for e in emitters]

        # Global rate estimate (high percentile = upper bound of single-emitter emission)
        lambda_global = max(1.0, quantile(photons, cfg.rate_percentile))

        # Estimate excitation field from per-bin percentiles
        actual_mode, field_params, field_r2, field_model = _estimate_excitation_field(
            xs, ys, photons, lambda_global, cfg, v)

        # Compute p-values and filter
        pvals = _compute_intensity_pvalues(xs, ys, photons, field_model)
        keep = pvals .>= cfg.cutoff

        filtered_emitters = emitters[keep]
        filtered_smld = BasicSMLD(filtered_emitters, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)
    end

    n_after = length(filtered_emitters)
    n_rejected = n_before - n_after

    # Estimate double-emitter fraction via mixture decomposition
    p2_est, p2_tail_obs, p2_tail_f2 = if cfg.estimate_p2
        _estimate_p2(xs, ys, photons, field_model, cfg)
    else
        (nothing, nothing, nothing)
    end

    info = IntensityFilterInfo(n_before, n_after, actual_mode, lambda_global, field_params, field_r2, t_total,
                                p2_est, cfg.p2_tail_threshold, p2_tail_obs, p2_tail_f2)

    if dir !== nothing
        _save_intensityfilter_outputs!(dir, cfg, v, info, xs, ys, photons, pvals, keep, field_model)
    end

    p2_str = p2_est !== nothing ? ", p₂=$(round(100*p2_est, digits=2))%" : ""
    v >= Verbosity.PROGRESS && @info "  → $n_rejected rejected ($actual_mode field, λ_global=$(round(lambda_global, digits=0)), R²=$(round(field_r2, digits=3))$p2_str) ($(round(t_total, digits=2))s)"
    (filtered_smld, info)
end

function _step_summary(info::IntensityFilterInfo)
    d = Dict{Symbol,Any}(
        :n_before => info.n_before,
        :n_after => info.n_after,
        :n_rejected => info.n_before - info.n_after,
        :field_mode => info.field_mode,
        :lambda_max_global => round(info.lambda_max_global, digits=1),
        :field_fit_r2 => round(info.field_fit_r2, digits=3),
    )
    if info.p2_estimate !== nothing
        d[:p2_estimate] = round(info.p2_estimate, sigdigits=3)
        d[:p2_tail_obs] = round(info.p2_tail_obs, sigdigits=3)
        d[:p2_tail_f2] = round(info.p2_tail_f2, sigdigits=3)
    end
    d
end

"""
    analyze(smld, cfg::IntensityFilterConfig; kwargs...) -> (filtered_smld, StepInfo)

Filter localizations by intensity-based multi-emitter rejection.
"""
function analyze(smld::BasicSMLD, cfg::IntensityFilterConfig;
                 outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD,
                 checkpoint::Int=Checkpoint.EXPENSIVE, kwargs...)
    t = @elapsed (filtered, if_info) = intensityfilter_step(smld, cfg;
        outdir=outdir, step_number=step_number, verbose=verbose)

    if checkpoint >= Checkpoint.ALL
        dir = step_outdir(outdir, step_number, cfg)
        _save_step_smld(dir, filtered; filename="smld_intensity.jld2")
    end

    (filtered, StepInfo(step_number, cfg, t, _step_summary(if_info); info=if_info))
end

# ============================================================
# Excitation field estimation
# ============================================================

"""
    _estimate_excitation_field(xs, ys, photons, lambda_global, cfg, verbose)

Estimate excitation field from per-bin percentiles of photon counts.

Each spatial bin computes `quantile(photons, rate_percentile)` independently, capturing
the local upper bound of single-emitter emission. A Gaussian field model is then fit
to these per-bin rates to produce a smooth `field_model(x, y)`.

Returns `(actual_mode, field_params, r2, field_model)`.
"""
function _estimate_excitation_field(xs, ys, photons, lambda_global, cfg::IntensityFilterConfig, verbose)
    if cfg.field_mode === :uniform
        return :uniform, nothing, 0.0, (x, y) -> lambda_global
    end

    # Spatial binning — compute percentile per bin
    bin_result = _spatial_bin_rates(xs, ys, photons, cfg)

    # Need enough valid bins for Gaussian fit
    n_valid = length(bin_result.rates)
    if n_valid < 5
        verbose >= Verbosity.PROGRESS && @warn "  → Only $n_valid valid bins, falling back to uniform field"
        return :uniform, nothing, 0.0, (x, y) -> lambda_global
    end

    # Fit Gaussian field model to per-bin rates
    params, r2 = _fit_gaussian_field(bin_result.centers_x, bin_result.centers_y, bin_result.rates, xs, ys)

    if params === nothing || r2 < 0.3
        verbose >= Verbosity.PROGRESS && @warn "  → Gaussian field fit poor (R²=$(round(r2, digits=3))), falling back to uniform"
        return :uniform, nothing, 0.0, (x, y) -> lambda_global
    end

    A, x0, y0, w, bg = params.A, params.x0, params.y0, params.w, params.bg
    field_model = (x, y) -> max(1.0, A * exp(-2.0 * ((x - x0)^2 + (y - y0)^2) / w^2) + bg)

    return :gaussian, params, r2, field_model
end

"""
    _compute_bin_stat(xs, ys, values, cfg, stat_fn) -> Matrix{Float64}

Bin emitters spatially and compute `stat_fn(bin_values)` per bin.
Returns `n_bins × n_bins` grid in Makie convention: `grid[ix, iy]` maps to
`(x_centers[ix], y_centers[iy])` when used with `heatmap!(ax, x_centers, y_centers, grid)`.

This is the fundamental convention: first array index = x-axis (horizontal),
second = y-axis (vertical). Julia arrays are (row, col) = (y, x) in image
convention, but Makie's heatmap maps first index to x. So we store as [ix, iy].
"""
function _compute_bin_stat(xs, ys, values, cfg::IntensityFilterConfig, stat_fn)
    n_bins = cfg.n_bins
    min_count = cfg.min_bin_count

    x_min, x_max = extrema(xs)
    y_min, y_max = extrema(ys)

    dx = (x_max - x_min) / n_bins
    dy = (y_max - y_min) / n_bins
    dx == 0 && (dx = 1.0)
    dy == 0 && (dy = 1.0)

    grid = fill(NaN, n_bins, n_bins)

    for ix in 1:n_bins, iy in 1:n_bins
        bx_lo = x_min + (ix - 1) * dx
        bx_hi = ix == n_bins ? x_max + eps(x_max) : x_min + ix * dx
        by_lo = y_min + (iy - 1) * dy
        by_hi = iy == n_bins ? y_max + eps(y_max) : y_min + iy * dy

        bin_vals = Float64[]
        for i in eachindex(xs)
            if bx_lo <= xs[i] < bx_hi && by_lo <= ys[i] < by_hi
                push!(bin_vals, values[i])
            end
        end

        if length(bin_vals) >= min_count
            grid[ix, iy] = stat_fn(bin_vals)
        end
    end

    grid
end

"""
    _spatial_bin_rates(xs, ys, photons, cfg)

Bin emitters spatially and compute `rate_percentile` of photons per bin.
Returns named tuple with `centers_x`, `centers_y`, `rates`, `counts`, and grid info.
"""
function _spatial_bin_rates(xs, ys, photons, cfg::IntensityFilterConfig)
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
    # Grid stored as [ix, iy] so heatmap(x, y, grid) maps correctly
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

    x_edges = range(x_min, x_max, length=n_bins+1)
    y_edges = range(y_min, y_max, length=n_bins+1)

    return (centers_x=centers_x, centers_y=centers_y, rates=rates, counts=counts,
            rate_grid=rate_grid, x_edges=x_edges, y_edges=y_edges)
end

"""
    _fit_gaussian_field(cx, cy, rates, xs, ys)

Fit `N(x,y) = A·exp(-2·((x-x₀)² + (y-y₀)²)/w²) + bg` to per-bin rates.

Returns `(params::NamedTuple, r2::Float64)` or `(nothing, 0.0)` on failure.
"""
function _fit_gaussian_field(cx, cy, rates, xs, ys)
    brightest = argmax(rates)
    A0 = maximum(rates) - minimum(rates)
    x0_0 = cx[brightest]
    y0_0 = cy[brightest]
    fov_diag = sqrt((maximum(xs) - minimum(xs))^2 + (maximum(ys) - minimum(ys))^2)
    w0 = fov_diag / 2
    bg0 = minimum(rates)

    function cost(p)
        A, x0, y0, w, bg = p
        w <= 0 && return 1e12
        s = 0.0
        for i in eachindex(cx)
            pred = A * exp(-2.0 * ((cx[i] - x0)^2 + (cy[i] - y0)^2) / w^2) + bg
            s += (rates[i] - pred)^2
        end
        s
    end

    result = try
        optimize(cost, [A0, x0_0, y0_0, w0, bg0], NelderMead(),
                 Optim.Options(iterations=5000, g_tol=1e-8))
    catch err
        err isa InterruptException && rethrow()
        @warn "intensityfilter: Gaussian excitation-field fit failed; falling back to uniform field" exception=err
        return nothing, 0.0
    end

    p = Optim.minimizer(result)
    A, x0, y0, w, bg = p

    ss_res = 0.0
    mean_rate = mean(rates)
    ss_tot = 0.0
    for i in eachindex(cx)
        pred = A * exp(-2.0 * ((cx[i] - x0)^2 + (cy[i] - y0)^2) / w^2) + bg
        ss_res += (rates[i] - pred)^2
        ss_tot += (rates[i] - mean_rate)^2
    end
    r2 = ss_tot > 0 ? 1.0 - ss_res / ss_tot : 0.0

    params = (A=A, x0=x0, y0=y0, w=w, bg=bg)
    return params, r2
end

# ============================================================
# P-value computation
# ============================================================

"""
    _compute_intensity_pvalues(xs, ys, photons, field_model)

Compute Poisson upper-tail p-values for each emitter.
`p = P(X ≥ N)` where `X ~ Poisson(λ)` and `λ = field_model(x, y)`.

λ is the per-bin percentile (upper bound of single-emitter emission at this position),
so emitters near or below λ get p ≈ 1 (kept), while emitters far above λ get small p.
"""
function _compute_intensity_pvalues(xs, ys, photons, field_model)
    n = length(xs)
    pvals = Vector{Float64}(undef, n)
    for i in 1:n
        λ = field_model(xs[i], ys[i])
        if λ <= 0
            pvals[i] = 1.0
        else
            λ = max(1.0, λ)
            N = round(Int, photons[i])
            N <= 0 && (pvals[i] = 1.0; continue)
            pvals[i] = ccdf(Poisson(λ), N - 1)  # P(X ≥ N)
        end
    end
    pvals
end

# ============================================================
# Double-emitter rate estimation (p₂)
# ============================================================

"""
    _convolve_1d(h::Vector{Float64}) -> Vector{Float64}

Self-convolution of histogram `h`. Returns vector of length `2n-1` where `n = length(h)`.
The result represents the distribution of the sum of two independent draws from `h`.
"""
function _convolve_1d(h::Vector{Float64})
    n = length(h)
    out = zeros(Float64, 2n - 1)
    @inbounds for i in 1:n, j in 1:n
        out[i + j - 1] += h[i] * h[j]
    end
    out
end

"""
    _field_normalized_photons(xs, ys, photons, field_model) -> Vector{Float64}

Field-normalized photon counts `n_i = photons_i / max(1, λ(x_i, y_i))`. Under a
well-fit field model, single emitters cluster near 1 and genuine doubles near 2.
"""
function _field_normalized_photons(xs, ys, photons, field_model)
    n = length(xs)
    normalized = Vector{Float64}(undef, n)
    for i in 1:n
        λ = field_model(xs[i], ys[i])
        normalized[i] = photons[i] / max(1.0, λ)
    end
    normalized
end

"""
    _p2_tail_masses(normalized, cfg) -> (tail_obs, tail_f2)

Diagnostic tail masses at `τ = cfg.p2_tail_threshold`: the observed mass above τ and
the self-convolution (double-emitter) mass above τ, using a `[0, p99.5]` histogram with
`cfg.p2_n_bins` bins. These populate `IntensityFilterInfo` and back the legacy `:tail`
estimator. Returns `(nothing, nothing)` on failure.
"""
function _p2_tail_masses(normalized, cfg::IntensityFilterConfig)
    upper = quantile(normalized, 0.995)
    upper <= 0 && return (nothing, nothing)

    n_bins = cfg.p2_n_bins
    bin_width = upper / n_bins
    τ = cfg.p2_tail_threshold

    h = zeros(Float64, n_bins)
    n_in_range = 0
    for v in normalized
        if 0 <= v < upper
            idx = clamp(floor(Int, v / bin_width) + 1, 1, n_bins)
            h[idx] += 1.0
            n_in_range += 1
        end
    end
    n_in_range == 0 && return (nothing, nothing)
    h ./= n_in_range

    f2 = _convolve_1d(h)

    τ_bin = clamp(floor(Int, τ / bin_width) + 1, 1, n_bins)
    tail_obs = sum(h[τ_bin:end])
    τ_bin_f2 = clamp(floor(Int, τ / bin_width) + 1, 1, length(f2))
    tail_f2 = sum(f2[τ_bin_f2:end])

    return (tail_obs, tail_f2)
end

"""
    _estimate_p2_mixture(normalized, cfg) -> Union{Float64, Nothing}

Threshold-free double-emitter fraction via a two-component mixture fit of the
field-normalized photon distribution:

    h(x) = (1 - p)·f_single(x) + p·f_double(x),   f_double = f_single ⊛ f_single

`f_single` is modelled as a Gamma density (2 shape/scale parameters); `f_double` is its
numerical self-convolution (`_convolve_1d`). The three parameters `(k, θ, p)` are fit by
maximum likelihood (multinomial NLL of the binned counts) via Nelder-Mead, and `p` — the
mixing weight — is returned directly as the double fraction, clamped to `[0, 1]`.

Unlike the legacy tail-ratio estimate this uses the whole distribution: the Gamma tail is
pinned by the bulk, so any excess near `x ≈ 2` (where singles are scarce) is attributed to
doubles rather than assumed away. It therefore removes the `~1 - rate_percentile` floor and
returns ≈0 when there are no doubles.

Note: identifiability requires a peaked single-emitter distribution (mode away from 0). For
a near-exponential single (Gamma shape ≈ 1) a small double fraction is partly absorbable by
the Gamma shape and may be under-estimated. Returns `nothing` if the fit fails.
"""
function _estimate_p2_mixture(normalized, cfg::IntensityFilterConfig)
    x = filter(v -> isfinite(v) && v > 0, normalized)
    length(x) < 100 && return nothing

    # Histogram grid extended past the double regime (~2×) so f_double fits on-grid.
    upper = clamp(quantile(x, 0.999), 3.0, 12.0)
    n_bins = cfg.p2_n_bins
    bin_width = upper / n_bins
    centers = [(i - 0.5) * bin_width for i in 1:n_bins]

    counts = zeros(Float64, n_bins)
    n_in_range = 0
    for v in x
        if 0 <= v < upper
            idx = clamp(floor(Int, v / bin_width) + 1, 1, n_bins)
            counts[idx] += 1.0
            n_in_range += 1
        end
    end
    n_in_range == 0 && return nothing

    # Method-of-moments init for the Gamma single-emitter model.
    μ = mean(x)
    σ2 = var(x)
    k0 = clamp(μ^2 / max(σ2, 1e-6), 0.5, 50.0)
    θ0 = clamp(σ2 / max(μ, 1e-6), 1e-3, 100.0)

    # Model density on the grid: (1-p)·f_single + p·f_single⊛f_single.
    model = function (k, θ, p)
        fs = [pdf(Gamma(k, θ), c) for c in centers]
        s = sum(fs)
        (s <= 0 || !isfinite(s)) && return nothing
        fs ./= s
        fd = _convolve_1d(fs)[1:n_bins]
        (1 - p) .* fs .+ p .* fd
    end

    # Multinomial negative log-likelihood of the binned counts.
    nll = function (par)
        k = exp(par[1]); θ = exp(par[2]); p = 1 / (1 + exp(-par[3]))
        (!isfinite(k) || !isfinite(θ) || k <= 0.1 || k > 200 || θ <= 0) && return 1e12
        m = model(k, θ, p)
        m === nothing && return 1e12
        s = 0.0
        @inbounds for i in 1:n_bins
            mi = m[i] < 1e-12 ? 1e-12 : m[i]
            s -= counts[i] * log(mi)
        end
        isfinite(s) ? s : 1e12
    end

    p0 = [log(k0), log(θ0), log(0.01 / 0.99)]
    res = try
        optimize(nll, p0, NelderMead(), Optim.Options(iterations=5000, g_tol=1e-9))
    catch err
        err isa InterruptException && rethrow()
        @warn "intensityfilter: p₂ mixture fit failed; reporting no estimate" exception=err
        return nothing
    end

    p = 1 / (1 + exp(-Optim.minimizer(res)[3]))
    isfinite(p) ? clamp(p, 0.0, 1.0) : nothing
end

"""
    _estimate_p2(xs, ys, photons, field_model, cfg) -> (p2, tail_obs, tail_f2)

Estimate the double-emitter fraction p₂ (info-only) from field-normalized photons.

`p2` is produced by `cfg.p2_method`:
- `:mixture` (default) — threshold-free mixture fit (`_estimate_p2_mixture`); unbiased.
- `:tail` — legacy tail-ratio `tail_obs / tail_f2` at `cfg.p2_tail_threshold`; biased high.

`tail_obs` / `tail_f2` are always returned as tail-mass diagnostics regardless of method.
Returns `(nothing, nothing, nothing)` if there are too few localizations to estimate.
"""
function _estimate_p2(xs, ys, photons, field_model, cfg::IntensityFilterConfig)
    n = length(xs)
    n < 100 && return (nothing, nothing, nothing)

    normalized = _field_normalized_photons(xs, ys, photons, field_model)

    tail_obs, tail_f2 = _p2_tail_masses(normalized, cfg)

    p2 = if cfg.p2_method === :mixture
        _estimate_p2_mixture(normalized, cfg)
    elseif cfg.p2_method === :tail
        (tail_obs === nothing || tail_f2 === nothing || tail_f2 <= 0) ? nothing :
            tail_obs / tail_f2
    else
        throw(ArgumentError("Unknown p2_method $(cfg.p2_method); expected :mixture or :tail"))
    end

    return (p2, tail_obs, tail_f2)
end

# ============================================================
# Diagnostic outputs
# ============================================================

function _save_intensityfilter_outputs!(dir, cfg, v, info, xs, ys, photons, pvals, keep, field_model)
    mkpath(dir)
    _save_config!(dir, cfg)
    _save_info!(dir, info)

    if v >= Verbosity.STANDARD
        _write_intensityfilter_stats(dir, cfg, info)
        _save_intensityfilter_figures(dir, xs, ys, photons, pvals, keep, field_model, cfg, info)
    end
end

function _write_intensityfilter_stats(dir, cfg, info)
    n_rejected = info.n_before - info.n_after
    rej_pct = info.n_before > 0 ? round(100 * n_rejected / info.n_before, digits=1) : 0.0

    filepath = joinpath(dir, "stats.md")
    open(filepath, "w") do io
        println(io, "# Intensity Filter Statistics\n")
        println(io, "## Summary")
        println(io, "- **Input**: $(info.n_before)")
        println(io, "- **Output**: $(info.n_after)")
        println(io, "- **Rejected**: $n_rejected ($rej_pct%)")
        println(io, "- **P-value cutoff**: $(cfg.cutoff)")
        println(io, "- **Rate percentile**: $(cfg.rate_percentile)")
        println(io, "- **Field mode**: $(info.field_mode)")
        println(io, "- **Global rate (p$(round(Int, cfg.rate_percentile*100)))**: $(round(info.lambda_max_global, digits=1)) photons")
        println(io, "- **Field fit R²**: $(round(info.field_fit_r2, digits=3))")
        println(io, "- **Time**: $(round(info.elapsed_s, digits=2))s")
        if info.field_params !== nothing
            p = info.field_params
            println(io, "")
            println(io, "## Gaussian Field Parameters")
            println(io, "- A: $(round(p.A, digits=1))")
            println(io, "- x₀: $(round(p.x0, digits=3)) μm")
            println(io, "- y₀: $(round(p.y0, digits=3)) μm")
            println(io, "- w: $(round(p.w, digits=3)) μm")
            println(io, "- bg: $(round(p.bg, digits=1))")
        end
        if info.p2_estimate !== nothing
            println(io, "")
            println(io, "## Double-Emitter Rate (p₂)")
            println(io, "- **p₂**: $(round(100 * info.p2_estimate, digits=2))%")
            println(io, "- **Tail threshold**: $(info.p2_tail_threshold)×λ")
            println(io, "- **Tail mass (observed)**: $(round(info.p2_tail_obs, sigdigits=3))")
            println(io, "- **Tail mass (double-conv)**: $(round(info.p2_tail_f2, sigdigits=3))")
            n_est_doubles = round(Int, info.p2_estimate * info.n_before)
            println(io, "- **Estimated doubles**: ~$n_est_doubles / $(info.n_before)")
        end
    end
end

function _save_intensityfilter_figures(dir, xs, ys, photons, pvals, keep, field_model, cfg, info)
    isempty(xs) && return

    # Recompute bin data for visualization
    bin_result = _spatial_bin_rates(xs, ys, photons, cfg)

    # --- Figure 1: Excitation field (scatter + binned rate + fitted field) ---
    fig1 = Figure(size=(1400, 500))

    x_centers = [(bin_result.x_edges[i] + bin_result.x_edges[i+1]) / 2 for i in 1:length(bin_result.x_edges)-1]
    y_centers = [(bin_result.y_edges[i] + bin_result.y_edges[i+1]) / 2 for i in 1:length(bin_result.y_edges)-1]

    # Left: scatter plot colored by photon count (subsampled for visibility)
    ax1a = Axis(fig1[1, 1], xlabel="x (μm)", ylabel="y (μm)", title="Photon Counts",
                aspect=DataAspect(), yreversed=true)
    p01 = quantile(photons, 0.01)
    p99_scatter = quantile(photons, 0.99)
    n_max_scatter = 10_000
    if length(xs) > n_max_scatter
        step = cld(length(xs), n_max_scatter)
        idx = 1:step:length(xs)
    else
        idx = eachindex(xs)
    end
    sc = scatter!(ax1a, xs[idx], ys[idx], color=photons[idx], colormap=:viridis,
                  colorrange=(p01, p99_scatter), markersize=2)
    Colorbar(fig1[1, 2], sc, label="Photons")

    # Middle: binned rate (percentile) — Makie convention: grid[ix, iy] → (x[ix], y[iy])
    valid_rates = filter(!isnan, vec(bin_result.rate_grid))
    ax1b = Axis(fig1[1, 3], xlabel="x (μm)", ylabel="y (μm)",
                title="p$(round(Int, cfg.rate_percentile*100)) Photons (λ estimate)",
                aspect=DataAspect(), yreversed=true)
    if !isempty(valid_rates)
        hm1 = heatmap!(ax1b, x_centers, y_centers, bin_result.rate_grid, colormap=:viridis,
                        colorrange=(minimum(valid_rates), maximum(valid_rates)), nan_color=:gray20)
        Colorbar(fig1[1, 4], hm1, label="λ (photons)")
    end

    # Right: fitted field — comprehension [f(x,y) for x in X, y in Y] → grid[ix, iy]
    ax1c = Axis(fig1[1, 5], xlabel="x (μm)", ylabel="y (μm)", title="Fitted Field λ(x,y)",
                aspect=DataAspect(), yreversed=true)
    x_range = range(minimum(xs), maximum(xs), length=50)
    y_range = range(minimum(ys), maximum(ys), length=50)
    field_grid = [field_model(x, y) for x in x_range, y in y_range]
    hm2 = heatmap!(ax1c, collect(x_range), collect(y_range), field_grid, colormap=:viridis)
    Colorbar(fig1[1, 6], hm2, label="λ (photons)")

    save(joinpath(dir, "excitation_field.png"), fig1)

    # --- Figure 2: P-value histogram ---
    fig2 = Figure(size=(700, 450))
    ax2 = Axis(fig2[1, 1], xlabel="log₁₀(p-value)", ylabel="Count",
               title="P-value Distribution (cutoff=$(cfg.cutoff))")
    pvals_pos = pvals[pvals .> 0]
    if !isempty(pvals_pos)
        log_pvals = log10.(pvals_pos)
        # Clip extreme values for readable histogram
        log_min = max(minimum(log_pvals), -50.0)
        log_pvals_clipped = clamp.(log_pvals, log_min, 0.0)
        hist!(ax2, log_pvals_clipped, bins=80, color=:steelblue)
        vlines!(ax2, [log10(cfg.cutoff)], color=:red, linewidth=2, linestyle=:dash,
                label="cutoff ($(cfg.cutoff))")
        n_rej = sum(pvals .< cfg.cutoff)
        n_tot = length(pvals)
        n_clipped = sum(log_pvals .< log_min)
        clip_text = n_clipped > 0 ? "\n($(n_clipped) clipped below $log_min)" : ""
        text!(ax2, 0.03, 0.95,
              text="Rejected: $n_rej / $n_tot ($(round(100*n_rej/n_tot, digits=1))%)$clip_text",
              align=(:left, :top), space=:relative, fontsize=12)
        axislegend(ax2, position=:lt, framevisible=false, offset=(0, -40))
    end
    save(joinpath(dir, "pvalue_histogram.png"), fig2)

    # --- Figure 3: Rejection map ---
    fig3 = Figure(size=(700, 550))
    ax3 = Axis(fig3[1, 1], xlabel="x (μm)", ylabel="y (μm)", title="Rejection Map",
               aspect=DataAspect(), yreversed=true)
    rejected = .!keep
    if any(keep)
        scatter!(ax3, xs[keep], ys[keep], color=:steelblue, markersize=1,
                 label="kept ($(sum(keep)))")
    end
    if any(rejected)
        scatter!(ax3, xs[rejected], ys[rejected], color=:red, markersize=2,
                 label="rejected ($(sum(rejected)))")
    end
    axislegend(ax3, position=:rb, framevisible=false, labelsize=10)
    save(joinpath(dir, "rejection_map.png"), fig3)

    # --- Figure 4: Photon distribution (raw + field-normalized) ---
    fig4 = Figure(size=(1400, 450))

    # Left: raw photon distribution
    p99 = quantile(photons, 0.99)
    ax4a = Axis(fig4[1, 1], xlabel="Photons", ylabel="Count",
                title="Raw Photon Distribution")
    hist!(ax4a, photons[photons .<= p99], bins=80, color=(:steelblue, 0.7))
    vlines!(ax4a, [info.lambda_max_global], color=:orange, linewidth=2, linestyle=:solid,
            label="λ global (p$(round(Int, cfg.rate_percentile*100))=$(round(Int, info.lambda_max_global)))")
    mode_val = _calculate_mode(photons)
    vlines!(ax4a, [mode_val], color=:green, linewidth=2, linestyle=:dash,
            label="mode=$(round(Int, mode_val))")
    λ_g = info.lambda_max_global
    if λ_g > 1
        n_thresh = round(Int, λ_g)
        while ccdf(Poisson(λ_g), n_thresh - 1) > cfg.cutoff && n_thresh < 10 * λ_g
            n_thresh += 1
        end
        if n_thresh < 10 * λ_g
            vlines!(ax4a, [n_thresh], color=:red, linewidth=2, linestyle=:dash,
                    label="reject above ≈$(n_thresh) (p<$(cfg.cutoff))")
        end
    end
    xlims!(ax4a, 0, p99)
    axislegend(ax4a, position=:rt, framevisible=false)

    # Right: field-normalized photon distribution (photons / λ(x,y))
    # Removes spatial excitation variation — remaining tail = genuine multi-emitter
    λ_local = [field_model(xs[i], ys[i]) for i in eachindex(xs)]
    normalized = photons ./ λ_local
    ax4b = Axis(fig4[1, 2], xlabel="Photons / λ(x,y)", ylabel="Count",
                title="Field-Normalized Distribution")
    norm_p99 = quantile(normalized, 0.99)
    hist!(ax4b, normalized[normalized .<= norm_p99], bins=80, color=(:steelblue, 0.7))
    vlines!(ax4b, [1.0], color=:orange, linewidth=2, linestyle=:solid,
            label="λ(x,y) (expected)")
    norm_mode = _calculate_mode(normalized)
    vlines!(ax4b, [norm_mode], color=:green, linewidth=2, linestyle=:dash,
            label="mode=$(round(norm_mode, digits=2))")
    # Rejection threshold in normalized units
    if λ_g > 1 && n_thresh < 10 * λ_g
        vlines!(ax4b, [n_thresh / λ_g], color=:red, linewidth=2, linestyle=:dash,
                label="reject ≈$(round(n_thresh / λ_g, digits=2))×λ")
    end
    xlims!(ax4b, 0, norm_p99)
    axislegend(ax4b, position=:rt, framevisible=false)

    save(joinpath(dir, "photon_distribution.png"), fig4)

    # --- Figure 5: Mixture decomposition (p₂ estimate) ---
    if info.p2_estimate !== nothing
        _save_mixture_decomposition_figure(dir, xs, ys, photons, field_model, cfg, info)
    end
end

function _save_mixture_decomposition_figure(dir, xs, ys, photons, field_model, cfg, info)
    n = length(xs)

    # Recompute field-normalized photons
    normalized = Vector{Float64}(undef, n)
    for i in 1:n
        λ = field_model(xs[i], ys[i])
        normalized[i] = photons[i] / max(1.0, λ)
    end

    upper = quantile(normalized, 0.995)
    upper <= 0 && return
    n_bins = cfg.p2_n_bins
    bin_width = upper / n_bins
    τ = cfg.p2_tail_threshold

    # Build histogram
    h = zeros(Float64, n_bins)
    n_in_range = 0
    for v in normalized
        if 0 <= v < upper
            idx = clamp(floor(Int, v / bin_width) + 1, 1, n_bins)
            h[idx] += 1.0
            n_in_range += 1
        end
    end
    n_in_range == 0 && return
    h_prob = h ./ n_in_range

    # Self-convolve
    f2 = _convolve_1d(h_prob)

    # Bin centers for f_obs and f₂
    centers_obs = [(i - 0.5) * bin_width for i in 1:n_bins]
    centers_f2 = [(i - 1) * bin_width for i in 1:length(f2)]

    p2 = info.p2_estimate

    fig5 = Figure(size=(1400, 500))

    # Left panel: full distribution with overlay
    ax5a = Axis(fig5[1, 1], xlabel="Photons / λ(x,y)", ylabel="Probability",
                title="Mixture Decomposition")
    barplot!(ax5a, centers_obs, h_prob, width=bin_width * 0.9,
             color=(:steelblue, 0.5), label="f_obs")
    # Scale f₂ by p₂
    lines!(ax5a, centers_f2, p2 .* f2, color=:red, linewidth=2,
           label="p₂ × f₂ (doubles)")
    # Scale f₁ by (1-p₂) — f₁ ≈ f_obs for p₂ << 1
    lines!(ax5a, centers_obs, (1 - p2) .* h_prob, color=:blue, linewidth=1.5,
           linestyle=:dash, label="(1-p₂) × f₁ (singles)")
    vlines!(ax5a, [τ], color=:orange, linewidth=2, linestyle=:dash,
            label="τ = $(τ)")
    xlims!(ax5a, 0, min(upper, 3.0))
    axislegend(ax5a, position=:rt, framevisible=false)
    text!(ax5a, 0.03, 0.85,
          text="p₂ = $(round(100 * p2, digits=2))%\nn_total = $(info.n_before)",
          align=(:left, :top), space=:relative, fontsize=13)

    # Right panel: tail zoom
    ax5b = Axis(fig5[1, 2], xlabel="Photons / λ(x,y)", ylabel="Probability",
                title="Tail Detail (> 0.5×λ)")
    # Only plot tail region
    tail_mask_obs = centers_obs .>= 0.5
    tail_mask_f2 = centers_f2 .>= 0.5
    if any(tail_mask_obs)
        barplot!(ax5b, centers_obs[tail_mask_obs], h_prob[tail_mask_obs],
                 width=bin_width * 0.9, color=(:steelblue, 0.5), label="f_obs (tail)")
    end
    if any(tail_mask_f2)
        lines!(ax5b, centers_f2[tail_mask_f2], p2 .* f2[tail_mask_f2],
               color=:red, linewidth=2, label="p₂ × f₂")
    end
    vlines!(ax5b, [τ], color=:orange, linewidth=2, linestyle=:dash, label="τ = $(τ)")

    # Count tail events
    τ_bin = clamp(floor(Int, τ / bin_width) + 1, 1, n_bins)
    n_above_τ = sum(v >= τ for v in normalized)
    n_est_doubles = round(Int, p2 * info.n_before)
    text!(ax5b, 0.97, 0.60,
          text="Above τ: $n_above_τ\nEst. doubles: ~$n_est_doubles",
          align=(:right, :top), space=:relative, fontsize=13)
    axislegend(ax5b, position=:lt, framevisible=false)

    save(joinpath(dir, "mixture_decomposition.png"), fig5)
end
