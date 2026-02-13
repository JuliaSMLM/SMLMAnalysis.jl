"""
Filtering step - filters localizations by various criteria
"""

"""
    FilterConfig <: AbstractSMLMConfig

Quality-based filtering of localizations. All criteria use `(min, max)` tuples.

# Keywords
- `photons`: Photon count range, e.g. `(500.0, Inf)`
- `precision`: Localization precision range in microns, e.g. `(0.0, 0.007)`
- `pvalue`: Goodness-of-fit p-value range, e.g. `(1e-3, 1.0)`
- `psf_sigma`: PSF width filter. `:auto` uses mode ± 10%, or explicit `(min, max)` in microns

All filters default to `nothing` (disabled).
"""
@kwdef struct FilterConfig <: SMLMData.AbstractSMLMConfig
    # All filters use (min, max) tuples. Use -Inf/Inf for unbounded.
    photons::Union{Tuple{Float64, Float64}, Nothing} = nothing      # (min, max)
    precision::Union{Tuple{Float64, Float64}, Nothing} = nothing    # (min, max) in microns
    pvalue::Union{Tuple{Float64, Float64}, Nothing} = nothing       # (min, max)
    # PSF sigma: :auto (mode ± 10%), or (min, max) tuple in microns
    psf_sigma::Union{Symbol, Tuple{Float64, Float64}, Nothing} = nothing
end

"""
    filter_step(smld, cfg; smld_raw=nothing, outdir=nothing, step_number=0, verbose=Verbosity.STANDARD)

Filter localizations by quality criteria. Returns `(filtered_smld, info)`.

# Arguments
- `smld::BasicSMLD`: Input localizations
- `cfg::FilterConfig`: Filter criteria

# Keyword Arguments
- `smld_raw`: Original unfiltered SMLD for detailed output diagnostics
- `outdir`: Output directory (nothing to skip file output)
- `step_number`: Step number for output directory naming
- `verbose`: Verbosity level

# Returns
`(filtered_smld, (step_record, n_before, n_after))`
"""
function filter_step(smld::BasicSMLD, cfg::FilterConfig;
                     smld_raw::Union{BasicSMLD,Nothing}=nothing,
                     outdir::Union{String,Nothing}=nothing,
                     step_number::Int=0,
                     verbose::Int=Verbosity.STANDARD)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)

    v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg))" photons=cfg.photons precision=cfg.precision

    n_before = length(smld.emitters)
    t = @elapsed filtered = _filter_smld(smld, cfg)
    n_after = length(filtered.emitters)

    summary = Dict{Symbol,Any}(
        :n_before => n_before,
        :n_after => n_after,
        :acceptance => round(n_after / n_before, digits=3)
    )
    record = StepRecord(step_number, cfg, t, summary)

    if dir !== nothing
        _save_filter_outputs!(dir, cfg, v, t, n_before, n_after, smld_raw, filtered)
    end

    v >= Verbosity.PROGRESS && @info "  → $n_after / $n_before ($(round(t, digits=2))s)"
    (filtered, (step_record=record, n_before=n_before, n_after=n_after))
end

"""
    analyze(smld, cfg::FilterConfig; kwargs...) -> (filtered_smld, info)

Filter localizations by quality criteria.
"""
analyze(smld::BasicSMLD, cfg::FilterConfig; kwargs...) = filter_step(smld, cfg; kwargs...)

function _filter_smld(smld::BasicSMLD, cfg::FilterConfig)
    emitters = smld.emitters
    mask = trues(length(emitters))

    if cfg.photons !== nothing
        lo, hi = cfg.photons
        mask .&= [lo <= e.photons <= hi for e in emitters]
    end

    if cfg.precision !== nothing
        lo, hi = cfg.precision
        mask .&= [lo <= max(e.σ_x, e.σ_y) <= hi for e in emitters]
    end

    if cfg.pvalue !== nothing
        lo, hi = cfg.pvalue
        mask .&= [lo <= e.pvalue <= hi for e in emitters]
    end

    if cfg.psf_sigma !== nothing && length(emitters) > 0
        # Determine bounds: :auto calculates mode ± 10%, or use explicit (min, max)
        if hasproperty(emitters[1], :σ)
            lo, hi = _get_psf_sigma_bounds(cfg.psf_sigma, [e.σ for e in emitters])
            if lo > 0 && hi > 0
                mask .&= [lo <= e.σ <= hi for e in emitters]
            end
        elseif hasproperty(emitters[1], :σx) && hasproperty(emitters[1], :σy)
            lo_x, hi_x = _get_psf_sigma_bounds(cfg.psf_sigma, [e.σx for e in emitters])
            lo_y, hi_y = _get_psf_sigma_bounds(cfg.psf_sigma, [e.σy for e in emitters])
            if lo_x > 0 && hi_x > 0 && lo_y > 0 && hi_y > 0
                mask .&= [lo_x <= e.σx <= hi_x && lo_y <= e.σy <= hi_y for e in emitters]
            end
        end
    end

    filtered = emitters[mask]
    BasicSMLD(filtered, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)
end

"""
    _get_psf_sigma_bounds(range_spec, values) -> (lo, hi)

Calculate PSF sigma filter bounds.
- `:auto` → mode ± 10%
- `(min, max)` → explicit bounds in microns
"""
function _get_psf_sigma_bounds(range_spec, values::Vector)
    if range_spec === :auto
        mode = _calculate_mode(values)
        mode > 0 || return (0.0, 0.0)
        return (mode * 0.90, mode * 1.10)
    elseif range_spec isa Tuple{Float64, Float64}
        return range_spec
    else
        error("psf_sigma_range must be :auto or (min, max) tuple, got: $range_spec")
    end
end

function _save_filter_outputs!(dir::String, cfg::FilterConfig, v::Int, t::Float64,
                               n_before::Int, n_after::Int,
                               smld_raw::Union{BasicSMLD,Nothing}, smld_filtered::BasicSMLD)
    mkpath(dir)
    _save_config!(dir, cfg)

    if v >= Verbosity.STANDARD
        _write_filter_stats(dir, cfg, n_before, n_after, t)

        # Fit quality figures showing distributions with filter thresholds
        if smld_raw !== nothing
            _save_filter_quality_figures(dir, smld_raw, cfg)
        end
    end

    if v >= Verbosity.DETAILED && smld_raw !== nothing
        _save_filter_detailed(dir, smld_raw, smld_filtered, cfg)
    end
end

function _write_filter_stats(dir, cfg, n_before, n_after, t)
    acceptance = n_after / n_before

    filepath = joinpath(dir, "stats.md")
    open(filepath, "w") do io
        println(io, "# Filter Statistics\n")
        println(io, "## Summary")
        println(io, "- **Input**: $n_before")
        println(io, "- **Output**: $n_after")
        println(io, "- **Acceptance**: $(round(100*acceptance, digits=1))%")
        println(io, "- **Time**: $(round(t*1000, digits=1))ms")
        println(io, "")
        println(io, "## Criteria Applied")
        if cfg.photons !== nothing
            lo, hi = cfg.photons
            println(io, "- photons: $(lo) - $(hi == Inf ? "∞" : hi)")
        end
        if cfg.precision !== nothing
            lo, hi = cfg.precision
            println(io, "- precision: $(round(lo*1000, digits=1)) - $(hi == Inf ? "∞" : round(hi*1000, digits=1)) nm")
        end
        if cfg.pvalue !== nothing
            lo, hi = cfg.pvalue
            println(io, "- pvalue: $(lo) - $(hi)")
        end
        if cfg.psf_sigma !== nothing
            if cfg.psf_sigma === :auto
                println(io, "- psf_sigma: :auto (mode ± 10%)")
            else
                lo, hi = cfg.psf_sigma
                println(io, "- psf_sigma: $(round(lo*1000, digits=1)) - $(round(hi*1000, digits=1)) nm")
            end
        end
    end
end

function _save_filter_detailed(dir, smld_raw, smld_filtered, cfg)
    # Show which criteria rejected what
    emitters = smld_raw.emitters
    n = length(emitters)

    filepath = joinpath(dir, "detailed_stats.md")
    open(filepath, "w") do io
        println(io, "# Filter Breakdown\n")
        println(io, "| Criterion | Pass | Fail | % Pass |")
        println(io, "|-----------|------|------|--------|")

        if cfg.photons !== nothing
            lo, hi = cfg.photons
            pass = sum(lo <= e.photons <= hi for e in emitters)
            hi_str = hi == Inf ? "∞" : string(hi)
            println(io, "| Photons ∈ [$lo, $hi_str] | $pass | $(n - pass) | $(round(100*pass/n, digits=1))% |")
        end

        if cfg.precision !== nothing
            lo, hi = cfg.precision
            pass = sum(lo <= max(e.σ_x, e.σ_y) <= hi for e in emitters)
            hi_str = hi == Inf ? "∞" : "$(round(hi*1000, digits=1))nm"
            println(io, "| Precision ∈ [$(round(lo*1000, digits=1))nm, $hi_str] | $pass | $(n - pass) | $(round(100*pass/n, digits=1))% |")
        end

        if cfg.pvalue !== nothing
            lo, hi = cfg.pvalue
            pass = sum(lo <= e.pvalue <= hi for e in emitters)
            println(io, "| P-value ∈ [$lo, $hi] | $pass | $(n - pass) | $(round(100*pass/n, digits=1))% |")
        end
    end
end

"""
    _save_filter_quality_figures(dir, smld_raw, cfg)

Generate fit quality distribution plots with filter thresholds highlighted.
Shows photons, background, precision, p-value, and PSF sigma distributions
from the raw (unfiltered) SMLD, with rejected regions marked per FilterConfig.
"""
function _save_filter_quality_figures(dir, smld_raw, cfg::FilterConfig)
    emitters = smld_raw.emitters
    isempty(emitters) && return

    photons = [e.photons for e in emitters]
    bg = [e.bg for e in emitters]
    σ_x = [e.σ_x for e in emitters] .* 1000  # precision in nm
    σ_y = [e.σ_y for e in emitters] .* 1000
    pvalue = [e.pvalue for e in emitters]

    # Check PSF model type
    has_psf_iso = hasproperty(emitters[1], :σ)
    has_psf_aniso = hasproperty(emitters[1], :σx)

    # Colors for consistent styling
    REJECTED_COLOR = (:gray30, 0.5)
    MEAN_COLOR = :blue
    MEDIAN_COLOR = :red
    THRESHOLD_COLOR = :black

    fig = Figure(size=(1200, 900))

    # Row 1: Photons and Background
    p98 = quantile(photons, 0.98)
    ax1 = Axis(fig[1, 1], xlabel="Photons", ylabel="Count", title="Photon Distribution")
    if cfg.photons !== nothing
        photons_lo = cfg.photons[1]
        vspan!(ax1, 0, photons_lo, color=REJECTED_COLOR)
        vlines!(ax1, [photons_lo], color=THRESHOLD_COLOR, linestyle=:dot, linewidth=2)
    end
    hist!(ax1, photons[photons .<= p98], bins=50)
    vlines!(ax1, [mean(photons)], color=MEAN_COLOR, linestyle=:solid, linewidth=2)
    vlines!(ax1, [median(photons)], color=MEDIAN_COLOR, linestyle=:dash, linewidth=2)
    xlims!(ax1, 0, p98)
    text!(ax1, 0.97, 0.95, text="mean: $(round(Int, mean(photons)))\nmedian: $(round(Int, median(photons)))",
          align=(:right, :top), space=:relative, fontsize=10)

    bg98 = quantile(bg, 0.98)
    ax2 = Axis(fig[1, 2], xlabel="Background", ylabel="Count", title="Background Distribution")
    hist!(ax2, bg[bg .<= bg98], bins=50)
    vlines!(ax2, [mean(bg)], color=MEAN_COLOR, linestyle=:solid, linewidth=2)
    vlines!(ax2, [median(bg)], color=MEDIAN_COLOR, linestyle=:dash, linewidth=2)
    xlims!(ax2, 0, bg98)
    text!(ax2, 0.97, 0.95, text="mean: $(round(mean(bg), digits=1))\nmedian: $(round(median(bg), digits=1))",
          align=(:right, :top), space=:relative, fontsize=10)

    # Row 2: Precision and P-value
    prec_data = vcat(σ_x, σ_y)
    prec98 = quantile(prec_data, 0.98)
    prec_xlim = prec98

    ax3 = Axis(fig[2, 1], xlabel="Localization Precision (nm)", ylabel="Count", title="Precision Distribution")
    if cfg.precision !== nothing
        prec_hi = cfg.precision[2] * 1000  # convert um to nm
        prec_xlim = max(prec98, prec_hi * 1.5)
        vspan!(ax3, prec_hi, prec_xlim, color=REJECTED_COLOR)
        vlines!(ax3, [prec_hi], color=THRESHOLD_COLOR, linestyle=:dot, linewidth=2)
    end
    hist!(ax3, σ_x[σ_x .<= prec98], bins=50, color=(:blue, 0.5), label="σ_x")
    hist!(ax3, σ_y[σ_y .<= prec98], bins=50, color=(:red, 0.5), label="σ_y")
    xlims!(ax3, 0, prec_xlim)
    axislegend(ax3, position=:rt, framevisible=false, labelsize=9)
    text!(ax3, 0.97, 0.70, text="σ_x: $(round(median(σ_x), digits=1)) nm\nσ_y: $(round(median(σ_y), digits=1)) nm",
          align=(:right, :top), space=:relative, fontsize=10)

    ax4 = Axis(fig[2, 2], xlabel="log10(p-value)", ylabel="Density", title="P-value Distribution")
    pval_nonzero = pvalue[pvalue .> 0]
    pval_thresh = cfg.pvalue !== nothing ? cfg.pvalue[1] : nothing

    if !isempty(pval_nonzero)
        log_pval = log10.(pval_nonzero)
        pval_lo = quantile(log_pval, 0.02)
        log_pval_filtered = log_pval[log_pval .>= pval_lo]
        if pval_thresh !== nothing
            vspan!(ax4, pval_lo - 1, log10(pval_thresh), color=REJECTED_COLOR)
        end
        hist!(ax4, log_pval_filtered, bins=50, normalization=:pdf, color=(:steelblue, 0.7))
        # Compute histogram max for y limits (based on data, not theory)
        nbins = 50
        bin_edges = range(pval_lo, 0, length=nbins+1)
        bin_width = step(bin_edges)
        counts = zeros(Int, nbins)
        for v in log_pval_filtered
            idx = clamp(floor(Int, (v - pval_lo) / bin_width) + 1, 1, nbins)
            counts[idx] += 1
        end
        max_density = maximum(counts) / (length(log_pval_filtered) * bin_width)
        # Theory curve (uniform p-values -> exponential in log space)
        u_range = range(0, -pval_lo, length=100)
        theory_pdf = log(10) .* (10.0 .^ (-u_range))
        lines!(ax4, -u_range, theory_pdf, color=:red, linewidth=2, label="Uniform theory")
        vlines!(ax4, [mean(log_pval)], color=MEAN_COLOR, linestyle=:solid, linewidth=1.5)
        vlines!(ax4, [median(log_pval)], color=MEDIAN_COLOR, linestyle=:dash, linewidth=1.5)
        if pval_thresh !== nothing
            vlines!(ax4, [log10(pval_thresh)], color=THRESHOLD_COLOR, linestyle=:dot, linewidth=2)
            pval_pass_pct = round(100 * sum(pvalue .> pval_thresh) / length(pvalue), digits=1)
            text!(ax4, 0.03, 0.95, text="pass: $(pval_pass_pct)%\nthreshold: $(pval_thresh)",
                  align=(:left, :top), space=:relative, fontsize=10)
        end
        xlims!(ax4, pval_lo, 0)
        ylims!(ax4, 0, max_density * 1.1)
    end

    # Row 3: PSF Sigma
    if has_psf_aniso
        psf_σx = [e.σx for e in emitters] .* 1000
        psf_σy = [e.σy for e in emitters] .* 1000
        psf_data = vcat(psf_σx, psf_σy)
        psf98 = quantile(psf_data, 0.98)
        psf02 = quantile(psf_data, 0.02)
        mode_x = _calculate_mode([e.σx for e in emitters]) * 1000
        mode_y = _calculate_mode([e.σy for e in emitters]) * 1000
        mode_avg = (mode_x + mode_y) / 2

        ax5 = Axis(fig[3, 1], xlabel="Fitted PSF σ (nm)", ylabel="Count", title="PSF Width Distribution")
        # Show PSF sigma filter bounds if configured
        if cfg.psf_sigma !== nothing
            if cfg.psf_sigma === :auto
                tol = 0.10
                vspan!(ax5, psf02 * 0.9, mode_avg * (1 - tol), color=REJECTED_COLOR)
                vspan!(ax5, mode_avg * (1 + tol), psf98 * 1.1, color=REJECTED_COLOR)
                vlines!(ax5, [mode_avg * (1 - tol), mode_avg * (1 + tol)], color=THRESHOLD_COLOR, linestyle=:dot, linewidth=2)
            elseif cfg.psf_sigma isa Tuple
                lo_nm, hi_nm = cfg.psf_sigma[1] * 1000, cfg.psf_sigma[2] * 1000
                vspan!(ax5, psf02 * 0.9, lo_nm, color=REJECTED_COLOR)
                vspan!(ax5, hi_nm, psf98 * 1.1, color=REJECTED_COLOR)
                vlines!(ax5, [lo_nm, hi_nm], color=THRESHOLD_COLOR, linestyle=:dot, linewidth=2)
            end
        end
        hist!(ax5, psf_σx[(psf_σx .>= psf02) .& (psf_σx .<= psf98)], bins=50, color=(:blue, 0.5), label="σx")
        hist!(ax5, psf_σy[(psf_σy .>= psf02) .& (psf_σy .<= psf98)], bins=50, color=(:red, 0.5), label="σy")
        xlims!(ax5, psf02, psf98)
        axislegend(ax5, position=:rt, framevisible=false, labelsize=9)
        text!(ax5, 0.97, 0.70, text="mode σx: $(round(mode_x, digits=1)) nm\nmode σy: $(round(mode_y, digits=1)) nm",
              align=(:right, :top), space=:relative, fontsize=10)

    elseif has_psf_iso
        psf_σ = [e.σ for e in emitters] .* 1000
        psf98 = quantile(psf_σ, 0.98)
        psf02 = quantile(psf_σ, 0.02)
        mode_σ = _calculate_mode([e.σ for e in emitters]) * 1000

        ax5 = Axis(fig[3, 1], xlabel="Fitted PSF σ (nm)", ylabel="Count", title="PSF Width Distribution")
        if cfg.psf_sigma !== nothing
            if cfg.psf_sigma === :auto
                tol = 0.10
                lo_bound = mode_σ * (1 - tol)
                hi_bound = mode_σ * (1 + tol)
            else
                lo_bound = cfg.psf_sigma[1] * 1000
                hi_bound = cfg.psf_sigma[2] * 1000
            end
            psf_xmin = min(psf02, lo_bound * 0.95)
            psf_xmax = max(psf98, hi_bound * 1.05)
            vspan!(ax5, psf_xmin, lo_bound, color=REJECTED_COLOR)
            vspan!(ax5, hi_bound, psf_xmax, color=REJECTED_COLOR)
            vlines!(ax5, [lo_bound, hi_bound], color=THRESHOLD_COLOR, linestyle=:dot, linewidth=2)
        else
            psf_xmin = psf02
            psf_xmax = psf98
        end
        hist!(ax5, psf_σ[(psf_σ .>= psf02) .& (psf_σ .<= psf98)], bins=50)
        vlines!(ax5, [mean(psf_σ)], color=MEAN_COLOR, linestyle=:solid, linewidth=2)
        vlines!(ax5, [median(psf_σ)], color=MEDIAN_COLOR, linestyle=:dash, linewidth=2)
        xlims!(ax5, psf_xmin, psf_xmax)
        text!(ax5, 0.97, 0.95, text="mode: $(round(mode_σ, digits=1)) nm\nmean: $(round(mean(psf_σ), digits=1)) nm",
              align=(:right, :top), space=:relative, fontsize=10)
    else
        ax5 = Axis(fig[3, 1], xlabel="PSF σ (nm)", ylabel="", title="PSF Width (Fixed)")
        text!(ax5, 0.5, 0.5, text="Fixed PSF model",
              align=(:center, :center), space=:relative, fontsize=14)
        hideydecorations!(ax5)
    end

    # Row 3, Col 2: Legend
    ax6 = Axis(fig[3, 2], title="Legend")
    hidedecorations!(ax6)
    hidespines!(ax6)
    lines!(ax6, [0.1, 0.25], [0.8, 0.8], color=MEAN_COLOR, linewidth=2)
    text!(ax6, 0.3, 0.8, text="Mean", fontsize=12)
    lines!(ax6, [0.1, 0.25], [0.6, 0.6], color=MEDIAN_COLOR, linewidth=2, linestyle=:dash)
    text!(ax6, 0.3, 0.6, text="Median", fontsize=12)
    lines!(ax6, [0.1, 0.25], [0.4, 0.4], color=THRESHOLD_COLOR, linewidth=2, linestyle=:dot)
    text!(ax6, 0.3, 0.4, text="Filter Threshold", fontsize=12)
    poly!(ax6, Point2f[(0.1, 0.15), (0.25, 0.15), (0.25, 0.25), (0.1, 0.25)], color=REJECTED_COLOR)
    text!(ax6, 0.3, 0.2, text="Rejected Region", fontsize=12)
    xlims!(ax6, 0, 1)
    ylims!(ax6, 0, 1)

    save(joinpath(dir, "fit_quality.png"), fig)
end
