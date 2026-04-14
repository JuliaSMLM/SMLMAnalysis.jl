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
    filter_step(smld, cfg; outdir=nothing, step_number=0, verbose=Verbosity.STANDARD)

Filter localizations by quality criteria. Returns `(filtered_smld, FilterInfo)`.

Diagnostic plots use `smld` (the input) for pre-filter distributions.

# Arguments
- `smld::BasicSMLD`: Input localizations
- `cfg::FilterConfig`: Filter criteria

# Keyword Arguments
- `outdir`: Output directory (nothing to skip file output)
- `step_number`: Step number for output directory naming
- `verbose`: Verbosity level

# Returns
`(filtered_smld, FilterInfo)`
"""
function filter_step(smld::BasicSMLD, cfg::FilterConfig;
                     outdir::Union{String,Nothing}=nothing,
                     step_number::Int=0,
                     verbose::Int=Verbosity.STANDARD)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)

    v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg))" photons=cfg.photons precision=cfg.precision

    n_before = length(smld.emitters)
    t = @elapsed filtered = _filter_smld(smld, cfg)
    n_after = length(filtered.emitters)

    if dir !== nothing
        _save_filter_outputs!(dir, outdir, cfg, v, t, n_before, n_after, smld, filtered)
    end

    v >= Verbosity.PROGRESS && @info "  → $n_after / $n_before ($(round(t, digits=2))s)"
    (filtered, FilterInfo(n_before, n_after, t))
end

_step_summary(info::FilterInfo) = Dict{Symbol,Any}(
    :n_before => info.n_before,
    :n_after => info.n_after,
    :acceptance => round(info.n_after / max(1, info.n_before), digits=3)
)

"""
    analyze(smld, cfg::FilterConfig; kwargs...) -> (filtered_smld, StepInfo)

Filter localizations by quality criteria.
"""
function analyze(smld::BasicSMLD, cfg::FilterConfig;
                 outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD,
                 checkpoint::Int=Checkpoint.EXPENSIVE, kwargs...)
    t = @elapsed (filtered, filter_info) = filter_step(smld, cfg;
        outdir=outdir, step_number=step_number, verbose=verbose)

    if checkpoint >= Checkpoint.ALL
        dir = step_outdir(outdir, step_number, cfg)
        _save_step_smld(dir, filtered; filename="smld_filtered.jld2")
    end

    (filtered, StepInfo(step_number, cfg, t, _step_summary(filter_info); info=filter_info))
end

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

function _save_filter_outputs!(dir::String, outdir::Union{String,Nothing}, cfg::FilterConfig, v::Int, t::Float64,
                               n_before::Int, n_after::Int,
                               smld_input::BasicSMLD, smld_filtered::BasicSMLD)
    mkpath(dir)
    _save_config!(dir, cfg)

    if v >= Verbosity.STANDARD
        _write_filter_stats(dir, cfg, n_before, n_after, t)
        _save_filter_quality_figures(dir, smld_input, cfg)
        _save_fit_overlay_from_cache(dir, outdir, smld_input, cfg)
        _save_loc_per_frame(dir, smld_filtered; title="Localizations per Frame (post-filter)")
    end

    if v >= Verbosity.DETAILED
        _save_filter_detailed(dir, smld_input, smld_filtered, cfg)
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
        prec_xlim = max(prec98, prec_hi * 1.2)
        vspan!(ax3, prec_hi, prec_xlim, color=REJECTED_COLOR)
        vlines!(ax3, [prec_hi], color=THRESHOLD_COLOR, linestyle=:dot, linewidth=2)
    end
    hist!(ax3, σ_x[σ_x .<= prec_xlim], bins=50, color=(:blue, 0.5), label="σ_x")
    hist!(ax3, σ_y[σ_y .<= prec_xlim], bins=50, color=(:red, 0.5), label="σ_y")
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

        # Compute per-axis bounds matching _filter_smld logic
        if cfg.psf_sigma !== nothing
            if cfg.psf_sigma === :auto
                tol = 0.10
                lo_x_nm, hi_x_nm = mode_x * (1 - tol), mode_x * (1 + tol)
                lo_y_nm, hi_y_nm = mode_y * (1 - tol), mode_y * (1 + tol)
            else
                lo_x_nm, hi_x_nm = cfg.psf_sigma[1] * 1000, cfg.psf_sigma[2] * 1000
                lo_y_nm, hi_y_nm = lo_x_nm, hi_x_nm
            end
        end

        ax5 = Axis(fig[3, 1], xlabel="Fitted PSF σ (nm)", ylabel="Count", title="PSF Width Distribution")
        # Show per-axis filter bounds (matching actual filter behavior)
        if cfg.psf_sigma !== nothing
            vlines!(ax5, [lo_x_nm, hi_x_nm], color=:blue, linestyle=:dot, linewidth=2, label="σx bounds")
            vlines!(ax5, [lo_y_nm, hi_y_nm], color=:red, linestyle=:dot, linewidth=2, label="σy bounds")
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

# ============================================================
# Fit overlay from detectfit cache (colored by FilterConfig thresholds)
# ============================================================

"""
Load detectfit sample cache and generate fit_overlay.png with FilterConfig thresholds.
Gracefully skips if cache is missing (e.g., standalone filter without prior detectfit).
"""
function _save_fit_overlay_from_cache(dir, outdir, smld_raw, cfg::FilterConfig)
    cache = load_cache(outdir, "detectfit_samples.jld2")
    cache === nothing && return

    # Reconstruct ROIBatch from cached components using smld_raw's camera
    roi_batch = ROIBatch(
        cache["sample_roi_data"],
        cache["sample_roi_x"],
        cache["sample_roi_y"],
        cache["sample_roi_frames"],
        smld_raw.camera,
    )

    _save_fit_overlay(dir, smld_raw, roi_batch,
        cache["sample_images"], cache["sample_original_frames"],
        cache["n_frames"], cache["n_datasets"], cfg)
end

"""
Generate fit overlay: boxes colored by fit quality using FilterConfig thresholds.

Colors: green=pass, red=photons fail, orange=precision fail, purple=pvalue fail, gray=no match.
"""
function _save_fit_overlay(dir, smld, sample_roi_batch, sample_images, sample_original_frames,
                           n_frames, n_datasets, cfg::FilterConfig)
    isempty(sample_roi_batch) && return

    # Map sample index (1:N) to absolute frame number
    sample_to_original = Dict(i => f for (i, f) in enumerate(sample_original_frames))
    original_frame_set = Set(sample_original_frames)

    # Build absolute frame for each emitter and filter to sampled frames
    emitter_abs_frames = [(e, (e.dataset - 1) * n_frames + e.frame) for e in smld.emitters]
    sample_emitters = [(e, af) for (e, af) in emitter_abs_frames if af in original_frame_set]
    isempty(sample_emitters) && return

    # Pre-resolve PSF sigma bounds (needed for box coloring)
    psf_bounds = _resolve_psf_bounds(smld.emitters, cfg)

    # Color each ROI box by matching emitter quality
    fit_colors = Symbol[]
    pix_size = smld.camera.pixel_edges_x[2] - smld.camera.pixel_edges_x[1]
    x_origin = smld.camera.pixel_edges_x[1]
    y_origin = smld.camera.pixel_edges_y[1]

    for i in 1:length(sample_roi_batch)
        roi_frame_idx = sample_roi_batch.frame_indices[i]
        original_frame = sample_to_original[roi_frame_idx]
        roi_x = sample_roi_batch.x_corners[i] + sample_roi_batch.roi_size ÷ 2
        roi_y = sample_roi_batch.y_corners[i] + sample_roi_batch.roi_size ÷ 2

        # Find closest matching emitter by position in this frame
        best_emitter = nothing
        best_dist = Inf
        for (e, af) in sample_emitters
            if af == original_frame
                ex_px = (e.x - x_origin) / pix_size
                ey_px = (e.y - y_origin) / pix_size
                dist = sqrt((ex_px - roi_x)^2 + (ey_px - roi_y)^2)
                if dist < best_dist
                    best_dist = dist
                    best_emitter = e
                end
            end
        end

        if best_emitter !== nothing && best_dist < sample_roi_batch.roi_size
            push!(fit_colors, _fit_box_color(best_emitter, cfg, psf_bounds))
        else
            push!(fit_colors, :gray)
        end
    end

    title = _fit_overlay_title(cfg)
    _save_box_overlay(dir, "fit_overlay.png", sample_images, sample_roi_batch, fit_colors;
                      title_prefix="Frame", frame_labels=sample_original_frames, suptitle=title)
end

"""
    _resolve_psf_bounds(emitters, cfg) -> Union{NamedTuple, Nothing}

Pre-resolve PSF sigma filter bounds from the emitter population.
Returns named tuple with :iso or :aniso bounds, or nothing if psf_sigma filter is off.
"""
function _resolve_psf_bounds(emitters, cfg::FilterConfig)
    cfg.psf_sigma === nothing && return nothing
    isempty(emitters) && return nothing

    if hasproperty(emitters[1], :σ)
        lo, hi = _get_psf_sigma_bounds(cfg.psf_sigma, [e.σ for e in emitters])
        (lo > 0 && hi > 0) || return nothing
        return (kind=:iso, lo=lo, hi=hi)
    elseif hasproperty(emitters[1], :σx) && hasproperty(emitters[1], :σy)
        lo_x, hi_x = _get_psf_sigma_bounds(cfg.psf_sigma, [e.σx for e in emitters])
        lo_y, hi_y = _get_psf_sigma_bounds(cfg.psf_sigma, [e.σy for e in emitters])
        (lo_x > 0 && hi_x > 0 && lo_y > 0 && hi_y > 0) || return nothing
        return (kind=:aniso, lo_x=lo_x, hi_x=hi_x, lo_y=lo_y, hi_y=hi_y)
    end
    return nothing
end

"""Determine box color based on FilterConfig thresholds. All-nothing filters = all green."""
function _fit_box_color(e, cfg::FilterConfig, psf_bounds)
    if cfg.photons !== nothing
        e.photons < cfg.photons[1] && return :red
        e.photons > cfg.photons[2] && return :red
    end
    if cfg.precision !== nothing
        max(e.σ_x, e.σ_y) > cfg.precision[2] && return :orange
        max(e.σ_x, e.σ_y) < cfg.precision[1] && return :orange
    end
    if cfg.pvalue !== nothing
        e.pvalue < cfg.pvalue[1] && return :purple
        e.pvalue > cfg.pvalue[2] && return :purple
    end
    if psf_bounds !== nothing
        if psf_bounds.kind === :iso
            (psf_bounds.lo <= e.σ <= psf_bounds.hi) || return :cyan
        elseif psf_bounds.kind === :aniso
            (psf_bounds.lo_x <= e.σx <= psf_bounds.hi_x && psf_bounds.lo_y <= e.σy <= psf_bounds.hi_y) || return :cyan
        end
    end
    return :green
end

"""Build title string from actual FilterConfig thresholds."""
function _fit_overlay_title(cfg::FilterConfig)
    parts = String["green=pass"]
    if cfg.photons !== nothing
        lo = cfg.photons[1]
        hi = cfg.photons[2]
        hi_str = hi == Inf ? "" : "/>$(round(Int, hi))"
        push!(parts, "red=photons<$(round(Int, lo))$hi_str")
    end
    if cfg.precision !== nothing
        hi_nm = round(cfg.precision[2] * 1000, digits=1)
        push!(parts, "orange=prec>$(hi_nm)nm")
    end
    if cfg.pvalue !== nothing
        push!(parts, "purple=pval<$(cfg.pvalue[1])")
    end
    if cfg.psf_sigma !== nothing
        push!(parts, "cyan=psf_σ fail")
    end
    push!(parts, "gray=no match")
    "Fit: " * join(parts, "  ")
end
