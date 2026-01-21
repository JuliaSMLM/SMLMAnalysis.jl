"""
Fitting step - wraps GaussMLE.fit
"""

@kwdef struct FitConfig <: StepConfig
    name::String = "fit"
    # GaussMLE kwargs
    psf_model::Symbol = :variable  # :fixed, :variable, :anisotropic
    psf_sigma::Float32 = 0.135f0   # For :fixed only
    iterations::Int = 20
    device::Union{Symbol,Nothing} = nothing  # nothing=auto, :cpu, :gpu
    # Extra
    verbose::Int = Verbosity.STANDARD
end

function run_step!(a::Analysis, cfg::FitConfig)
    a.roi_batch === nothing && error("Must run Detect first")
    a.step_counter += 1
    v = _get_verbose(a, cfg)
    dir = _stepdir(a, cfg)

    v >= Verbosity.PROGRESS && @info "[$(a.step_counter)] $(cfg.name)" psf_model=cfg.psf_model iterations=cfg.iterations

    psf = if cfg.psf_model == :fixed
        GaussianXYNB(cfg.psf_sigma)
    elseif cfg.psf_model == :variable
        GaussianXYNBS()
    elseif cfg.psf_model == :anisotropic
        GaussianXYNBSXSY()
    else
        error("Unknown psf_model: $(cfg.psf_model). Use :fixed, :variable, or :anisotropic")
    end

    fitter = GaussMLEFitter(psf_model=psf, iterations=cfg.iterations, device=cfg.device)
    t = @elapsed smld = GaussMLE.fit(fitter, a.roi_batch)

    # Assign dataset from roi_datasets (set during detection)
    # Frames are already per-dataset (1:n_frames_per_dataset) from looped detection
    if a.roi_datasets !== nothing && a.n_datasets > 1
        smld = _assign_datasets_from_roi(smld, a.roi_datasets, a.n_frames_per_dataset, a.n_datasets)
    end

    a.smld_raw = smld
    a.smld = smld

    n_fits = length(smld.emitters)
    summary = Dict{Symbol,Any}(:n_fits => n_fits, :psf_model => cfg.psf_model)
    _record!(a, cfg, t, summary)
    _checkpoint!(a)  # Auto-checkpoint after expensive step

    if dir !== nothing
        _save_step_outputs!(dir, a, cfg, v, t)
    end

    v >= Verbosity.PROGRESS && @info "  → $n_fits fits ($(round(t, digits=2))s)"
    a
end

function _save_step_outputs!(dir::String, a::Analysis, cfg::FitConfig, v::Int, t::Float64)
    mkpath(dir)
    _save_config!(dir, cfg)

    images = get_images(a.data)

    if v >= Verbosity.STANDARD
        _write_fit_stats(dir, a.smld, cfg, t)
        _save_fit_figures(dir, a.smld, a.roi_batch, images, cfg)
    end

    if v >= Verbosity.DETAILED
        _save_fit_detailed(dir, a.smld, a.roi_batch, cfg)
    end

    if v >= Verbosity.DEBUG
        _save_fit_debug(dir, a.smld, a.roi_batch, cfg)
    end
end

function _write_fit_stats(dir, smld, cfg, t)
    emitters = smld.emitters
    n = length(emitters)

    photons = [e.photons for e in emitters]
    bg = [e.bg for e in emitters]
    σ_x = [e.σ_x for e in emitters]
    σ_y = [e.σ_y for e in emitters]
    pvalue = [e.pvalue for e in emitters]

    filepath = joinpath(dir, "stats.md")
    open(filepath, "w") do io
        println(io, "# Fitting Statistics\n")
        println(io, "## Summary")
        println(io, "- **Fits**: $n")
        println(io, "- **Time**: $(round(t, digits=2))s ($(round(n/t/1000, digits=1))k fits/s)")
        println(io, "- **Model**: $(cfg.psf_model)")
        println(io, "")
        println(io, "## Distributions\n")
        println(io, "| Parameter | Median | 5% | 95% |")
        println(io, "|-----------|--------|-----|-----|")
        println(io, "| Photons | $(round(median(photons), digits=0)) | $(round(quantile(photons, 0.05), digits=0)) | $(round(quantile(photons, 0.95), digits=0)) |")
        println(io, "| Background | $(round(median(bg), digits=1)) | $(round(quantile(bg, 0.05), digits=1)) | $(round(quantile(bg, 0.95), digits=1)) |")
        println(io, "| σ_x (nm) | $(round(median(σ_x)*1000, digits=1)) | $(round(quantile(σ_x, 0.05)*1000, digits=1)) | $(round(quantile(σ_x, 0.95)*1000, digits=1)) |")
        println(io, "| σ_y (nm) | $(round(median(σ_y)*1000, digits=1)) | $(round(quantile(σ_y, 0.05)*1000, digits=1)) | $(round(quantile(σ_y, 0.95)*1000, digits=1)) |")
        println(io, "")
        println(io, "## P-value")
        pval_pass = sum(pvalue .> 0.001) / n
        println(io, "- pvalue > 0.001: $(round(100*pval_pass, digits=1))%")
        println(io, "- pvalue > 0.01: $(round(100*sum(pvalue .> 0.01)/n, digits=1))%")
    end
end

function _save_fit_figures(dir, smld, roi_batch, images, cfg)
    emitters = smld.emitters
    isempty(emitters) && return

    photons = [e.photons for e in emitters]
    bg = [e.bg for e in emitters]
    σ_x = [e.σ_x for e in emitters] .* 1000  # precision in nm
    σ_y = [e.σ_y for e in emitters] .* 1000
    pvalue = [e.pvalue for e in emitters]

    # Check for z precision
    has_z = hasproperty(emitters[1], :σ_z) && emitters[1].σ_z !== nothing
    σ_z = has_z ? [e.σ_z for e in emitters] .* 1000 : nothing

    # Check PSF model type
    has_psf_iso = hasproperty(emitters[1], :σ)
    has_psf_aniso = hasproperty(emitters[1], :σx)

    # Colors for consistent styling
    REJECTED_COLOR = (:gray85, 0.7)
    MEAN_COLOR = :blue
    MEDIAN_COLOR = :red
    THRESHOLD_COLOR = :black

    fig = Figure(size=(1200, 900))

    # ============================================================
    # Row 1: Photons and Background
    # ============================================================

    # Photons histogram
    p98 = quantile(photons, 0.98)
    ax1 = Axis(fig[1, 1], xlabel="Photons", ylabel="Count", title="Photon Distribution")

    # Gray rejected region (below typical photons threshold preview)
    # Preview default: (500, Inf) - reject below 500
    photons_lo_preview = 500.0
    vspan!(ax1, 0, photons_lo_preview, color=REJECTED_COLOR)

    hist!(ax1, photons[photons .<= p98], bins=50)
    vlines!(ax1, [mean(photons)], color=MEAN_COLOR, linestyle=:solid, linewidth=2, label="mean")
    vlines!(ax1, [median(photons)], color=MEDIAN_COLOR, linestyle=:dash, linewidth=2, label="median")
    vlines!(ax1, [photons_lo_preview], color=THRESHOLD_COLOR, linestyle=:dot, linewidth=2, label="threshold")
    xlims!(ax1, 0, p98)

    # Stats annotation
    text!(ax1, 0.97, 0.95, text="mean: $(round(Int, mean(photons)))\nmedian: $(round(Int, median(photons)))",
          align=(:right, :top), space=:relative, fontsize=10)

    # Background histogram
    bg98 = quantile(bg, 0.98)
    ax2 = Axis(fig[1, 2], xlabel="Background", ylabel="Count", title="Background Distribution")
    hist!(ax2, bg[bg .<= bg98], bins=50)
    vlines!(ax2, [mean(bg)], color=MEAN_COLOR, linestyle=:solid, linewidth=2)
    vlines!(ax2, [median(bg)], color=MEDIAN_COLOR, linestyle=:dash, linewidth=2)
    xlims!(ax2, 0, bg98)

    text!(ax2, 0.97, 0.95, text="mean: $(round(mean(bg), digits=1))\nmedian: $(round(median(bg), digits=1))",
          align=(:right, :top), space=:relative, fontsize=10)

    # ============================================================
    # Row 2: Localization Precision (combined) and P-value
    # ============================================================

    # Combined precision histogram (σ_x, σ_y, and σ_z if present)
    # Preview default: (0, 0.015) - reject above 15nm
    prec_hi_preview = 15.0  # nm
    prec_data = vcat(σ_x, σ_y)
    prec_labels = ["σ_x", "σ_y"]
    if has_z && σ_z !== nothing
        prec_data = vcat(prec_data, σ_z)
        push!(prec_labels, "σ_z")
    end
    prec98 = quantile(prec_data, 0.98)

    ax3 = Axis(fig[2, 1], xlabel="Localization Precision (nm)", ylabel="Count",
               title="Precision Distribution")

    # Ensure rejected region is visible by extending xlims if needed
    prec_xlim = max(prec98, prec_hi_preview * 1.2)

    # Gray rejected region (above max_precision)
    vspan!(ax3, prec_hi_preview, prec_xlim, color=REJECTED_COLOR)

    hist!(ax3, σ_x[σ_x .<= prec98], bins=50, color=(:blue, 0.5), label="σ_x")
    hist!(ax3, σ_y[σ_y .<= prec98], bins=50, color=(:red, 0.5), label="σ_y")
    if has_z && σ_z !== nothing
        hist!(ax3, σ_z[σ_z .<= prec98], bins=50, color=(:green, 0.5), label="σ_z")
    end

    vlines!(ax3, [prec_hi_preview], color=THRESHOLD_COLOR, linestyle=:dot, linewidth=2)
    xlims!(ax3, 0, prec_xlim)
    axislegend(ax3, position=:rt, framevisible=false, labelsize=9)

    text!(ax3, 0.97, 0.70, text="σ_x: $(round(median(σ_x), digits=1)) nm\nσ_y: $(round(median(σ_y), digits=1)) nm" *
          (has_z ? "\nσ_z: $(round(median(σ_z), digits=1)) nm" : ""),
          align=(:right, :top), space=:relative, fontsize=10)

    # P-value histogram
    # Preview default: (1e-3, 1) - reject below 0.001
    # For perfect data with correct model, p-values are uniform on [0,1]
    # In -log10 space (u = -log10(p)), the PDF is: f(u) = ln(10) * 10^(-u)
    ax4 = Axis(fig[2, 2], xlabel="log₁₀(p-value)", ylabel="Density", title="P-value Distribution")
    pval_nonzero = pvalue[pvalue .> 0]
    pval_lo_preview = 1e-3

    if !isempty(pval_nonzero)
        log_pval = log10.(pval_nonzero)
        pval_lo = quantile(log_pval, 0.02)

        # Gray rejected region (below min_pvalue)
        vspan!(ax4, pval_lo - 1, log10(pval_lo_preview), color=REJECTED_COLOR)

        # Histogram normalized to density (for comparison with theory curve)
        hist!(ax4, log_pval[log_pval .>= pval_lo], bins=50, normalization=:pdf, color=(:steelblue, 0.7))

        # Theoretical curve: if p ~ Uniform[0,1], then u = -log10(p) has PDF = ln(10) * 10^(-u)
        # Note: log_pval is log10(p), so u = -log_pval
        u_range = range(0, -pval_lo, length=100)
        theory_pdf = log(10) .* (10.0 .^ (-u_range))
        lines!(ax4, -u_range, theory_pdf, color=:red, linewidth=2, label="Uniform theory")

        vlines!(ax4, [mean(log_pval)], color=MEAN_COLOR, linestyle=:solid, linewidth=1.5)
        vlines!(ax4, [median(log_pval)], color=MEDIAN_COLOR, linestyle=:dash, linewidth=1.5)
        vlines!(ax4, [log10(pval_lo_preview)], color=THRESHOLD_COLOR, linestyle=:dot, linewidth=2)
        xlims!(ax4, pval_lo, 0)

        pval_pass_pct = round(100 * sum(pvalue .> pval_lo_preview) / length(pvalue), digits=1)
        text!(ax4, 0.03, 0.95, text="pass: $(pval_pass_pct)%\nthreshold: $(pval_lo_preview)",
              align=(:left, :top), space=:relative, fontsize=10)
    end

    # ============================================================
    # Row 3: PSF Sigma
    # ============================================================

    if has_psf_aniso
        # Anisotropic: overlay σx and σy
        psf_σx = [e.σx for e in emitters] .* 1000
        psf_σy = [e.σy for e in emitters] .* 1000
        psf_data = vcat(psf_σx, psf_σy)
        psf98 = quantile(psf_data, 0.98)
        psf02 = quantile(psf_data, 0.02)

        # Calculate mode for threshold preview
        mode_x = _calculate_mode([e.σx for e in emitters]) * 1000
        mode_y = _calculate_mode([e.σy for e in emitters]) * 1000
        mode_avg = (mode_x + mode_y) / 2
        tol = 0.10  # default psf_sigma_range :auto uses mode ± 10%

        ax5 = Axis(fig[3, 1], xlabel="Fitted PSF σ (nm)", ylabel="Count", title="PSF Width Distribution")

        # Gray rejected regions (outside mode ± tol)
        vspan!(ax5, psf02 * 0.9, mode_avg * (1 - tol), color=REJECTED_COLOR)
        vspan!(ax5, mode_avg * (1 + tol), psf98 * 1.1, color=REJECTED_COLOR)

        hist!(ax5, psf_σx[(psf_σx .>= psf02) .& (psf_σx .<= psf98)], bins=50, color=(:blue, 0.5), label="σx")
        hist!(ax5, psf_σy[(psf_σy .>= psf02) .& (psf_σy .<= psf98)], bins=50, color=(:red, 0.5), label="σy")

        vlines!(ax5, [mode_avg * (1 - tol), mode_avg * (1 + tol)], color=THRESHOLD_COLOR, linestyle=:dot, linewidth=2)
        xlims!(ax5, psf02, psf98)
        axislegend(ax5, position=:rt, framevisible=false, labelsize=9)

        text!(ax5, 0.97, 0.70, text="mode σx: $(round(mode_x, digits=1)) nm\nmode σy: $(round(mode_y, digits=1)) nm\ntol: ±$(round(Int, tol*100))%",
              align=(:right, :top), space=:relative, fontsize=10)

    elseif has_psf_iso
        # Isotropic: single σ
        psf_σ = [e.σ for e in emitters] .* 1000
        psf98 = quantile(psf_σ, 0.98)
        psf02 = quantile(psf_σ, 0.02)

        mode_σ = _calculate_mode([e.σ for e in emitters]) * 1000
        tol = 0.10  # default psf_sigma_range :auto uses mode ± 10%

        # Ensure xlims include both threshold bounds
        lo_bound = mode_σ * (1 - tol)
        hi_bound = mode_σ * (1 + tol)
        psf_xmin = min(psf02, lo_bound * 0.95)
        psf_xmax = max(psf98, hi_bound * 1.05)

        ax5 = Axis(fig[3, 1], xlabel="Fitted PSF σ (nm)", ylabel="Count", title="PSF Width Distribution")

        # Gray rejected regions (outside mode ± tol)
        vspan!(ax5, psf_xmin, lo_bound, color=REJECTED_COLOR)
        vspan!(ax5, hi_bound, psf_xmax, color=REJECTED_COLOR)

        hist!(ax5, psf_σ[(psf_σ .>= psf02) .& (psf_σ .<= psf98)], bins=50)
        vlines!(ax5, [mean(psf_σ)], color=MEAN_COLOR, linestyle=:solid, linewidth=2)
        vlines!(ax5, [median(psf_σ)], color=MEDIAN_COLOR, linestyle=:dash, linewidth=2)
        vlines!(ax5, [lo_bound, hi_bound], color=THRESHOLD_COLOR, linestyle=:dot, linewidth=2)
        xlims!(ax5, psf_xmin, psf_xmax)

        text!(ax5, 0.97, 0.95, text="mode: $(round(mode_σ, digits=1)) nm\nmean: $(round(mean(psf_σ), digits=1)) nm\ntol: ±$(round(Int, tol*100))%",
              align=(:right, :top), space=:relative, fontsize=10)
    else
        # Fixed PSF: show vertical line at fixed sigma
        ax5 = Axis(fig[3, 1], xlabel="PSF σ (nm)", ylabel="", title="PSF Width (Fixed)")
        fixed_σ = cfg.psf_sigma * 1000
        vlines!(ax5, [fixed_σ], color=:blue, linewidth=3)
        text!(ax5, 0.5, 0.5, text="Fixed: $(round(fixed_σ, digits=1)) nm",
              align=(:center, :center), space=:relative, fontsize=14)
        hideydecorations!(ax5)
    end

    # Row 3, Col 2: Legend / Summary
    ax6 = Axis(fig[3, 2], title="Legend")
    hidedecorations!(ax6)
    hidespines!(ax6)

    # Draw legend items manually
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

    # Frame overlay with boxes colored by fit status (like detect panel)
    _save_fit_overlay(dir, smld, roi_batch, images)
end

"""
Determine box color based on fit status.
Precedence (highest to lowest): photons (red) > precision (orange) > pvalue (purple) > accepted (green)
"""
function _fit_box_color(e; min_photons=500.0, max_precision=0.015, min_pvalue=1e-6)
    # Check in precedence order (highest first)
    e.photons < min_photons && return :red           # failed photons
    prec = sqrt(e.σ_x^2 + e.σ_y^2) / sqrt(2)
    prec > max_precision && return :orange           # failed precision
    e.pvalue < min_pvalue && return :purple          # failed pvalue
    return :green                                     # accepted
end

"""
Save fit overlay showing frame images with boxes colored by fit status.
Colors: green=accepted, purple=failed pvalue, orange=failed precision, red=failed photons
"""
function _save_fit_overlay(dir, smld, roi_batch, images)
    # Skip if roi_batch is not available (e.g., when resuming from checkpoint)
    roi_batch === nothing && return

    emitters = smld.emitters

    # Compute box color for each emitter based on fit status
    box_colors = [_fit_box_color(e) for e in emitters]

    _save_box_overlay(dir, "fit_overlay.png", images, roi_batch, box_colors; title_prefix="Frame")
end

function _save_fit_detailed(dir, smld, roi_batch, cfg)
    # Reserved for future detailed diagnostics
end

function _save_fit_debug(dir, smld, roi_batch, cfg)
    # Save individual ROI fits montage with residuals
    # TODO: implement when needed
end

"""
Assign dataset field to emitters using roi_datasets vector from detection.
Frames are already per-dataset (1:n_frames_per_dataset) from looped detection.
"""
function _assign_datasets_from_roi(smld::BasicSMLD, roi_datasets::Vector{Int},
                                    n_frames_per_dataset::Int, n_datasets::Int)
    length(roi_datasets) == length(smld.emitters) ||
        error("roi_datasets length ($(length(roi_datasets))) != emitters length ($(length(smld.emitters)))")

    new_emitters = map(enumerate(smld.emitters)) do (i, e)
        ds = roi_datasets[i]
        _with_dataset(e, ds)
    end

    BasicSMLD(new_emitters, smld.camera, n_frames_per_dataset, n_datasets, smld.metadata)
end

"""Update emitter's dataset field using struct reconstruction"""
function _with_dataset(e::Emitter2DFit{T}, ds::Int) where T
    Emitter2DFit{T}(
        e.x, e.y, e.photons, e.bg, e.σ_x, e.σ_y,
        e.frame, ds, e.track_id, e.pvalue, e.σ, e.σx, e.σy
    )
end

function _with_dataset(e::Emitter2D{T}, ds::Int) where T
    Emitter2D{T}(e.x, e.y, e.photons, e.σ_x, e.σ_y, e.frame, ds, e.track_id)
end

function _with_dataset(e::GaussMLE.Emitter2DFitSigma{T}, ds::Int) where T
    GaussMLE.Emitter2DFitSigma{T}(
        e.x, e.y, e.photons, e.bg, e.σ, e.σ_x, e.σ_y, e.σ_photons, e.σ_bg, e.σ_σ,
        e.pvalue, e.frame, ds, e.track_id, e.id
    )
end
