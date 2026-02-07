"""
    writers.jl

Statistics writing functions for SMLM analysis pipeline.
All functions write markdown files for human and LLM diagnostics.
"""

using Statistics

# =============================================================================
# Detection Stats
# =============================================================================

"""Write detection statistics markdown file."""
function write_detection_stats(roi_batch, data, config, elapsed_time)
    nframes = size(data, 3)
    n_rois = length(roi_batch)

    rois_per_frame = [sum(roi_batch.frame_indices .== f) for f in 1:nframes]

    roi_intensities = [maximum(roi_batch.data[:,:,i]) for i in 1:min(n_rois, 10000)]
    bg_estimates = [minimum(roi_batch.data[:,:,i]) for i in 1:min(n_rois, 10000)]
    signal_above_bg = roi_intensities .- bg_estimates

    filepath = joinpath(config.outdir, "01_detection", "detection_stats.md")
    open(filepath, "w") do io
        println(io, "# Detection Statistics\n")
        println(io, "## Summary")
        println(io, "- **Total ROIs detected**: $(n_rois)")
        println(io, "- **Frames analyzed**: $(nframes)")
        println(io, "- **Time**: $(round(elapsed_time, digits=2))s ($(round(nframes/elapsed_time, digits=0)) frames/s)")
        println(io, "- **ROIs per frame**: mean=$(round(mean(rois_per_frame), digits=1)), std=$(round(std(rois_per_frame), digits=1))")
        println(io, "- **ROIs/frame range**: $(minimum(rois_per_frame)) - $(maximum(rois_per_frame))")
        println(io, "")
        println(io, "## Detection Parameters")
        println(io, "- Box size: $(config.boxsize)")
        println(io, "- Min photons threshold: $(config.detect_min_photons)")
        println(io, "- PSF sigma: $(round(config.psf_sigma * 1000, digits=0)) nm")
        println(io, "- Overlap: $(config.overlap)")
        println(io, "")
        println(io, "## ROI Intensity (sampled)")
        println(io, "- Peak intensity: median=$(round(median(roi_intensities), digits=0)) ADU")
        println(io, "- Background estimate: median=$(round(median(bg_estimates), digits=0)) ADU")
        println(io, "- Signal above background: median=$(round(median(signal_above_bg), digits=0)) ADU")
        println(io, "")
        println(io, "## Health Check")
        rpf_mean = mean(rois_per_frame)
        rpf_cv = std(rois_per_frame) / rpf_mean
        println(io, "- ROIs/frame CV: $(round(rpf_cv, digits=2)) ", rpf_cv < 0.3 ? "✓" : "⚠ (high variation)")
        println(io, "- Signal/background: $(round(median(signal_above_bg)/median(bg_estimates), digits=1))x ",
                median(signal_above_bg) > median(bg_estimates) ? "✓" : "⚠")
    end
end

# =============================================================================
# Fitting Stats
# =============================================================================

"""Write fitting statistics markdown file."""
function write_fitting_stats(smld_raw, config, elapsed_time)
    emitters = smld_raw.emitters
    n = length(emitters)

    photons = [e.photons for e in emitters]
    bg = [e.bg for e in emitters]
    σ_x = [e.σ_x for e in emitters]
    σ_y = [e.σ_y for e in emitters]
    pvalue = [e.pvalue for e in emitters]

    has_sigma = hasfield(typeof(emitters[1]), :σ)
    has_sigma_xy = hasfield(typeof(emitters[1]), :σx)

    filepath = joinpath(config.outdir, "02_fitting", "fitting_stats.md")
    open(filepath, "w") do io
        println(io, "# Fitting Statistics\n")
        println(io, "## Summary")
        println(io, "- **Total fits**: $(n)")
        println(io, "- **Time**: $(round(elapsed_time, digits=2))s ($(round(n/elapsed_time/1000, digits=0))k fits/s)")
        println(io, "- **Model**: $(config.fit_model)")
        println(io, "- **Iterations**: $(config.iterations)")
        println(io, "")

        println(io, "## Parameter Distributions\n")
        println(io, "| Parameter | 5% | 25% | 50% (median) | 75% | 95% |")
        println(io, "|-----------|-----|-----|--------------|-----|-----|")

        for (name, vals, unit, scale) in [
            ("Photons", photons, "", 1),
            ("Background", bg, "e⁻", 1),
            ("Precision σ_x", σ_x, "nm", 1000),
            ("Precision σ_y", σ_y, "nm", 1000),
        ]
            p = quantile(vals .* scale, [0.05, 0.25, 0.50, 0.75, 0.95])
            println(io, "| $name ($unit) | $(round(p[1], digits=1)) | $(round(p[2], digits=1)) | $(round(p[3], digits=1)) | $(round(p[4], digits=1)) | $(round(p[5], digits=1)) |")
        end

        if has_sigma
            σ = [e.σ for e in emitters]
            p = quantile(σ .* 1000, [0.05, 0.25, 0.50, 0.75, 0.95])
            println(io, "| PSF σ (nm) | $(round(p[1], digits=1)) | $(round(p[2], digits=1)) | $(round(p[3], digits=1)) | $(round(p[4], digits=1)) | $(round(p[5], digits=1)) |")
        end

        if has_sigma_xy
            σx = [e.σx for e in emitters]
            σy = [e.σy for e in emitters]
            px = quantile(σx .* 1000, [0.05, 0.25, 0.50, 0.75, 0.95])
            py = quantile(σy .* 1000, [0.05, 0.25, 0.50, 0.75, 0.95])
            println(io, "| PSF σx (nm) | $(round(px[1], digits=1)) | $(round(px[2], digits=1)) | $(round(px[3], digits=1)) | $(round(px[4], digits=1)) | $(round(px[5], digits=1)) |")
            println(io, "| PSF σy (nm) | $(round(py[1], digits=1)) | $(round(py[2], digits=1)) | $(round(py[3], digits=1)) | $(round(py[4], digits=1)) | $(round(py[5], digits=1)) |")
        end

        println(io, "")
        println(io, "## P-value Distribution")
        n_zero = sum(pvalue .== 0)
        println(io, "- pvalue = 0: $(n_zero) ($(round(100*n_zero/n, digits=1))%)")
        println(io, "- pvalue > 0: $(sum(pvalue .> 0)) ($(round(100*sum(pvalue .> 0)/n, digits=1))%)")
        println(io, "- pvalue > 0.001: $(sum(pvalue .> 0.001)) ($(round(100*sum(pvalue .> 0.001)/n, digits=1))%)")
        println(io, "- pvalue > 0.01: $(sum(pvalue .> 0.01)) ($(round(100*sum(pvalue .> 0.01)/n, digits=1))%)")
        println(io, "- pvalue > 0.05: $(sum(pvalue .> 0.05)) ($(round(100*sum(pvalue .> 0.05)/n, digits=1))%)")
        println(io, "")
        println(io, "## Health Check")
        pval_pass = sum(pvalue .> 0.001) / n
        prec_med = median(σ_x) * 1000
        println(io, "- Precision σ_x median: $(round(prec_med, digits=1)) nm ", prec_med < 20 ? "✓" : "⚠")
        println(io, "- pvalue > 0.001: $(round(100*pval_pass, digits=1))% ", pval_pass > 0.05 ? "✓" : "⚠ (low)")
        println(io, "- Photons median: $(round(median(photons), digits=0)) ", 1000 < median(photons) < 100000 ? "✓" : "⚠")
    end
end

# =============================================================================
# Filter Stats
# =============================================================================

"""Write filter statistics markdown file."""
function write_filter_stats(smld_raw, smld_filtered, config, elapsed_time, calculate_mode_fn)
    n_raw = length(smld_raw.emitters)
    n_filtered = length(smld_filtered.emitters)

    emitters = smld_raw.emitters
    photons = [e.photons for e in emitters]
    σ_x = [e.σ_x for e in emitters]
    σ_y = [e.σ_y for e in emitters]
    pvalue = [e.pvalue for e in emitters]

    has_sigma = hasfield(typeof(emitters[1]), :σ)
    has_sigma_xy = hasfield(typeof(emitters[1]), :σx) && hasfield(typeof(emitters[1]), :σy)

    psf_sigma_mode = 0.0
    psf_sigma_mode_x = 0.0
    psf_sigma_mode_y = 0.0

    if has_sigma
        psf_sigmas = [e.σ for e in emitters]
        psf_sigma_mode = calculate_mode_fn(psf_sigmas)
    elseif has_sigma_xy
        psf_sigmas_x = [e.σx for e in emitters]
        psf_sigmas_y = [e.σy for e in emitters]
        psf_sigma_mode_x = calculate_mode_fn(psf_sigmas_x)
        psf_sigma_mode_y = calculate_mode_fn(psf_sigmas_y)
    end

    photon_ok = config.min_photons === nothing ? trues(n_raw) : photons .> config.min_photons
    precision_ok = config.max_precision === nothing ? trues(n_raw) :
                   [max(e.σ_x, e.σ_y) < config.max_precision for e in emitters]

    if has_sigma && config.psf_sigma_mode_tolerance !== nothing && psf_sigma_mode > 0
        tol = config.psf_sigma_mode_tolerance
        lo = psf_sigma_mode * (1 - tol)
        hi = psf_sigma_mode * (1 + tol)
        psf_sigma_ok = [lo <= e.σ <= hi for e in emitters]
    elseif has_sigma_xy && config.psf_sigma_mode_tolerance !== nothing && psf_sigma_mode_x > 0 && psf_sigma_mode_y > 0
        tol = config.psf_sigma_mode_tolerance
        lo_x = psf_sigma_mode_x * (1 - tol)
        hi_x = psf_sigma_mode_x * (1 + tol)
        lo_y = psf_sigma_mode_y * (1 - tol)
        hi_y = psf_sigma_mode_y * (1 + tol)
        psf_sigma_ok = [lo_x <= e.σx <= hi_x && lo_y <= e.σy <= hi_y for e in emitters]
    else
        psf_sigma_ok = trues(n_raw)
    end

    pvalue_ok = config.min_pvalue === nothing ? trues(n_raw) : pvalue .> config.min_pvalue

    n_photon_pass = sum(photon_ok)
    n_precision_pass = sum(precision_ok)
    n_psf_sigma_pass = sum(psf_sigma_ok)
    n_pvalue_pass = sum(pvalue_ok)

    filepath = joinpath(config.outdir, "03_filtered", "filter_stats.md")
    open(filepath, "w") do io
        println(io, "# Filter Statistics\n")
        println(io, "## Summary")
        println(io, "- **Input**: $(n_raw) localizations")
        println(io, "- **Output**: $(n_filtered) localizations")
        println(io, "- **Acceptance rate**: $(round(100*n_filtered/n_raw, digits=1))%")
        println(io, "- **Time**: $(round(elapsed_time*1000, digits=1))ms")
        println(io, "")

        if has_sigma && config.psf_sigma_mode_tolerance !== nothing && psf_sigma_mode > 0
            psf_sigma_mode_nm = psf_sigma_mode * 1000
            println(io, "## PSF Sigma Mode Analysis")
            println(io, "- **PSF sigma mode**: $(round(psf_sigma_mode_nm, digits=1)) nm")
            tol_pct = round(config.psf_sigma_mode_tolerance * 100, digits=0)
            lo_nm = round(psf_sigma_mode_nm * (1 - config.psf_sigma_mode_tolerance), digits=1)
            hi_nm = round(psf_sigma_mode_nm * (1 + config.psf_sigma_mode_tolerance), digits=1)
            println(io, "- **Acceptance range**: $(lo_nm) - $(hi_nm) nm (mode ± $(tol_pct)%)")
            println(io, "")
        elseif has_sigma_xy && config.psf_sigma_mode_tolerance !== nothing && psf_sigma_mode_x > 0 && psf_sigma_mode_y > 0
            psf_sigma_mode_x_nm = psf_sigma_mode_x * 1000
            psf_sigma_mode_y_nm = psf_sigma_mode_y * 1000
            println(io, "## PSF Sigma Mode Analysis (Anisotropic)")
            println(io, "- **PSF σx mode**: $(round(psf_sigma_mode_x_nm, digits=1)) nm")
            println(io, "- **PSF σy mode**: $(round(psf_sigma_mode_y_nm, digits=1)) nm")
            println(io, "")
        end

        println(io, "## Per-Filter Results\n")
        println(io, "| Filter | Threshold | Pass | Fail | % Pass |")
        println(io, "|--------|-----------|------|------|--------|")

        if config.min_photons !== nothing
            println(io, "| Photons | >$(config.min_photons) | $(n_photon_pass) | $(n_raw - n_photon_pass) | $(round(100*n_photon_pass/n_raw, digits=1))% |")
        end
        if config.max_precision !== nothing
            println(io, "| Precision | <$(config.max_precision*1000)nm | $(n_precision_pass) | $(n_raw - n_precision_pass) | $(round(100*n_precision_pass/n_raw, digits=1))% |")
        end
        if config.psf_sigma_mode_tolerance !== nothing && (psf_sigma_mode > 0 || (psf_sigma_mode_x > 0 && psf_sigma_mode_y > 0))
            tol_pct = round(config.psf_sigma_mode_tolerance * 100, digits=0)
            println(io, "| PSF σ | mode ±$(tol_pct)% | $(n_psf_sigma_pass) | $(n_raw - n_psf_sigma_pass) | $(round(100*n_psf_sigma_pass/n_raw, digits=1))% |")
        end
        if config.min_pvalue !== nothing
            println(io, "| P-value | >$(config.min_pvalue) | $(n_pvalue_pass) | $(n_raw - n_pvalue_pass) | $(round(100*n_pvalue_pass/n_raw, digits=1))% |")
        end

        println(io, "")
        println(io, "## Health Check")
        acc_rate = n_filtered / n_raw
        println(io, "- Acceptance rate: $(round(100*acc_rate, digits=1))% ",
                0.01 < acc_rate < 0.5 ? "✓" : (acc_rate < 0.01 ? "⚠ (too strict)" : "⚠ (too loose)"))
    end
end

# =============================================================================
# Calibration Stats
# =============================================================================

"""Write calibration statistics to markdown file."""
function write_calibration_stats(σ_motion_nm::Float64, k_scale::Float64, cal, config, elapsed::Float64)
    filepath = joinpath(config.outdir, "05_calibration", "calibration_stats.md")

    open(filepath, "w") do io
        println(io, "# Uncertainty Calibration Applied")
        println(io)
        println(io, "## Calibration Parameters")
        println(io, "| Parameter | Value | Description |")
        println(io, "|-----------|-------|-------------|")
        println(io, "| σ_motion | $(round(σ_motion_nm, digits=1)) nm | Frame-to-frame motion/vibration |")
        println(io, "| k (CRLB scale) | $(round(k_scale, digits=2)) | Multiply σ_CRLB by this |")
        println(io)
        println(io, "## Correction Formula")
        println(io, "```")
        println(io, "σ_corrected = √(σ_motion² + k² × σ_CRLB²)")
        println(io, "```")
        println(io)
        println(io, "## Effect")
        println(io, "- All localization uncertainties have been adjusted using this model")
        println(io, "- Frame-connected track uncertainties recalculated with corrected weights")
        println(io, "- Time: $(round(elapsed, digits=2))s")
    end
end

# =============================================================================
# Drift Stats
# =============================================================================

"""Write drift correction statistics markdown file."""
function write_drift_stats(drift_model, smld, config, elapsed_time, DC)
    n_datasets = length(drift_model.intra)
    n_frames = smld.n_frames
    frames = collect(1:n_frames)

    filepath = joinpath(config.outdir, "06_drift", "drift_stats.md")
    open(filepath, "w") do io
        println(io, "# Drift Correction Statistics\n")
        println(io, "## Summary")
        println(io, "- **Time**: $(round(elapsed_time, digits=2))s")
        println(io, "- **Frames**: $n_frames")
        println(io, "- **Datasets**: $n_datasets")

        if n_datasets == 1
            drift_x = [DC.applydrift(0.0, f, drift_model.intra[1].dm[1]) for f in frames]
            drift_y = [DC.applydrift(0.0, f, drift_model.intra[1].dm[2]) for f in frames]
            drift_x_nm = drift_x .* 1000
            drift_y_nm = drift_y .* 1000

            max_x = maximum(abs.(drift_x_nm))
            max_y = maximum(abs.(drift_y_nm))
            total_x = drift_x_nm[end] - drift_x_nm[1]
            total_y = drift_y_nm[end] - drift_y_nm[1]
            total_dist = sqrt(total_x^2 + total_y^2)

            println(io, "- **Max X drift**: $(round(max_x, digits=1)) nm")
            println(io, "- **Max Y drift**: $(round(max_y, digits=1)) nm")
            println(io, "- **Total displacement**: $(round(total_dist, digits=1)) nm")
            println(io, "")
            println(io, "## Parameters")
            println(io, "- Model: $(config.drift_model)")
            println(io, "- Degree: $(config.drift_degree)")
            println(io, "- Cost function: $(config.drift_cost_fun)")
        end
    end
end

# =============================================================================
# Frame Connection Stats
# =============================================================================

"""Write frame connection statistics markdown file."""
function write_frameconnect_stats(n_before, n_after, smld_connected, fc_params, config, elapsed_time; drift_analysis=nothing)
    emitters = smld_connected.emitters
    track_ids = [e.track_id for e in emitters]
    counts = Dict{Int, Int}()
    for tid in track_ids
        counts[tid] = get(counts, tid, 0) + 1
    end
    locs_per_track = collect(values(counts))

    n_tracks = length(locs_per_track)
    compression = n_before / n_after

    filepath = joinpath(config.outdir, "04_frameconnect", "frameconnect_stats.md")
    open(filepath, "w") do io
        println(io, "# Frame Connection Statistics\n")
        println(io, "## Summary")
        println(io, "- **Input localizations**: $(n_before)")
        println(io, "- **Output tracks**: $(n_tracks)")
        println(io, "- **Compression**: $(round(compression, digits=1))x")
        println(io, "- **Time**: $(round(elapsed_time, digits=2))s")
        println(io, "")

        println(io, "## Track Statistics")
        if !isempty(locs_per_track)
            println(io, "- Median locs/track: $(median(locs_per_track))")
            println(io, "- Max locs/track: $(maximum(locs_per_track))")
            println(io, "- Single-loc tracks: $(sum(locs_per_track .== 1)) ($(round(100*sum(locs_per_track .== 1)/n_tracks, digits=1))%)")
        end

        if drift_analysis !== nothing
            println(io, "")
            println(io, "## Uncertainty Validation")
            println(io, "- Mean χ²: $(round(drift_analysis.mean_chi2, digits=2)) (expected: 2.0)")
            println(io, "- Pairs analyzed: $(drift_analysis.n_pairs_total)")

            cal = drift_analysis.calibration
            if !isnan(cal.A)
                σ_motion = cal.A > 0 ? sqrt(cal.A) : 0.0
                k_scale = sqrt(cal.B)
                println(io, "")
                println(io, "## Calibration Model")
                println(io, "- σ_motion: $(round(σ_motion, digits=1)) nm")
                println(io, "- k (CRLB scale): $(round(k_scale, digits=2))")
                println(io, "- R²: $(round(cal.r_squared, digits=3))")
            end
        end
    end
end

# =============================================================================
# Render Stats
# =============================================================================

"""Write render statistics markdown file."""
function write_render_stats(smld, config, elapsed_time)
    emitters = smld.emitters
    n_locs = length(emitters)

    σ_x = [e.σ_x for e in emitters]
    σ_y = [e.σ_y for e in emitters]
    mean_precision = mean([sqrt(σ_x[i]^2 + σ_y[i]^2)/sqrt(2) for i in 1:n_locs]) * 1000
    nyquist_resolution = 2 * mean_precision

    x_range = maximum(e.x for e in emitters) - minimum(e.x for e in emitters)
    y_range = maximum(e.y for e in emitters) - minimum(e.y for e in emitters)
    area_um2 = x_range * y_range
    density = n_locs / area_um2

    filepath = joinpath(config.outdir, "08_superres", "render_stats.md")
    open(filepath, "w") do io
        println(io, "# Render Statistics\n")
        println(io, "## Summary")
        println(io, "- **Localizations**: $(n_locs)")
        println(io, "- **Area**: $(round(area_um2, digits=2)) μm²")
        println(io, "- **Density**: $(round(density, digits=1)) /μm²")
        println(io, "- **Time**: $(round(elapsed_time, digits=2))s")
        println(io, "")
        println(io, "## Renders Generated\n")
        println(io, "| Type | Zoom | Pixel Size | Purpose |")
        println(io, "|------|------|------------|---------|")
        if config.render_gaussian
            ps = round(1000 / config.render_gaussian_zoom * 0.078, digits=1)
            println(io, "| Gaussian (inferno) | $(config.render_gaussian_zoom)x | $(ps) nm | Primary super-res |")
        end
        if config.render_histogram
            ps = round(1000 / config.render_histogram_zoom * 0.078, digits=1)
            println(io, "| Histogram (time) | $(config.render_histogram_zoom)x | $(ps) nm | Temporal coverage |")
        end
        if config.render_circles
            ps = round(1000 / config.render_circles_zoom * 0.078, digits=1)
            println(io, "| Circles (time) | $(config.render_circles_zoom)x | $(ps) nm | Individual locs |")
        end
        println(io, "")
        println(io, "## Resolution Estimate")
        println(io, "- Mean precision: $(round(mean_precision, digits=1)) nm")
        println(io, "- Nyquist resolution (2×precision): $(round(nyquist_resolution, digits=1)) nm")
        println(io, "")
        println(io, "## Health Check")
        println(io, "- Localization density: $(round(density, digits=1)) /μm² ", density > 10 ? "✓" : "⚠ (sparse)")
    end
end

# =============================================================================
# Summary Stats
# =============================================================================

"""Write summary markdown file with health check."""
function write_summary(roi_batch, smld_raw, smld_filtered, data, config, timings)
    nframes = size(data, 3)
    n_rois = length(roi_batch)
    n_raw = length(smld_raw.emitters)
    n_filtered = length(smld_filtered.emitters)

    emitters_raw = smld_raw.emitters
    photons = [e.photons for e in emitters_raw]
    σ_x = [e.σ_x for e in emitters_raw]
    pvalue = [e.pvalue for e in emitters_raw]

    has_sigma = hasfield(typeof(emitters_raw[1]), :σ)
    has_sigma_xy = hasfield(typeof(emitters_raw[1]), :σx) && hasfield(typeof(emitters_raw[1]), :σy)
    psf_sigma = has_sigma ? median([e.σ for e in emitters_raw]) * 1000 : nothing
    psf_sigma_x = has_sigma_xy ? median([e.σx for e in emitters_raw]) * 1000 : nothing
    psf_sigma_y = has_sigma_xy ? median([e.σy for e in emitters_raw]) * 1000 : nothing

    rois_per_frame = n_rois / nframes
    pval_pass_rate = sum(pvalue .> 0.001) / n_raw
    acceptance_rate = n_filtered / n_raw

    filepath = joinpath(config.outdir, "summary.md")
    open(filepath, "w") do io
        println(io, "# SMLM Analysis Summary\n")
        println(io, "## Quick Health Check\n")
        println(io, "```")
        println(io, "Detection:    $(round(rois_per_frame, digits=1)) ROIs/frame     ", 50 < rois_per_frame < 500 ? "✓" : "⚠")
        println(io, "Fitting:      $(round(100*pval_pass_rate, digits=1))% pvalue>0.001  ", pval_pass_rate > 0.05 ? "✓" : "⚠")
        if psf_sigma !== nothing
            println(io, "PSF sigma:    $(round(psf_sigma, digits=0)) nm            ", 80 < psf_sigma < 250 ? "✓" : "⚠")
        elseif psf_sigma_x !== nothing && psf_sigma_y !== nothing
            avg_sigma = (psf_sigma_x + psf_sigma_y) / 2
            println(io, "PSF σx,σy:    $(round(psf_sigma_x, digits=0)),$(round(psf_sigma_y, digits=0)) nm     ", 80 < avg_sigma < 250 ? "✓" : "⚠")
        end
        println(io, "Precision:    $(round(median(σ_x)*1000, digits=1)) nm median      ", median(σ_x)*1000 < 20 ? "✓" : "⚠")
        println(io, "Photons:      $(round(median(photons), digits=0)) median      ", 1000 < median(photons) < 100000 ? "✓" : "⚠")
        println(io, "Filtering:    $(round(100*acceptance_rate, digits=1))% accepted      ", 0.01 < acceptance_rate < 0.5 ? "✓" : "⚠")
        println(io, "```")

        psf_sigma_ok = if psf_sigma !== nothing
            80 < psf_sigma < 250
        elseif psf_sigma_x !== nothing && psf_sigma_y !== nothing
            avg = (psf_sigma_x + psf_sigma_y) / 2
            80 < avg < 250
        else
            true
        end

        all_ok = (50 < rois_per_frame < 500) &&
                 (pval_pass_rate > 0.05) &&
                 psf_sigma_ok &&
                 (median(σ_x)*1000 < 20) &&
                 (1000 < median(photons) < 100000) &&
                 (0.01 < acceptance_rate < 0.5)
        println(io, "\n**STATUS: $(all_ok ? "HEALTHY" : "NEEDS ATTENTION")**\n")

        println(io, "## Pipeline Results\n")
        println(io, "| Step | Count | Time |")
        println(io, "|------|-------|------|")
        println(io, "| Detection | $(n_rois) ROIs | $(round(get(timings, "detection", 0.0), digits=1))s |")
        println(io, "| Fitting | $(n_raw) fits | $(round(get(timings, "fitting", 0.0), digits=1))s |")
        println(io, "| Filtering | $(n_filtered) kept ($(round(100*acceptance_rate, digits=1))%) | $(round(get(timings, "filtering", 0.0), digits=1))s |")
        println(io, "| **Total** | **$(n_filtered) localizations** | **$(round(sum(values(timings)), digits=1))s** |")
    end
end
