"""
Standardized SMLM Analysis Workflow

Results saved to {source_dir}/Results_JuliaSMLM/{basename}/{model}/

Supports:
- 1-channel files: fit right side with both sigma and sxsy models
- 2-channel files: fit whole image with sxsy, analyze left/right polarizations
"""

using SMLMAnalysis
using SMLMDriftCorrection
using CairoMakie
using Statistics
using Printf
using Dates
using FileIO
using Images: RGB

"""
    estimate_mode(x; nbins=100)

Estimate the mode of a distribution using histogram peak finding.
More robust than median for skewed distributions (e.g., PSF sigma).
"""
function estimate_mode(x; nbins=100)
    lo, hi = extrema(x)
    edges = range(lo, hi, length=nbins+1)
    bin_width = (hi - lo) / nbins
    counts = zeros(Int, nbins)
    for v in x
        bin = clamp(floor(Int, (v - lo) / bin_width) + 1, 1, nbins)
        counts[bin] += 1
    end
    max_bin = argmax(counts)
    return (edges[max_bin] + edges[max_bin+1]) / 2
end

# =============================================================================
# INPUT FILE - from command line or default
# =============================================================================
if length(ARGS) >= 1
    h5file = ARGS[1]
else
    # Default file for interactive use
    h5file = "../data/gatta_ruler/2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5"
end

# =============================================================================
# AUTO-DERIVE OUTPUT PATH AND DETECT FILE TYPE
# =============================================================================
h5_abspath = abspath(h5file)
h5_dir = dirname(h5_abspath)
h5_basename = basename(h5_abspath)
h5_name = replace(h5_basename, ".h5" => "")

# Create Results directory alongside source file (MATLAB SMITE style)
results_base = joinpath(h5_dir, "Results_JuliaSMLM")
mkpath(results_base)

# Detect channel configuration from filename
is_2channel = occursin("-2ch", h5_name)
is_1channel = occursin("-1ch", h5_name)

if !is_1channel && !is_2channel
    @warn "Could not detect channel configuration from filename: $h5_name. Assuming 1-channel."
    is_1channel = true
end

# Configure based on channel type
if is_1channel
    # 1-channel: fit right side only (image splitter setup)
    config_roi_x_pixels = (430, 860)
    config_roi_y_pixels = nothing
    # Run both models for 1ch
    run_modes = [:sigma, :sxsy]
    println("Detected: 1-channel file -> right side ROI, both models")
else  # 2-channel
    # 2-channel: fit whole image, analyze polarizations after drift
    config_roi_x_pixels = nothing
    config_roi_y_pixels = nothing
    # Only sxsy model for 2ch (needed for polarization analysis)
    run_modes = [:sxsy]
    println("Detected: 2-channel file -> full image, sxsy model only")
end

# =============================================================================
# FIXED ANALYSIS PARAMETERS
# =============================================================================
const PARAMS = (
    # Data loading
    max_frames = nothing,  # Use all frames

    # Camera (ORCA-Fusion with 60x 1.2NA objective)
    pixel_size_um = 0.078f0,  # 78nm effective pixel size
    readnoise = 0.7f0,        # e- RMS
    offset_adu = 100.0f0,
    gain = 0.24f0,            # e-/ADU
    qe = 0.80f0,

    # Detection (SMLMBoxer) - PSF-aware interface
    boxsize = 9,
    overlap = 2.0,
    psf_sigma = 0.13f0,       # PSF sigma in microns (~130nm typical TIRF)
    min_photons = 400.0f0,    # Detection threshold

    # Fitting (GaussMLE) - model set per run_mode
    iterations = 20,

    # Filtering
    pvalue_threshold = 1e-3,
    max_sigma_xy_nm = 10.0,      # Max precision (sigma_x or sigma_y) in nm
    psf_sigma_tolerance = 0.20,  # ±20% of median PSF sigma

    # Spatial ROI - set based on channel type
    roi_x_pixels = config_roi_x_pixels,
    roi_y_pixels = config_roi_y_pixels,

    # Drift correction
    drift_degree = 3,              # Reduced from 5 for speed
    drift_model = "LegendrePoly",  # LegendrePoly or Polynomial
    drift_cost_fun = "Kdtree",     # Fast (154s vs 3520s for Entropy)

    # Rendering
    zoom = 20,
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

"""Compute empirical CDF values"""
function ecdf_values(x)
    sorted = sort(x)
    n = length(sorted)
    return sorted, (1:n) ./ n
end

"""Convert pixel ROI to micron ROI"""
function pixels_to_microns(pixel_range, pixel_size_um)
    if pixel_range === nothing
        return nothing
    end
    p1, p2 = pixel_range
    return (p1 * pixel_size_um, p2 * pixel_size_um)
end

"""Filter SMLD by spatial ROI (in microns)"""
function filter_roi(smld; x_range=nothing, y_range=nothing)
    emitters = smld.emitters
    if x_range !== nothing
        x_min, x_max = x_range
        if x_min !== nothing
            emitters = filter(e -> e.x >= x_min, emitters)
        end
        if x_max !== nothing
            emitters = filter(e -> e.x <= x_max, emitters)
        end
    end
    if y_range !== nothing
        y_min, y_max = y_range
        if y_min !== nothing
            emitters = filter(e -> e.y >= y_min, emitters)
        end
        if y_max !== nothing
            emitters = filter(e -> e.y <= y_max, emitters)
        end
    end
    return BasicSMLD(collect(emitters), smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)
end

# =============================================================================
# MAIN WORKFLOW - loops over run_modes for 1ch files
# =============================================================================

for run_mode in run_modes

# Set up output directory for this run
if length(run_modes) > 1
    # Multiple modes: create subdirectory for each
    outdir = joinpath(results_base, h5_name, string(run_mode))
else
    # Single mode: use main directory
    outdir = joinpath(results_base, h5_name)
end
mkpath(outdir)

# Select PSF model based on run_mode
if run_mode == :sigma
    psf_model_type = GaussianXYNBS()
    model_name = "GaussianXYNBS (isotropic σ)"
else  # :sxsy
    psf_model_type = GaussianXYNBSXSY()
    model_name = "GaussianXYNBSXSY (anisotropic σx,σy)"
end

println("="^80)
println("SMLM Analysis Workflow")
println("="^80)
println("Started: $(now())")
println("File: $h5_name")
println("Mode: $run_mode ($model_name)")
println("Output: $outdir/")

# Convert pixel ROI to microns
roi_x_um = pixels_to_microns(PARAMS.roi_x_pixels, PARAMS.pixel_size_um)
roi_y_um = pixels_to_microns(PARAMS.roi_y_pixels, PARAMS.pixel_size_um)
println("\nROI (pixels): x=$(PARAMS.roi_x_pixels), y=$(PARAMS.roi_y_pixels)")
println("ROI (microns): x=$roi_x_um, y=$roi_y_um")

# =============================================================================
# STEP 1: LOAD DATA
# =============================================================================
println("\n" * "-"^60)
println("STEP 1: Load Data")
println("-"^60)

info = load_smart_h5_info(h5file)
println("File: $(basename(h5file))")
println("  Full: $(info.width) x $(info.height) pixels, $(info.nframes) frames")

t_load = @elapsed data, _ = smart_h5_to_array(h5file, max_frames=PARAMS.max_frames, verbose=true)
println("  Loaded: $(size(data,3)) frames in $(round(t_load, digits=1))s")

# Create 2D readnoise array (required for GPU fitting - scalar causes issues)
# SCMOSCamera expects (nx, ny) = (cols, rows), so swap dimensions
readnoise_array = fill(PARAMS.readnoise, size(data,2), size(data,1))
camera = SCMOSCamera(size(data,2), size(data,1), PARAMS.pixel_size_um, readnoise_array;
    offset=PARAMS.offset_adu, gain=PARAMS.gain, qe=PARAMS.qe)

# Raw frames grid
all_pixels = vec(data)
pmin, pmax = quantile(all_pixels, [0.01, 0.99])

fig = Figure(size=(1600, 900))
nframes = size(data, 3)
frame_indices = [1, 2, 3, 4,
    round(Int, nframes*0.25), round(Int, nframes*0.33),
    round(Int, nframes*0.50), round(Int, nframes*0.67),
    nframes-3, nframes-2, nframes-1, nframes]
for (idx, frame_num) in enumerate(frame_indices)
    row = div(idx - 1, 4) + 1
    col = mod(idx - 1, 4) + 1
    local ax = Axis(fig[row, col], title="Frame $frame_num", aspect=DataAspect(), yreversed=true)
    heatmap!(ax, data[:, :, frame_num]', colormap=:grays, colorrange=(pmin, pmax))
    hidedecorations!(ax)
end
save(joinpath(outdir, "01_frames_grid.png"), fig)
println("  Saved: 01_frames_grid.png")

# =============================================================================
# STEP 2: DETECTION
# =============================================================================
println("\n" * "-"^60)
println("STEP 2: Detection")
println("-"^60)

# SMLMBoxer handles GPU memory batching internally (PSF-aware detection)
t_detect = @elapsed roi_batch = getboxes(data, camera;
    boxsize=PARAMS.boxsize, overlap=PARAMS.overlap,
    psf_sigma=PARAMS.psf_sigma, min_photons=PARAMS.min_photons,
    use_gpu=true)
num_rois = size(roi_batch.data, 3)
println("  Detected: $num_rois ROIs in $(round(t_detect, digits=1))s ($(round(num_rois/nframes, digits=1))/frame)")

# Detection overlay
fig = Figure(size=(1600, 900))
for (idx, frame_num) in enumerate(frame_indices)
    row = div(idx - 1, 4) + 1
    col = mod(idx - 1, 4) + 1
    local ax = Axis(fig[row, col], title="Frame $frame_num", aspect=DataAspect(), yreversed=true)
    heatmap!(ax, data[:, :, frame_num]', colormap=:grays, colorrange=(pmin, pmax))
    frame_mask = roi_batch.frame_indices .== frame_num
    if any(frame_mask)
        det_x = roi_batch.x_corners[frame_mask]
        det_y = roi_batch.y_corners[frame_mask]
        bs = roi_batch.roi_size
        for (x, y) in zip(det_x, det_y)
            lines!(ax, [x, x+bs, x+bs, x, x], [y, y, y+bs, y+bs, y], color=:yellow, linewidth=0.5)
        end
    end
    hidedecorations!(ax)
end
save(joinpath(outdir, "02_detection_overlay.png"), fig)
println("  Saved: 02_detection_overlay.png")

# =============================================================================
# STEP 3: FITTING
# =============================================================================
println("\n" * "-"^60)
println("STEP 3: Fitting")
println("-"^60)

fitter = GaussMLEFitter(psf_model=psf_model_type, iterations=PARAMS.iterations, device=:gpu)
t_fit = @elapsed smld = fit(fitter, roi_batch)
println("  Fitted $num_rois ROIs in $(round(t_fit, digits=1))s (model: $run_mode)")

# =============================================================================
# STEP 4: SPATIAL ROI FILTERING (before quality diagnostics)
# =============================================================================
println("\n" * "-"^60)
println("STEP 4: Spatial ROI Filtering")
println("-"^60)

# Extract all fit data (before filtering) for comparison plot
x_all = [e.x for e in smld.emitters]
y_all = [e.y for e in smld.emitters]

# Apply spatial ROI filter first (important for image splitter setups)
smld_roi = filter_roi(smld; x_range=roi_x_um, y_range=roi_y_um)
n_roi = length(smld_roi.emitters)
println("  Spatial ROI: $n_roi / $(length(smld.emitters)) ($(round(100*n_roi/length(smld.emitters), digits=1))%)")

# Spatial filter visualization
x_roi = [e.x for e in smld_roi.emitters]
y_roi = [e.y for e in smld_roi.emitters]

fig = Figure(size=(1200, 500))
ax1 = Axis(fig[1, 1], xlabel="X (um)", ylabel="Y (um)", title="All localizations (n=$(length(x_all)))", aspect=DataAspect())
scatter!(ax1, x_all, y_all, markersize=1, alpha=0.1, color=:gray)
if roi_x_um !== nothing
    vlines!(ax1, [roi_x_um[1]], color=:red, linewidth=2, linestyle=:dash)
end
if roi_y_um !== nothing
    hlines!(ax1, [roi_y_um[1]], color=:red, linewidth=2, linestyle=:dash)
end

ax2 = Axis(fig[1, 2], xlabel="X (um)", ylabel="Y (um)", title="After spatial filter (n=$n_roi)", aspect=DataAspect())
scatter!(ax2, x_roi, y_roi, markersize=1, alpha=0.3, color=:blue)

save(joinpath(outdir, "03_spatial_filter.png"), fig)
println("  Saved: 03_spatial_filter.png")

# =============================================================================
# STEP 5: FIT QUALITY DIAGNOSTICS (on spatially filtered data only)
# =============================================================================
println("\n" * "-"^60)
println("STEP 5: Fit Quality Diagnostics (ROI only)")
println("-"^60)

# Extract fit data from ROI-filtered localizations
photons_roi = [e.photons for e in smld_roi.emitters]
bg_roi = [e.bg for e in smld_roi.emitters]
sigma_x_roi = [e.σ_x for e in smld_roi.emitters]  # Position uncertainty (CRLB)
sigma_y_roi = [e.σ_y for e in smld_roi.emitters]  # Position uncertainty (CRLB)
pvalue_roi = [e.pvalue for e in smld_roi.emitters]

# PSF sigma extraction depends on model type
if run_mode == :sigma
    # Isotropic model: single σ field
    sigma_roi = [e.σ for e in smld_roi.emitters]
    psf_sx_roi = sigma_roi  # Same as σ for isotropic
    psf_sy_roi = sigma_roi  # Same as σ for isotropic
else  # :sxsy
    # Anisotropic model: separate σx, σy fields for PSF widths
    psf_sx_roi = [e.σx for e in smld_roi.emitters]  # Fitted PSF width in x
    psf_sy_roi = [e.σy for e in smld_roi.emitters]  # Fitted PSF width in y
    sigma_roi = sqrt.(psf_sx_roi .* psf_sy_roi)  # Geometric mean for filtering
end

# Fit quality diagnostic panel (using ROI data only)
fig = Figure(size=(1800, 1200))

# P-value
ax1 = Axis(fig[1, 1], xlabel="log10(p-value)", ylabel="Density", title="P-value Distribution (ROI)")
ax1_cdf = Axis(fig[1, 1], ylabel="CDF", yaxisposition=:right, yticklabelcolor=:red)
hidespines!(ax1_cdf); hidexdecorations!(ax1_cdf)
pval_nz = pvalue_roi[pvalue_roi .> 0]
if !isempty(pval_nz)
    log_pval = log10.(pval_nz)
    hist!(ax1, log_pval, bins=50, color=(:steelblue, 0.7), normalization=:pdf)
    cx, cy = ecdf_values(log_pval)
    lines!(ax1_cdf, cx, cy, color=:red, linewidth=2)
    vlines!(ax1, [log10(PARAMS.pvalue_threshold)], color=:orange, linestyle=:dash, linewidth=2)
end

# Photons
ax2 = Axis(fig[1, 2], xlabel="Photons", ylabel="Density", title="Photon Distribution (ROI)", xscale=log10)
ax2_cdf = Axis(fig[1, 2], ylabel="CDF", yaxisposition=:right, yticklabelcolor=:red, xscale=log10)
hidespines!(ax2_cdf); hidexdecorations!(ax2_cdf)
ph_pos = photons_roi[photons_roi .> 0]
hist!(ax2, ph_pos, bins=50, color=(:green, 0.7), normalization=:pdf)
cx, cy = ecdf_values(ph_pos)
lines!(ax2_cdf, cx, cy, color=:red, linewidth=2)
vlines!(ax2, [median(ph_pos)], color=:black, linestyle=:dash)

# PSF sigma
ax3 = Axis(fig[1, 3], xlabel="PSF sigma (nm)", ylabel="Density", title="Fitted PSF Sigma (ROI)")
ax3_cdf = Axis(fig[1, 3], ylabel="CDF", yaxisposition=:right, yticklabelcolor=:red)
hidespines!(ax3_cdf); hidexdecorations!(ax3_cdf)
sig_nm = sigma_roi .* 1000
hist!(ax3, sig_nm, bins=50, color=(:purple, 0.7), normalization=:pdf)
cx, cy = ecdf_values(sig_nm)
lines!(ax3_cdf, cx, cy, color=:red, linewidth=2)
vlines!(ax3, [median(sig_nm)], color=:black, linestyle=:dash)

# Background
ax4 = Axis(fig[2, 1], xlabel="Background", ylabel="Density", title="Background Distribution (ROI)", xscale=log10)
ax4_cdf = Axis(fig[2, 1], ylabel="CDF", yaxisposition=:right, yticklabelcolor=:red, xscale=log10)
hidespines!(ax4_cdf); hidexdecorations!(ax4_cdf)
bg_pos = bg_roi[bg_roi .> 0]
hist!(ax4, bg_pos, bins=50, color=(:orange, 0.7), normalization=:pdf)
cx, cy = ecdf_values(bg_pos)
lines!(ax4_cdf, cx, cy, color=:red, linewidth=2)
vlines!(ax4, [median(bg_pos)], color=:black, linestyle=:dash)

# Precision (linear scale, 0-50nm range with 50 bins)
ax5 = Axis(fig[2, 2], xlabel="sigma_x (nm)", ylabel="Density", title="X Precision CRLB (ROI)")
ax5_cdf = Axis(fig[2, 2], ylabel="CDF", yaxisposition=:right, yticklabelcolor=:red)
hidespines!(ax5_cdf); hidexdecorations!(ax5_cdf)
prec_nm = sigma_x_roi .* 1000
prec_pos = prec_nm[(prec_nm .> 0) .& (prec_nm .<= 50)]  # Filter to 0-50nm range
hist!(ax5, prec_pos, bins=range(0, 50, length=51), color=(:teal, 0.7), normalization=:pdf)
xlims!(ax5, 0, 50)
cx, cy = ecdf_values(prec_pos)
lines!(ax5_cdf, cx, cy, color=:red, linewidth=2)
xlims!(ax5_cdf, 0, 50)
vlines!(ax5, [median(prec_pos)], color=:black, linestyle=:dash)

# Photons vs Precision scatter
ax6 = Axis(fig[2, 3], xlabel="Photons", ylabel="sigma_x (nm)", title="Photons vs Precision (ROI)", xscale=log10, yscale=log10)
scatter!(ax6, photons_roi[1:min(5000,end)], sigma_x_roi[1:min(5000,end)] .* 1000, markersize=2, alpha=0.3)

save(joinpath(outdir, "04_fit_quality.png"), fig)
println("  Saved: 04_fit_quality.png")

# Print stats for ROI
println("\n  Fit Statistics (ROI: $n_roi localizations):")
println("    Photons: median=$(round(median(photons_roi), digits=0))")
println("    PSF sigma: median=$(round(median(sigma_roi)*1000, digits=1)) nm")
println("    Precision: median=$(round(median(sigma_x_roi)*1000, digits=1)) nm")
println("    P-value > $(PARAMS.pvalue_threshold): $(round(100*count(p -> p > PARAMS.pvalue_threshold, pvalue_roi)/length(pvalue_roi), digits=1))%")

# PSF σx/σy CDF plot for sxsy model (1ch or 2ch)
if run_mode == :sxsy
    fig_cdf = Figure(size=(600, 400))
    ax_cdf = Axis(fig_cdf[1, 1], xlabel="PSF width (nm)", ylabel="CDF",
                  title="PSF σx/σy CDF (n=$n_roi)")

    # Convert to nm
    sx_nm = psf_sx_roi .* 1000
    sy_nm = psf_sy_roi .* 1000

    # Plot CDFs
    cx, cy = ecdf_values(sx_nm)
    lines!(ax_cdf, cx, cy, color=:blue, linewidth=2, label="σx (med=$(round(median(sx_nm), digits=1)))")
    cx, cy = ecdf_values(sy_nm)
    lines!(ax_cdf, cx, cy, color=:red, linewidth=2, label="σy (med=$(round(median(sy_nm), digits=1)))")

    axislegend(ax_cdf, position=:rb)
    save(joinpath(outdir, "04b_psf_sxsy_cdf.png"), fig_cdf)
    println("  Saved: 04b_psf_sxsy_cdf.png")
end

# =============================================================================
# STEP 6: QUALITY FILTERING
# =============================================================================
println("\n" * "-"^60)
println("STEP 6: Quality Filtering")
println("-"^60)

# Convert nm thresholds to microns
max_sigma_xy_um = PARAMS.max_sigma_xy_nm / 1000.0

# PSF sigma filter: ±20% of mode (more robust than median for skewed distributions)
mode_psf_sigma = estimate_mode(sigma_roi)
psf_sigma_min = mode_psf_sigma * (1 - PARAMS.psf_sigma_tolerance)
psf_sigma_max = mode_psf_sigma * (1 + PARAMS.psf_sigma_tolerance)

# Apply all filters
good_idx = (pvalue_roi .> PARAMS.pvalue_threshold) .&
           (sigma_x_roi .<= max_sigma_xy_um) .&
           (sigma_y_roi .<= max_sigma_xy_um) .&
           (sigma_roi .>= psf_sigma_min) .&
           (sigma_roi .<= psf_sigma_max)

smld_filtered = BasicSMLD(smld_roi.emitters[good_idx], smld_roi.camera, smld_roi.n_frames, smld_roi.n_datasets, smld_roi.metadata)
n_final = length(smld_filtered.emitters)

# Print filter statistics
n_pval = count(pvalue_roi .> PARAMS.pvalue_threshold)
n_sigx = count(sigma_x_roi .<= max_sigma_xy_um)
n_sigy = count(sigma_y_roi .<= max_sigma_xy_um)
n_psf = count((sigma_roi .>= psf_sigma_min) .& (sigma_roi .<= psf_sigma_max))

println("  Filter criteria:")
println("    P-value > $(PARAMS.pvalue_threshold): $n_pval / $n_roi ($(round(100*n_pval/n_roi, digits=1))%)")
println("    sigma_x <= $(PARAMS.max_sigma_xy_nm) nm: $n_sigx / $n_roi ($(round(100*n_sigx/n_roi, digits=1))%)")
println("    sigma_y <= $(PARAMS.max_sigma_xy_nm) nm: $n_sigy / $n_roi ($(round(100*n_sigy/n_roi, digits=1))%)")
println("    PSF sigma: $(round(psf_sigma_min*1000, digits=1))-$(round(psf_sigma_max*1000, digits=1)) nm (±$(Int(PARAMS.psf_sigma_tolerance*100))% of mode $(round(mode_psf_sigma*1000, digits=1)) nm): $n_psf / $n_roi ($(round(100*n_psf/n_roi, digits=1))%)")
println("  Combined: $n_final / $n_roi ($(round(100*n_final/n_roi, digits=1))%)")

# Plot accepted fits per frame
frames_filtered = [e.frame for e in smld_filtered.emitters]
frame_counts = Dict{Int, Int}()
for f in frames_filtered
    frame_counts[f] = get(frame_counts, f, 0) + 1
end
all_frames = 1:maximum(frames_filtered)
counts_per_frame = [get(frame_counts, f, 0) for f in all_frames]

fig = Figure(size=(800, 400))
ax = Axis(fig[1, 1], xlabel="Frame", ylabel="Accepted Fits", title="Accepted Fits per Frame (after p-value filter)")
barplot!(ax, collect(all_frames), counts_per_frame, color=:steelblue)
hlines!(ax, [mean(counts_per_frame)], color=:red, linestyle=:dash, linewidth=2, label="Mean: $(round(mean(counts_per_frame), digits=1))")
axislegend(ax, position=:rt)
save(joinpath(outdir, "07_fits_per_frame.png"), fig)
println("  Saved: 07_fits_per_frame.png")

# =============================================================================
# STEP 7: DRIFT CORRECTION
# =============================================================================
println("\n" * "-"^60)
println("STEP 7: Drift Correction")
println("-"^60)

# Store original coordinates for drift visualization
x_orig = [e.x for e in smld_filtered.emitters]
y_orig = [e.y for e in smld_filtered.emitters]
frames_orig = [e.frame for e in smld_filtered.emitters]

# Debug: print input coordinate ranges
println("  Input coordinates: X=$(round(minimum(x_orig)*1000, digits=1))-$(round(maximum(x_orig)*1000, digits=1))nm, Y=$(round(minimum(y_orig)*1000, digits=1))-$(round(maximum(y_orig)*1000, digits=1))nm")

# Drift correction with PARAMS settings
t_drift = @elapsed begin
    smld_corrected = driftcorrect(smld_filtered;
                                  degree=PARAMS.drift_degree,
                                  cost_fun=PARAMS.drift_cost_fun,
                                  intramodel=PARAMS.drift_model, verbose=1)
end
println("  Drift correction completed in $(round(t_drift, digits=1))s")
println("  Localizations: $(length(smld_corrected.emitters))")

# Debug: print output coordinate ranges
x_corr_dbg = [e.x for e in smld_corrected.emitters]
y_corr_dbg = [e.y for e in smld_corrected.emitters]
println("  Output coordinates: X=$(round(minimum(x_corr_dbg)*1000, digits=1))-$(round(maximum(x_corr_dbg)*1000, digits=1))nm, Y=$(round(minimum(y_corr_dbg)*1000, digits=1))-$(round(maximum(y_corr_dbg)*1000, digits=1))nm")

# Compute drift trajectory (difference between original and corrected)
x_corr = [e.x for e in smld_corrected.emitters]
y_corr = [e.y for e in smld_corrected.emitters]
dx = x_orig .- x_corr  # drift applied = original - corrected
dy = y_orig .- y_corr

# Bin drift by frame to get trajectory
unique_frames = sort(unique(frames_orig))
drift_x = Float64[]
drift_y = Float64[]
for f in unique_frames
    mask = frames_orig .== f
    push!(drift_x, mean(dx[mask]))
    push!(drift_y, mean(dy[mask]))
end

# Plot drift trajectory
fig = Figure(size=(1400, 400))
ax1 = Axis(fig[1, 1], xlabel="Frame", ylabel="X Drift (nm)", title="X Drift vs Frame")
lines!(ax1, unique_frames, drift_x .* 1000, color=:blue)
ax2 = Axis(fig[1, 2], xlabel="Frame", ylabel="Y Drift (nm)", title="Y Drift vs Frame")
lines!(ax2, unique_frames, drift_y .* 1000, color=:red)
# Equal aspect ratio for XY path
ax3 = Axis(fig[1, 3], xlabel="X Drift (nm)", ylabel="Y Drift (nm)", title="XY Drift Path", aspect=DataAspect())
lines!(ax3, drift_x .* 1000, drift_y .* 1000, color=:black, linewidth=1.5)
scatter!(ax3, [drift_x[1]*1000], [drift_y[1]*1000], color=:green, markersize=12, label="Start")
scatter!(ax3, [drift_x[end]*1000], [drift_y[end]*1000], color=:red, markersize=12, label="End")
axislegend(ax3, position=:lt)
save(joinpath(outdir, "06_drift_trajectory.png"), fig)
println("  Saved: 06_drift_trajectory.png")
println("  Max drift: X=$(round(maximum(abs.(drift_x))*1000, digits=1))nm, Y=$(round(maximum(abs.(drift_y))*1000, digits=1))nm")

# Use corrected data for rendering
smld_filtered = smld_corrected

# Save drift-corrected SMLD to HDF5
h5_output_file = joinpath(outdir, "smld_corrected.h5")
save_smld(h5_output_file, smld_corrected; source_file=h5_abspath)
println("  Saved: smld_corrected.h5")

# Extract final filtered data
x_filt = [e.x for e in smld_filtered.emitters]
y_filt = [e.y for e in smld_filtered.emitters]
photons_filt = [e.photons for e in smld_filtered.emitters]
sigma_x_filt = [e.σ_x for e in smld_filtered.emitters]

# =============================================================================
# STEP 8: SUPER-RESOLUTION RENDERING (restricted to ROI)
# =============================================================================
println("\n" * "-"^60)
println("STEP 8: Super-Resolution Rendering")
println("-"^60)

# Build roi tuple from PARAMS - use camera pixel ranges directly
# SMLMVis accepts: nothing (full), (range, range), (range, :), (:, range) but NOT (:, :)
render_roi = if PARAMS.roi_x_pixels === nothing && PARAMS.roi_y_pixels === nothing
    nothing  # Full image
elseif PARAMS.roi_x_pixels !== nothing && PARAMS.roi_y_pixels !== nothing
    (PARAMS.roi_x_pixels[1]:PARAMS.roi_x_pixels[2], PARAMS.roi_y_pixels[1]:PARAMS.roi_y_pixels[2])
elseif PARAMS.roi_x_pixels !== nothing
    (PARAMS.roi_x_pixels[1]:PARAMS.roi_x_pixels[2], :)
else
    (:, PARAMS.roi_y_pixels[1]:PARAMS.roi_y_pixels[2])
end
println("  ROI (camera pixels): x=$(PARAMS.roi_x_pixels), y=$(PARAMS.roi_y_pixels)")

# Zoom levels - reduced circle_zoom from 100 to 50 to avoid OOM on large images
circle_zoom = 50
pixel_size_gauss = PARAMS.pixel_size_um / PARAMS.zoom
pixel_size_circle = PARAMS.pixel_size_um / circle_zoom

"""Add scale bar to an image file (bottom-left corner)"""
function add_scalebar!(filename; scalebar_um=1.0, pixel_size_um=0.005, bar_color=:white, max_pixels=50_000_000)
    # Check file size first to avoid loading huge images
    file_size_mb = filesize(filename) / 1e6

    # Try to get dimensions without loading full image
    # PNG dimensions are in header - use ImageMagick identify if available
    dims_result = try
        dims_str = read(`identify -format "%w %h" $filename`, String)
        parts = split(strip(dims_str))
        (parse(Int, parts[1]), parse(Int, parts[2]))  # (width, height)
    catch
        nothing
    end

    if dims_result !== nothing
        w, h = dims_result
        total_pixels = w * h
        if total_pixels > max_pixels
            println("    Skipping scalebar for large image ($(w)x$(h) = $(round(total_pixels/1e6, digits=1))M pixels)")
            return 0
        end
    elseif file_size_mb > 10  # Fallback: skip if file > 10MB
        println("    Skipping scalebar for large file ($(round(file_size_mb, digits=1)) MB)")
        return 0
    end

    # Safe to load
    img = load(filename)
    h, w = size(img)

    # Scale bar dimensions
    bar_length_px = round(Int, scalebar_um / pixel_size_um)
    bar_height_px = max(5, round(Int, h * 0.01))
    margin_px = round(Int, min(h, w) * 0.03)

    # Bar position (bottom-left)
    bar_x = margin_px
    bar_y = h - margin_px - bar_height_px

    # Draw scale bar
    for y in bar_y:(bar_y + bar_height_px)
        for x in bar_x:(bar_x + bar_length_px)
            if 1 <= y <= h && 1 <= x <= w
                img[y, x] = bar_color == :white ? RGB(1.0, 1.0, 1.0) : RGB(0.0, 0.0, 0.0)
            end
        end
    end

    save(filename, img)
    return bar_length_px
end

# Gaussian render
println("  Rendering Gaussian (inferno) at $(PARAMS.zoom)x...")
t_render = @elapsed result = render(smld_filtered;
    strategy=GaussianRender(),
    zoom=PARAMS.zoom,
    roi=render_roi,
    colormap=:inferno,
    clip_percentile=0.9999,
    filename=joinpath(outdir, "05_gaussian_inferno.png"))
add_scalebar!(joinpath(outdir, "05_gaussian_inferno.png"); scalebar_um=1.0, pixel_size_um=pixel_size_gauss)
println("    Saved: 05_gaussian_inferno.png + 1um scalebar ($(round(t_render, digits=1))s, $(size(result.image)))")

# Gaussian render with 99% percentile clip (more contrast)
println("  Rendering Gaussian (inferno, 99% clip) at $(PARAMS.zoom)x...")
t_render = @elapsed result = render(smld_filtered;
    strategy=GaussianRender(),
    zoom=PARAMS.zoom,
    roi=render_roi,
    colormap=:inferno,
    clip_percentile=0.99,
    filename=joinpath(outdir, "05_gaussian_inferno_99p.png"))
add_scalebar!(joinpath(outdir, "05_gaussian_inferno_99p.png"); scalebar_um=1.0, pixel_size_um=pixel_size_gauss)
println("    Saved: 05_gaussian_inferno_99p.png + 1um scalebar ($(round(t_render, digits=1))s, $(size(result.image)))")

# Histogram colored by frame (turbo) at 20x zoom
println("  Rendering Histogram by frame (turbo) at $(PARAMS.zoom)x...")
t_render = @elapsed result = render(smld_filtered;
    strategy=HistogramRender(),
    zoom=PARAMS.zoom,
    roi=render_roi,
    color_by=:frame,
    colormap=:turbo,
    filename=joinpath(outdir, "05_histogram_frame.png"))
add_scalebar!(joinpath(outdir, "05_histogram_frame.png"); scalebar_um=1.0, pixel_size_um=pixel_size_gauss)
println("    Saved: 05_histogram_frame.png + 1um scalebar ($(round(t_render, digits=1))s, $(size(result.image)))")

# Circles colored by frame (turbo) at 100x zoom
println("  Rendering Circles by frame (turbo) at $(circle_zoom)x...")
t_render = @elapsed result = render(smld_filtered;
    strategy=CircleRender(),
    zoom=circle_zoom,
    roi=render_roi,
    color_by=:frame,
    colormap=:turbo,
    filename=joinpath(outdir, "05_circles_frame.png"))
add_scalebar!(joinpath(outdir, "05_circles_frame.png"); scalebar_um=1.0, pixel_size_um=pixel_size_circle)
println("    Saved: 05_circles_frame.png + 1um scalebar ($(round(t_render, digits=1))s, $(size(result.image)))")

# =============================================================================
# STEP 9: POLARIZATION ANALYSIS (2-channel only)
# =============================================================================
if is_2channel && run_mode == :sxsy
    println("\n" * "-"^60)
    println("STEP 9: Polarization Analysis (2-channel)")
    println("-"^60)

    # Split at image center (x = 43.0 µm = pixel 430)
    polarization_split_um = 43.0

    # Extract PSF widths from corrected/filtered data
    psf_sx_final = [e.σx for e in smld_filtered.emitters]
    psf_sy_final = [e.σy for e in smld_filtered.emitters]
    x_final = [e.x for e in smld_filtered.emitters]

    # Split by polarization (left = pol1, right = pol2)
    left_mask = x_final .< polarization_split_um
    right_mask = x_final .>= polarization_split_um

    # Convert to nm for plotting
    sx_left_nm = psf_sx_final[left_mask] .* 1000
    sy_left_nm = psf_sy_final[left_mask] .* 1000
    sx_right_nm = psf_sx_final[right_mask] .* 1000
    sy_right_nm = psf_sy_final[right_mask] .* 1000

    n_left = sum(left_mask)
    n_right = sum(right_mask)
    println("  Localizations: Left (pol1)=$n_left, Right (pol2)=$n_right")

    # Create 4-panel CDF plot
    fig = Figure(size=(1200, 800))

    # Plot 1: σx Left (Polarization 1)
    ax1 = Axis(fig[1, 1], xlabel="PSF σx (nm)", ylabel="CDF",
               title="σx Left - Polarization 1 (n=$n_left)")
    if !isempty(sx_left_nm)
        cx, cy = ecdf_values(sx_left_nm)
        lines!(ax1, cx, cy, color=:blue, linewidth=2)
        vlines!(ax1, [median(sx_left_nm)], color=:red, linestyle=:dash, linewidth=1.5,
                label="median=$(round(median(sx_left_nm), digits=1)) nm")
        axislegend(ax1, position=:rb)
    end

    # Plot 2: σy Left (Polarization 1)
    ax2 = Axis(fig[1, 2], xlabel="PSF σy (nm)", ylabel="CDF",
               title="σy Left - Polarization 1 (n=$n_left)")
    if !isempty(sy_left_nm)
        cx, cy = ecdf_values(sy_left_nm)
        lines!(ax2, cx, cy, color=:green, linewidth=2)
        vlines!(ax2, [median(sy_left_nm)], color=:red, linestyle=:dash, linewidth=1.5,
                label="median=$(round(median(sy_left_nm), digits=1)) nm")
        axislegend(ax2, position=:rb)
    end

    # Plot 3: σx Right (Polarization 2)
    ax3 = Axis(fig[2, 1], xlabel="PSF σx (nm)", ylabel="CDF",
               title="σx Right - Polarization 2 (n=$n_right)")
    if !isempty(sx_right_nm)
        cx, cy = ecdf_values(sx_right_nm)
        lines!(ax3, cx, cy, color=:blue, linewidth=2)
        vlines!(ax3, [median(sx_right_nm)], color=:red, linestyle=:dash, linewidth=1.5,
                label="median=$(round(median(sx_right_nm), digits=1)) nm")
        axislegend(ax3, position=:rb)
    end

    # Plot 4: σy Right (Polarization 2)
    ax4 = Axis(fig[2, 2], xlabel="PSF σy (nm)", ylabel="CDF",
               title="σy Right - Polarization 2 (n=$n_right)")
    if !isempty(sy_right_nm)
        cx, cy = ecdf_values(sy_right_nm)
        lines!(ax4, cx, cy, color=:green, linewidth=2)
        vlines!(ax4, [median(sy_right_nm)], color=:red, linestyle=:dash, linewidth=1.5,
                label="median=$(round(median(sy_right_nm), digits=1)) nm")
        axislegend(ax4, position=:rb)
    end

    # Add overall title
    Label(fig[0, :], "PSF Width Comparison: Polarization 1 (Left) vs Polarization 2 (Right)",
          fontsize=16, font=:bold)

    save(joinpath(outdir, "08_polarization_cdf.png"), fig)
    println("  Saved: 08_polarization_cdf.png")

    # Print polarization summary
    println("\n  Polarization PSF Summary:")
    println("    Left (Pol1):  σx=$(round(median(sx_left_nm), digits=1))nm, σy=$(round(median(sy_left_nm), digits=1))nm")
    println("    Right (Pol2): σx=$(round(median(sx_right_nm), digits=1))nm, σy=$(round(median(sy_right_nm), digits=1))nm")
    println("    Difference:   Δσx=$(round(median(sx_right_nm) - median(sx_left_nm), digits=1))nm, Δσy=$(round(median(sy_right_nm) - median(sy_left_nm), digits=1))nm")

    # Overlay CDF comparison plot
    fig2 = Figure(size=(1000, 500))

    ax_sx = Axis(fig2[1, 1], xlabel="PSF σx (nm)", ylabel="CDF",
                 title="σx Comparison: Pol1 vs Pol2")
    if !isempty(sx_left_nm) && !isempty(sx_right_nm)
        cx1, cy1 = ecdf_values(sx_left_nm)
        cx2, cy2 = ecdf_values(sx_right_nm)
        lines!(ax_sx, cx1, cy1, color=:blue, linewidth=2, label="Left (Pol1)")
        lines!(ax_sx, cx2, cy2, color=:orange, linewidth=2, label="Right (Pol2)")
        axislegend(ax_sx, position=:rb)
    end

    ax_sy = Axis(fig2[1, 2], xlabel="PSF σy (nm)", ylabel="CDF",
                 title="σy Comparison: Pol1 vs Pol2")
    if !isempty(sy_left_nm) && !isempty(sy_right_nm)
        cx1, cy1 = ecdf_values(sy_left_nm)
        cx2, cy2 = ecdf_values(sy_right_nm)
        lines!(ax_sy, cx1, cy1, color=:blue, linewidth=2, label="Left (Pol1)")
        lines!(ax_sy, cx2, cy2, color=:orange, linewidth=2, label="Right (Pol2)")
        axislegend(ax_sy, position=:rb)
    end

    save(joinpath(outdir, "08_polarization_comparison.png"), fig2)
    println("  Saved: 08_polarization_comparison.png")
end

# =============================================================================
# SUMMARY
# =============================================================================
println("\n" * "="^80)
println("WORKFLOW COMPLETE")
println("="^80)
println("Finished: $(now())")
println()
println("Results in $outdir/:")
for f in sort(readdir(outdir))
    if endswith(f, ".png") || endswith(f, ".txt")
        println("  $f")
    end
end
println()
println("Summary:")
println("  Input: $(size(data,3)) frames")
println("  Detected: $num_rois ROIs")
println("  Final: $n_final localizations")
println("  Median precision: $(round(median(sigma_x_filt)*1000, digits=1)) nm")
println("  Median photons: $(round(median(photons_filt), digits=0))")

# Save stats to file
open(joinpath(outdir, "analysis_stats.txt"), "w") do io
    println(io, "SMLM Analysis Statistics")
    println(io, "="^60)
    println(io, "Date: $(now())")
    println(io, "File: $(basename(h5file))")
    println(io)
    println(io, "Camera:")
    println(io, "  Pixel size: $(PARAMS.pixel_size_um) um")
    println(io, "  Gain: $(PARAMS.gain) e-/ADU")
    println(io, "  Offset: $(PARAMS.offset_adu) ADU")
    println(io)
    println(io, "ROI (pixels): x=$(PARAMS.roi_x_pixels), y=$(PARAMS.roi_y_pixels)")
    println(io, "ROI (microns): x=$roi_x_um, y=$roi_y_um")
    println(io)
    println(io, "Data Loading: ($(round(t_load, digits=1))s)")
    println(io, "  Frames: $(size(data,3))")
    println(io)
    println(io, "Detection: ($(round(t_detect, digits=1))s)")
    println(io, "  Box size: $(PARAMS.boxsize)")
    println(io, "  PSF sigma: $(PARAMS.psf_sigma) μm")
    println(io, "  Min photons: $(PARAMS.min_photons)")
    println(io, "  Detected: $num_rois ROIs")
    println(io)
    println(io, "Fitting: ($(round(t_fit, digits=1))s)")
    println(io, "  Model: $model_name")
    println(io, "  Iterations: $(PARAMS.iterations)")
    println(io)
    println(io, "Filtering (spatial first, then quality):")
    println(io, "  After spatial ROI: $n_roi ($(round(100*n_roi/length(smld.emitters), digits=1))% of all fits)")
    println(io, "  P-value threshold: $(PARAMS.pvalue_threshold)")
    println(io, "  After p-value: $n_final ($(round(100*n_final/n_roi, digits=1))% of ROI)")
    println(io)
    println(io, "Drift Correction: ($(round(t_drift, digits=1))s)")
    println(io, "  Model: LegendrePoly (degree $(PARAMS.drift_degree))")
    println(io, "  Cost function: $(PARAMS.drift_cost_fun)")
    println(io)
    println(io, "Final Results:")
    println(io, "  Localizations: $n_final")
    println(io, "  X range: $(round(minimum(x_filt), digits=2)) - $(round(maximum(x_filt), digits=2)) um")
    println(io, "  Y range: $(round(minimum(y_filt), digits=2)) - $(round(maximum(y_filt), digits=2)) um")
    println(io, "  Photons (median): $(round(median(photons_filt), digits=0))")
    println(io, "  Precision (median): $(round(median(sigma_x_filt)*1000, digits=1)) nm")
end
println("  Saved: analysis_stats.txt")

end  # for run_mode in run_modes

println("\n" * "="^80)
println("ALL WORKFLOWS COMPLETE")
println("="^80)
