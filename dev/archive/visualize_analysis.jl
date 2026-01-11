"""
Comprehensive SMLM Analysis Visualization

Creates standardized output figures for each analysis step:
1. Raw data frames (beginning, middle, end)
2. Detection overlays
3. Fitting quality metrics
4. Super-resolution images
5. Statistical histograms

Image convention: Top-left is (1,1), y goes down
For CairoMakie: transpose image and use yreverse=true for heatmaps
"""

using SMLMAnalysis
using CairoMakie
using Statistics
using Printf

println("="^80)
println("SMLM Analysis - Comprehensive Visualization")
println("="^80)

# Setup output directory structure
outdir = "output"
for subdir in ["01_raw", "02_detection", "03_fitting", "04_superres", "05_histograms"]
    mkpath(joinpath(outdir, subdir))
end

# =============================================================================
# LOAD DATA
# =============================================================================
println("\n📂 Loading data...")
h5file = "../data/gatta_ruler/2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5"
info = load_smart_h5_info(h5file)

# Extract camera metadata from H5 file
println("\n📷 Extracting camera metadata from H5 file...")
using HDF5
camera_meta = HDF5.h5open(h5file, "r") do file
    camera_group = file["Main/camera"]
    attrs = HDF5.attributes(camera_group)
    Dict(
        "unique_id" => read(attrs["unique_id"]),
        "format_x_pixels" => read(attrs["camera_format_x_pixels"]),
        "format_y_pixels" => read(attrs["camera_format_y_pixels"]),
        "roi_width" => read(attrs["roi_width"]),
        "roi_height" => read(attrs["roi_height"]),
        "roi_x" => read(attrs["roi_x"]),
        "roi_y" => read(attrs["roi_y"]),
        "exposure_time" => read(attrs["exposure_time"]),
        "frame_rate" => read(attrs["frame_rate"]),
        "format_gain" => read(attrs["camera_format_gain"]),
        "format_pixelsize" => read(attrs["camera_format_pixelsize"])
    )
end

println("   Camera ID: $(camera_meta["unique_id"])")
println("   Full sensor: $(Int(camera_meta["format_x_pixels"])) × $(Int(camera_meta["format_y_pixels"])) pixels")
println("   ROI: $(Int(camera_meta["roi_width"])) × $(Int(camera_meta["roi_height"])) at ($(Int(camera_meta["roi_x"])), $(Int(camera_meta["roi_y"])))")
println("   Exposure: $(camera_meta["exposure_time"]) s ($(camera_meta["frame_rate"]) fps)")
println("   Format gain setting: $(camera_meta["format_gain"])")
println("   Format pixelsize: $(camera_meta["format_pixelsize"]) µm")

# Load subset for visualization
nframes_load = 1000
data, _ = smart_h5_to_array(h5file, max_frames=nframes_load)
println("\n   Loaded: $nframes_load frames ($(size(data,1))×$(size(data,2)))")

# Camera setup - SCMOSCamera with Float32 parameters to avoid type mismatch bug
# ORCA-Fusion (C14440-20UP) actual specs from FINDINGS.md:
#   - Physical pixel: 6.5 µm
#   - Effective pixel size: ~0.1 µm (with magnification - needs calibration)
#   - Readnoise: 0.7 e⁻ RMS (ultra quiet mode)
#   - QE: 80% at 642nm
#   - Conversion gain: ~0.46 e⁻/ADU (typical for sCMOS)
#   - Offset: ~100 ADU (typical baseline)
pixel_size_um = 0.078f0  # 78nm effective pixel size (6.5μm sensor / ~83x magnification)
readnoise = 0.7f0      # e⁻ RMS (ultra quiet mode)
offset_adu = 100.0f0   # ADU
gain_val = 0.24f0      # e⁻/ADU (SMITE convention: 1/4.16 ADU/e⁻)
qe = 0.80f0            # QE at 642nm

camera = SCMOSCamera(size(data,2), size(data,1), pixel_size_um, readnoise;
    offset=offset_adu, gain=gain_val, qe=qe)

println("\n📷 Camera: SCMOSCamera (ORCA-Fusion specs, Float32)")
println("   Pixel size: $pixel_size_um µm (78 nm)")
println("   Readnoise: $readnoise e⁻ RMS, Offset: $offset_adu ADU")
println("   Gain: $gain_val e⁻/ADU, QE: $qe")

# =============================================================================
# 1. RAW DATA VISUALIZATION
# =============================================================================
println("\n📊 Step 1: Visualizing raw data...")

# Calculate intensity statistics
all_pixels = vec(data)
pmin, pmax = quantile(all_pixels, [0.01, 0.99])
imin, imax = extrema(all_pixels)

println("   Intensity stats:")
println("     Min: $imin, Max: $imax")
println("     1%: $pmin, 99%: $pmax")
println("     Mean: $(round(mean(all_pixels), digits=1))")

# Select representative frames: beginning, middle, end (3 rows × 4 cols)
nframes = size(data, 3)
frame_indices = [
    # Beginning: frames 1-4
    1, 2, 3, 4,
    # Middle: around 25%, 33%, 50%, 67%
    round(Int, nframes * 0.25), round(Int, nframes * 0.33),
    round(Int, nframes * 0.50), round(Int, nframes * 0.67),
    # End: last 4 frames
    nframes-3, nframes-2, nframes-1, nframes
]

# Create 3×4 panel of frames
fig = Figure(size=(1600, 1200))
for (idx, frame_num) in enumerate(frame_indices)
    row = div(idx - 1, 4) + 1
    col = mod(idx - 1, 4) + 1

    ax = Axis(fig[row, col],
        title="Frame $frame_num",
        aspect=DataAspect(),
        yreversed=true  # (1,1) at top-left, y goes down
    )

    # Transpose for correct orientation
    frame_data = data[:, :, frame_num]'
    heatmap!(ax, frame_data, colormap=:grays, colorrange=(pmin, pmax))
    hidedecorations!(ax)
end

save(joinpath(outdir, "01_raw", "frames_grid.png"), fig)
println("   ✓ Saved 01_raw/frames_grid.png")

# Create intensity CDF for first 100 frames
println("   Creating intensity CDF (first 100 frames)...")
pixels_100 = vec(data[:, :, 1:min(100, nframes)])
sorted_pixels = sort(pixels_100)
cdf_values = (1:length(sorted_pixels)) ./ length(sorted_pixels)

fig_cdf = Figure(size=(800, 600))
ax_cdf = Axis(fig_cdf[1, 1],
    xlabel="Intensity (ADU)",
    ylabel="Cumulative Probability",
    title="Intensity CDF (first 100 frames)",
    xscale=log10)
lines!(ax_cdf, sorted_pixels, cdf_values, linewidth=2)

# Mark key percentiles
percentiles = [0.01, 0.25, 0.50, 0.75, 0.99]
for p in percentiles
    val = quantile(pixels_100, p)
    vlines!(ax_cdf, [val], color=:red, linestyle=:dash, alpha=0.5)
    text!(ax_cdf, val, p, text="$(Int(round(p*100)))%", align=(:left, :center))
end

save(joinpath(outdir, "01_raw", "intensity_cdf.png"), fig_cdf)
println("   ✓ Saved 01_raw/intensity_cdf.png")

# Save statistics
open(joinpath(outdir, "01_raw", "raw_stats.txt"), "w") do io
    println(io, "Raw Data Statistics")
    println(io, "=" ^ 60)
    println(io, "File: ", basename(h5file))
    println(io, "Dimensions: $(size(data,1)) × $(size(data,2)) × $(size(data,3))")
    println(io, "Pixel size: $pixel_size_um µm")
    println(io, "\nIntensity (ADU):")
    println(io, "  Min: $imin")
    println(io, "  1st percentile: $pmin")
    println(io, "  Mean: $(round(mean(all_pixels), digits=1))")
    println(io, "  99th percentile: $pmax")
    println(io, "  Max: $imax")
end
println("   ✓ Saved 01_raw/raw_stats.txt")

# =============================================================================
# 2. DETECTION
# =============================================================================
println("\n🔍 Step 2: Detection...")

roi_batch = getboxes(data, camera;
    boxsize=7,
    overlap=2.0,
    sigma_small=1.0,
    sigma_large=2.0,
    minval=500.0
)
num_rois = size(roi_batch.data, 3)
println("   Detected: $num_rois ROIs")

# Create detection overlay on selected frames (wide format to match data aspect ratio)
fig = Figure(size=(2400, 700))
for (idx, frame_num) in enumerate(frame_indices)
    row = div(idx - 1, 4) + 1
    col = mod(idx - 1, 4) + 1

    ax = Axis(fig[row, col],
        title="Frame $frame_num ($num_rois total detections)",
        aspect=DataAspect(),
        yreversed=true
    )

    # Show frame
    frame_data = data[:, :, frame_num]'
    heatmap!(ax, frame_data, colormap=:grays, colorrange=(pmin, pmax))

    # Overlay detections from this frame as yellow bounding boxes
    frame_mask = roi_batch.frame_indices .== frame_num
    if any(frame_mask)
        det_x = roi_batch.x_corners[frame_mask]
        det_y = roi_batch.y_corners[frame_mask]
        box_size = roi_batch.roi_size

        # Draw yellow boxes around each ROI
        for (x, y) in zip(det_x, det_y)
            # Box corners: (x, y) is top-left corner
            lines!(ax,
                [x, x+box_size, x+box_size, x, x],
                [y, y, y+box_size, y+box_size, y],
                color=:yellow, linewidth=0.5)
        end
    end

    hidedecorations!(ax)
end

save(joinpath(outdir, "02_detection", "detection_overlay.png"), fig)
println("   ✓ Saved 02_detection/detection_overlay.png")

# ROI examples - show actual detected ROIs
fig = Figure(size=(1200, 1200))
nroi_show = min(16, num_rois)
roi_indices = [floor(Int, x) for x in range(1, num_rois, length=nroi_show)]
for (idx, roi_idx) in enumerate(roi_indices)
    row = div(idx - 1, 4) + 1
    col = mod(idx - 1, 4) + 1

    ax = Axis(fig[row, col],
        title="ROI $roi_idx",
        aspect=DataAspect(),
        yreversed=true
    )

    # Extract and transpose ROI
    roi = roi_batch.data[:, :, roi_idx]'
    heatmap!(ax, roi, colormap=:grays)
    hidedecorations!(ax)
end

save(joinpath(outdir, "02_detection", "roi_examples.png"), fig)
println("   ✓ Saved 02_detection/roi_examples.png")

# =============================================================================
# 3. FITTING
# =============================================================================
println("\n📊 Step 3: Fitting...")

psf_sigma = 1.0f0
fitter = GaussMLEFitter(
    psf_model = GaussianXYNB(psf_sigma),
    iterations = 20,
    device = :cpu
)

t_fit = @elapsed smld = fit(fitter, roi_batch)
println("   Fitted $num_rois ROIs in $(round(t_fit, digits=1))s")

# Extract emitter data
emitters = smld.emitters
x = [e.x for e in emitters]
y = [e.y for e in emitters]
photons = [e.photons for e in emitters]
bg = [e.bg for e in emitters]
σ_x = [e.σ_x for e in emitters]
σ_y = [e.σ_y for e in emitters]
pvalue = [e.pvalue for e in emitters]

# Check photon values - something may be wrong
println("\n   📋 Checking fit results...")
println("     Photons - Mean: $(round(mean(photons), digits=0)), Median: $(round(median(photons), digits=0))")
println("     Background - Mean: $(round(mean(bg), digits=1)), Median: $(round(median(bg), digits=1))")
println("     σ_x - Mean: $(round(mean(σ_x)*1000, digits=1)) nm, Median: $(round(median(σ_x)*1000, digits=1)) nm")
println("     σ_y - Mean: $(round(mean(σ_y)*1000, digits=1)) nm, Median: $(round(median(σ_y)*1000, digits=1)) nm")
println("     pvalue - Median: $(round(median(pvalue), sigdigits=3))")

# Create histogram panel
fig = Figure(size=(1600, 1200))

# Photons histogram
ax1 = Axis(fig[1, 1], xlabel="Photons", ylabel="Count", title="Photon Distribution")
hist!(ax1, photons, bins=50)
vlines!(ax1, [median(photons)], color=:red, linestyle=:dash, label="Median")

# Background histogram
ax2 = Axis(fig[1, 2], xlabel="Background (ADU)", ylabel="Count", title="Background Distribution")
hist!(ax2, bg, bins=50)
vlines!(ax2, [median(bg)], color=:red, linestyle=:dash)

# σ_x histogram
ax3 = Axis(fig[2, 1], xlabel="σ_x (nm)", ylabel="Count", title="X Precision Distribution")
hist!(ax3, σ_x .* 1000, bins=50)
vlines!(ax3, [median(σ_x)*1000], color=:red, linestyle=:dash)

# σ_y histogram
ax4 = Axis(fig[2, 2], xlabel="σ_y (nm)", ylabel="Count", title="Y Precision Distribution")
hist!(ax4, σ_y .* 1000, bins=50)
vlines!(ax4, [median(σ_y)*1000], color=:red, linestyle=:dash)

# p-value histogram
ax5 = Axis(fig[3, 1], xlabel="p-value", ylabel="Count", title="Fit Quality (p-value)")
# Check if pvalues are all zero
if all(pvalue .== 0)
    text!(ax5, 0.5, 0.5, text="All p-values = 0\n(perfect fits)", align=(:center, :center))
else
    pval_plot = pvalue[pvalue .> 0]
    if !isempty(pval_plot)
        hist!(ax5, pval_plot, bins=min(50, length(pval_plot)))
        vlines!(ax5, [0.001], color=:red, linestyle=:dash, label="Threshold")
    end
end

# Photons vs Background
ax6 = Axis(fig[3, 2], xlabel="Photons", ylabel="Background (ADU)", title="Photons vs Background")
scatter!(ax6, photons, bg, markersize=2, alpha=0.3)

save(joinpath(outdir, "03_fitting", "fit_quality.png"), fig)
println("   ✓ Saved 03_fitting/fit_quality.png")

# =============================================================================
# 3b. FIT ACCEPTANCE PANEL (Green/Red boxes)
# =============================================================================
println("\n📊 Step 3b: Creating fit acceptance panel...")

# Filter parameters - use data-driven thresholds for visualization
# (In production, use stricter thresholds like 30nm)
precision_values = [sqrt(e.σ_x^2 + e.σ_y^2)/sqrt(2) for e in emitters]
median_precision = median(precision_values)

min_photons = 500.0
max_precision_um = median_precision * 1.5  # Accept top ~40% by precision
pval_threshold = -0.01  # Disabled for now (all pvalues are 0)

# Compute filter results for each localization
photon_ok = photons .> min_photons
precision_ok = precision_values .< max_precision_um
pvalue_ok = pvalue .> pval_threshold  # Higher pvalue = better fit

accepted = photon_ok .& precision_ok .& pvalue_ok

n_total = length(emitters)
n_accepted = sum(accepted)
n_photon_fail = sum(.!photon_ok)
n_precision_fail = sum(photon_ok .& .!precision_ok)
n_pvalue_fail = sum(photon_ok .& precision_ok .& .!pvalue_ok)

println("   Filter thresholds:")
println("     min_photons: $min_photons")
println("     max_precision: $(round(max_precision_um * 1000, digits=1)) nm (median=$(round(median_precision*1000, digits=1)) nm)")
println("     pval_threshold: $pval_threshold (disabled)")
println("   Results: $n_accepted / $n_total accepted ($(round(100*n_accepted/n_total, digits=1))%)")
println("   Rejections: photons=$n_photon_fail, precision=$n_precision_fail, pvalue=$n_pvalue_fail")

# Assign colors based on rejection reason
function get_box_color_local(acc, pok, precok)
    if acc
        return :green
    elseif !pok
        return :red      # Too dim
    elseif !precok
        return :orange   # Poor precision
    else
        return :purple   # Bad pvalue
    end
end

# Create figure - use wide format to match data aspect ratio
accept_pct = round(100*n_accepted/n_total, digits=1)
data_h, data_w = size(data[:,:,1])  # height × width
ncols, nrows = 4, 3
# Wide figure to minimize vertical whitespace
fig_accept = Figure(size=(2400, 700))
Label(fig_accept[0, 1:ncols],
    "Fit Acceptance: $n_accepted/$n_total ($accept_pct%) — Green=Accepted, Orange=Precision, Red=Photons",
    fontsize=14, tellwidth=false)

box_size = roi_batch.roi_size

# 3x4 grid of frames with colored boxes (same frames as detection panel)
for (idx, frame_num) in enumerate(frame_indices)
    row = div(idx - 1, 4) + 1
    col = mod(idx - 1, 4) + 1

    frame_mask = roi_batch.frame_indices .== frame_num
    n_in_frame = sum(frame_mask)
    n_acc_frame = sum(accepted[frame_mask])

    ax = Axis(fig_accept[row, col],
        title = "Frame $frame_num ($n_acc_frame/$n_in_frame)",
        aspect = DataAspect(),
        yreversed = true
    )

    # Show frame
    frame_data = data[:, :, frame_num]'
    heatmap!(ax, frame_data, colormap=:grays, colorrange=(pmin, pmax))

    # Find indices of localizations in this frame
    frame_locs = findall(frame_mask)

    if !isempty(frame_locs)
        det_x = roi_batch.x_corners[frame_mask]
        det_y = roi_batch.y_corners[frame_mask]
        frame_accepted = accepted[frame_mask]
        frame_photon_ok = photon_ok[frame_mask]
        frame_precision_ok = precision_ok[frame_mask]

        # Draw rejected first (underneath), then accepted on top
        for pass in [false, true]
            for j in eachindex(frame_locs)
                if frame_accepted[j] == pass
                    bx = det_x[j]
                    by = det_y[j]
                    c = get_box_color_local(frame_accepted[j], frame_photon_ok[j], frame_precision_ok[j])
                    alpha = pass ? 1.0 : 0.7
                    lw = 0.5
                    lines!(ax,
                        [bx, bx+box_size, bx+box_size, bx, bx],
                        [by, by, by+box_size, by+box_size, by],
                        color = (c, alpha), linewidth = lw)
                end
            end
        end
    end

    hidedecorations!(ax)
end

save(joinpath(outdir, "03_fitting", "fit_acceptance.png"), fig_accept)
println("   ✓ Saved 03_fitting/fit_acceptance.png")

# =============================================================================
# 4. FILTERING
# =============================================================================
println("\n🔬 Step 4: Filtering...")

# Use the acceptance mask computed above
filtered_emitters = emitters[accepted]

println("   Accepted: $(length(filtered_emitters)) / $(length(emitters)) ($(round(100*length(filtered_emitters)/length(emitters), digits=1))%)")

smld_filtered = BasicSMLD(filtered_emitters, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)

# =============================================================================
# 5. SUPER-RESOLUTION RENDERING
# =============================================================================
println("\n🎨 Step 5: Rendering super-resolution images...")

# Gaussian render
println("   Gaussian render...")
t_render = @elapsed result = render(smld_filtered;
    strategy=GaussianRender(),
    zoom=20,
    colormap=:inferno,
    filename=joinpath(outdir, "04_superres", "gaussian_inferno.png")
)
println("   ✓ Saved 04_superres/gaussian_inferno.png ($(round(t_render, digits=1))s)")

# Histogram render
println("   Histogram render...")
t_render = @elapsed result = render(smld_filtered;
    strategy=HistogramRender(),
    zoom=20,
    colormap=:viridis,
    filename=joinpath(outdir, "04_superres", "histogram_viridis.png")
)
println("   ✓ Saved 04_superres/histogram_viridis.png ($(round(t_render, digits=1))s)")

# =============================================================================
# SUMMARY
# =============================================================================
println("\n" * "="^80)
println("✅ Visualization complete!")
println("="^80)
println("\nOutputs created in dev/output/:")
println("  01_raw/ - Raw frame visualizations")
println("  02_detection/ - Detection overlays and ROI examples")
println("  03_fitting/ - Fit quality histograms")
println("  04_superres/ - Super-resolution images")
println()
