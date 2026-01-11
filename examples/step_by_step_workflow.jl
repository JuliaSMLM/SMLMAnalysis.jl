"""
    step_by_step_workflow.jl

Step-by-step composition of JuliaSMLM packages for troubleshooting.
Shows explicit calls to each package with results printed after each step.
"""

# Activate environment
import Pkg
Pkg.activate(@__DIR__)

using SMLMAnalysis
using MicroscopePSFs
using Statistics
using SMLMData: @filter

# Setup output
output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)
log_file = joinpath(output_dir, "step_by_step_log.txt")

# Open log file
io = open(log_file, "w")

function log_and_print(msg)
    println(msg)
    println(io, msg)
    flush(io)
end

log_and_print("="^70)
log_and_print("STEP-BY-STEP SMLM WORKFLOW")
log_and_print("="^70)
log_and_print("")

# ============================================================================
# STEP 1: Setup Camera
# ============================================================================
log_and_print("STEP 1: Setup Camera")
log_and_print("-"^70)

camera = IdealCamera(64, 128, 0.1)  # Rectangular: 64 wide × 128 tall
log_and_print(repr(camera))
log_and_print("  FOV: 6.4 μm × 12.8 μm (rectangular)")
log_and_print("")

# ============================================================================
# STEP 2: Setup Simulation Parameters
# ============================================================================
log_and_print("STEP 2: Setup Simulation")
log_and_print("-"^70)

sim_params = StaticSMLMParams(
    density = 2.0,
    σ_psf = 0.13,
    nframes = 500,      # More frames for better sampling
    ndatasets = 2
)
log_and_print("Simulation params:")
log_and_print("  density = $(sim_params.density) patterns/μm²")
log_and_print("  σ_psf = $(sim_params.σ_psf) μm")
log_and_print("  nframes = $(sim_params.nframes)")
log_and_print("  ndatasets = $(sim_params.ndatasets)")
log_and_print("")

pattern = Nmer2D(n=8, d=0.15)
log_and_print("Pattern: 8-mer, 150nm diameter")
log_and_print("")

fluor = GenericFluor(
    photons = 200000.0,  # 200k photons/sec → ~10k photons/frame @ 20fps (high SNR)
    k_off = 20.0,        # Switch off at 20 Hz (50ms on-time)
    k_on = 0.06          # Lower on-rate for sparser labeling
)
log_and_print("Fluorophore: photons=200000 Hz, k_off=20 Hz, k_on=0.06 Hz")
log_and_print("  Duty cycle: ~0.3%, avg on-time: 50ms, avg off-time: 16.7s")
log_and_print("")

# ============================================================================
# STEP 3: Simulate
# ============================================================================
log_and_print("STEP 3: Simulate SMLM Data")
log_and_print("-"^70)

t_sim = @elapsed begin
    pattern_result, smld_true, smld_noisy = simulate(sim_params;
        pattern=pattern,
        molecule=fluor,
        camera=camera
    )
end

log_and_print("Simulation time: $(round(t_sim, digits=2))s")
log_and_print("Ground truth SMLD:")
log_and_print("  Emitters: $(length(smld_noisy.emitters))")
log_and_print("  Frames: $(smld_noisy.n_frames)")
log_and_print("  Datasets: $(smld_noisy.n_datasets)")
if length(smld_noisy.emitters) > 0
    log_and_print("  Mean photons: $(round(mean([e.photons for e in smld_noisy.emitters]), digits=1))")
end
log_and_print("")

# ============================================================================
# STEP 4: Generate Camera Images
# ============================================================================
log_and_print("STEP 4: Generate Camera Images")
log_and_print("-"^70)

psf = MicroscopePSFs.GaussianPSF(sim_params.σ_psf)

t_img = @elapsed begin
    images = gen_images(smld_noisy, psf;
        bg=10.0,
        poisson_noise=true,
        camera_noise=false
    )
end

log_and_print("Image generation time: $(round(t_img, digits=2))s")
log_and_print("Image stack: $(size(images))")
log_and_print("Intensity range (frame 1): $(round(minimum(images[:,:,1]), digits=1)) - $(round(maximum(images[:,:,1]), digits=1))")
log_and_print("")

# Save example frames
using FileIO, ImageCore
for i in 1:min(3, size(images, 3))
    frame = images[:, :, i]
    fmin, fmax = extrema(frame)
    frame_norm = fmax > fmin ? (frame .- fmin) ./ (fmax - fmin) : zeros(Float32, size(frame))
    save(joinpath(output_dir, "raw_frame_$(lpad(i,3,'0')).png"), Gray.(frame_norm))
end
log_and_print("Saved raw_frame_001-003.png")
log_and_print("")

# ============================================================================
# STEP 5: Detect Particles
# ============================================================================
log_and_print("STEP 5: Detect Particles")
log_and_print("-"^70)

t_detect = @elapsed begin
    roi_batch = getboxes(images, camera;
        boxsize=7,
        overlap=2.0,
        sigma_small=1.0,
        sigma_large=2.0,
        minval=5.0,
        use_gpu=false
    )
end

log_and_print("Detection time: $(round(t_detect, digits=2))s")
log_and_print("Detections: $(length(roi_batch))")
log_and_print("ROI size: $(roi_batch.roi_size)×$(roi_batch.roi_size)")
log_and_print("Detection rate: $(round(length(roi_batch)/length(smld_noisy.emitters)*100, digits=1))%")

# Check ROI corners (in pixel coordinates)
if length(roi_batch) > 0
    x_corners = roi_batch.x_corners
    y_corners = roi_batch.y_corners
    log_and_print("")
    log_and_print("ROI Corners (pixel coordinates):")
    log_and_print("  x_corners: min=$(minimum(x_corners)), max=$(maximum(x_corners))")
    log_and_print("  y_corners: min=$(minimum(y_corners)), max=$(maximum(y_corners))")
    log_and_print("  Frame indices: min=$(minimum(roi_batch.frame_indices)), max=$(maximum(roi_batch.frame_indices))")
end
log_and_print("")

# Save ROI examples
if length(roi_batch) > 0
    n_show = min(9, length(roi_batch))
    roi_montage = zeros(Float32, roi_batch.roi_size * 3, roi_batch.roi_size * 3)
    for i in 1:n_show
        row_idx = div(i-1, 3)
        col_idx = mod(i-1, 3)
        roi_data = roi_batch.data[:, :, i]
        roi_norm = (roi_data .- minimum(roi_data)) ./ (maximum(roi_data) - minimum(roi_data) + 1e-10)

        r_start = row_idx * roi_batch.roi_size + 1
        r_end = (row_idx + 1) * roi_batch.roi_size
        c_start = col_idx * roi_batch.roi_size + 1
        c_end = (col_idx + 1) * roi_batch.roi_size

        roi_montage[r_start:r_end, c_start:c_end] = roi_norm
    end
    save(joinpath(output_dir, "roi_examples.png"), Gray.(roi_montage))
    log_and_print("Saved roi_examples.png ($n_show ROIs)")
end
log_and_print("")

# ============================================================================
# STEP 6: Fit Gaussian PSFs
# ============================================================================
if length(roi_batch) == 0
    log_and_print("STEP 6: Fit Gaussian PSFs")
    log_and_print("-"^70)
    log_and_print("SKIPPED - No detections to fit")
    log_and_print("")
    close(io)
    error("No detections - cannot continue. Check detection threshold or simulation params.")
end

log_and_print("STEP 6: Fit Gaussian PSFs")
log_and_print("-"^70)

fitter = GaussMLEFitter(
    psf_model = GaussianXYNB(Float32(sim_params.σ_psf)),
    iterations = 20,
    device = :cpu
)

t_fit = @elapsed begin
    smld_fitted = fit(fitter, roi_batch)
end

log_and_print("Fitting time: $(round(t_fit, digits=2))s")
log_and_print("Fitted SMLD:")
log_and_print("  Emitters: $(length(smld_fitted.emitters))")
log_and_print("  Frames: $(smld_fitted.n_frames)")

# Direct inspection of first emitter
if length(smld_fitted.emitters) > 0
    log_and_print("")
    log_and_print("First fitted emitter (direct inspection):")
    e1 = smld_fitted.emitters[1]
    log_and_print("  Position: ($(round(e1.x, digits=3)), $(round(e1.y, digits=3))) μm")
    log_and_print("  Photons: $(round(e1.photons, digits=1))")
    log_and_print("  Background: $(round(e1.bg, digits=1))")
    log_and_print("  σ_x: $(round(e1.σ_x*1000, digits=1)) nm")
    log_and_print("  σ_y: $(round(e1.σ_y*1000, digits=1)) nm")
    log_and_print("  Frame: $(e1.frame)")
end

if length(smld_fitted.emitters) > 0
    using Statistics
    log_and_print("")
    log_and_print("Statistics across all emitters:")
    log_and_print("  Mean photons: $(round(mean([e.photons for e in smld_fitted.emitters]), digits=1))")
    log_and_print("  Mean background: $(round(mean([e.bg for e in smld_fitted.emitters]), digits=1))")
    log_and_print("  Mean σ_x: $(round(mean([e.σ_x for e in smld_fitted.emitters])*1000, digits=1)) nm")
    log_and_print("  Mean σ_y: $(round(mean([e.σ_y for e in smld_fitted.emitters])*1000, digits=1)) nm")

    # Check fitted emitter positions
    xs = [e.x for e in smld_fitted.emitters]
    ys = [e.y for e in smld_fitted.emitters]
    frames = [e.frame for e in smld_fitted.emitters]

    log_and_print("")
    log_and_print("FITTED POSITIONS:")
    log_and_print("  x: min=$(round(minimum(xs), digits=2)) μm, max=$(round(maximum(xs), digits=2)) μm")
    log_and_print("  y: min=$(round(minimum(ys), digits=2)) μm, max=$(round(maximum(ys), digits=2)) μm")
    log_and_print("  x span: $(round(maximum(xs) - minimum(xs), digits=2)) μm")
    log_and_print("  y span: $(round(maximum(ys) - minimum(ys), digits=2)) μm")
    log_and_print("  Frames: min=$(minimum(frames)), max=$(maximum(frames))")

    # Check if coordinates are within camera FOV
    camera_x_min = minimum(camera.pixel_edges_x)
    camera_x_max = maximum(camera.pixel_edges_x)
    camera_y_min = minimum(camera.pixel_edges_y)
    camera_y_max = maximum(camera.pixel_edges_y)

    log_and_print("")
    log_and_print("COORDINATE CHECK vs CAMERA FOV:")
    log_and_print("  Camera: x=[$(camera_x_min), $(camera_x_max)] μm, y=[$(camera_y_min), $(camera_y_max)] μm")

    n_outside_x = count(x -> x < camera_x_min || x > camera_x_max, xs)
    n_outside_y = count(y -> y < camera_y_min || y > camera_y_max, ys)

    if n_outside_x > 0 || n_outside_y > 0
        log_and_print("  ⚠️  WARNING: $(n_outside_x) emitters outside camera x-range ($(round(n_outside_x/length(xs)*100, digits=1))%)")
        log_and_print("  ⚠️  WARNING: $(n_outside_y) emitters outside camera y-range ($(round(n_outside_y/length(ys)*100, digits=1))%)")
    else
        log_and_print("  ✓ All emitters within camera FOV")
    end

    # Check for edge bias - compare positions near edges
    log_and_print("")
    log_and_print("EDGE BIAS CHECK:")
    # Find emitters from edge ROIs (corners near boundary)
    edge_threshold = 10  # pixels from edge
    cam_width = length(camera.pixel_edges_x) - 1
    cam_height = length(camera.pixel_edges_y) - 1
    edge_indices = findall(i ->
        roi_batch.x_corners[i] <= edge_threshold ||
        roi_batch.x_corners[i] >= (cam_width - edge_threshold) ||
        roi_batch.y_corners[i] <= edge_threshold ||
        roi_batch.y_corners[i] >= (cam_height - edge_threshold),
        1:length(roi_batch)
    )

    if length(edge_indices) > 0
        log_and_print("  Edge ROIs (within 10 pixels of boundary): $(length(edge_indices))")
        # Show first 3 edge ROIs
        for (idx, i) in enumerate(edge_indices[1:min(3, length(edge_indices))])
            corner_x = roi_batch.x_corners[i]
            corner_y = roi_batch.y_corners[i]
            fitted_x = smld_fitted.emitters[i].x
            fitted_y = smld_fitted.emitters[i].y
            log_and_print("    ROI[$i]: corner=($corner_x, $corner_y) px → fitted=($(round(fitted_x, digits=2)), $(round(fitted_y, digits=2))) μm")
        end
    end
end
log_and_print("")

# ============================================================================
# STEP 7: Filter Localizations
# ============================================================================
log_and_print("STEP 7: Filter Localizations (p > 0.01 && σ < 20nm)")
log_and_print("-"^70)

# Filter using p-value and precision thresholds
smld_filtered = @filter(smld_fitted, pvalue > 0.01 && σ_x < 0.02 && σ_y < 0.02)

n_before = length(smld_fitted.emitters)
n_after = length(smld_filtered.emitters)
n_removed = n_before - n_after

log_and_print("Filtering results (p > 0.01 && σ_x,σ_y < 20nm):")
log_and_print("  Before: $n_before emitters")
log_and_print("  After: $n_after emitters")
log_and_print("  Removed: $n_removed ($(round(n_removed/n_before*100, digits=1))%)")

if n_after > 0
    using Statistics
    pvalues = [e.pvalue for e in smld_filtered.emitters]
    log_and_print("  Mean p-value: $(round(mean(pvalues), digits=6))")
    log_and_print("  Max p-value: $(round(maximum(pvalues), digits=6))")
    log_and_print("  Mean σ_x (filtered): $(round(mean([e.σ_x for e in smld_filtered.emitters])*1000, digits=1)) nm")
    log_and_print("  Mean σ_y (filtered): $(round(mean([e.σ_y for e in smld_filtered.emitters])*1000, digits=1)) nm")

    # Create p-value histogram
    using CairoMakie
    all_pvalues = [e.pvalue for e in smld_fitted.emitters]
    log_pvalues = log10.(all_pvalues .+ 1e-300)  # Add epsilon to avoid log(0)

    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1],
        xlabel = "log₁₀(p-value)",
        ylabel = "Count",
        title = "P-value Distribution (reject p < 0.01)"
    )

    hist!(ax, log_pvalues, bins=50, color=(:blue, 0.5))
    vlines!(ax, [log10(0.01)], color=:red, linewidth=2, linestyle=:dash,
            label="Threshold (p=0.01)")
    axislegend(ax, position=:lt)

    save(joinpath(output_dir, "pvalue_histogram.png"), fig)
    log_and_print("  Saved pvalue_histogram.png")
end
log_and_print("")

# ============================================================================
# STEP 8: Compare GT vs Filtered (Jaccard Index)
# ============================================================================
log_and_print("STEP 8: Compare GT vs Filtered (Jaccard Index)")
log_and_print("-"^70)

# Calculate Jaccard index with 100nm cutoff
if length(smld_filtered.emitters) > 0 && length(smld_noisy.emitters) > 0
    cutoff = 0.1  # 100nm in microns

    # Extract coordinates
    gt_coords = [(e.x, e.y) for e in smld_noisy.emitters]
    fit_coords = [(e.x, e.y) for e in smld_filtered.emitters]

    # Count matches within cutoff
    local n_matched = 0
    for (fx, fy) in fit_coords
        for (gx, gy) in gt_coords
            dist = sqrt((fx - gx)^2 + (fy - gy)^2)
            if dist < cutoff
                n_matched += 1
                break
            end
        end
    end

    # Jaccard = matched / (GT + Fitted - matched)
    jaccard = n_matched / (length(gt_coords) + length(fit_coords) - n_matched)

    log_and_print("Jaccard Index (100nm cutoff):")
    log_and_print("  Matched: $n_matched")
    log_and_print("  GT total: $(length(gt_coords))")
    log_and_print("  Fitted total: $(length(fit_coords))")
    log_and_print("  Jaccard: $(round(jaccard, digits=3))")
    log_and_print("")

    # Sample coordinate comparison
    log_and_print("Sample coordinate comparison (first 5):")
    for i in 1:min(5, length(smld_filtered.emitters))
        fx, fy = smld_filtered.emitters[i].x, smld_filtered.emitters[i].y
        log_and_print("  Filtered[$i]: x=$(round(fx, digits=3)) μm, y=$(round(fy, digits=3)) μm")
    end
    log_and_print("")
    for i in 1:min(5, length(smld_noisy.emitters))
        gx, gy = smld_noisy.emitters[i].x, smld_noisy.emitters[i].y
        log_and_print("  GT[$i]: x=$(round(gx, digits=3)) μm, y=$(round(gy, digits=3)) μm")
    end
end
log_and_print("")

# ============================================================================
# STEP 9: Render Results
# ============================================================================
log_and_print("STEP 9: Render Filtered Results")
log_and_print("-"^70)

# Render ground truth
log_and_print("Rendering ground truth (Gaussian, inferno)...")
t_render_gt = @elapsed begin
    result_gt = render(smld_noisy;
        strategy=GaussianRender(),
        zoom=20,
        colormap=:inferno,
        filename=joinpath(output_dir, "ground_truth_gaussian.png")
    )
end
log_and_print("  Saved ground_truth_gaussian.png ($(round(t_render_gt*1000, digits=1))ms)")
log_and_print("  Image size: $(size(result_gt.image))")

# Render filtered
log_and_print("Rendering filtered localizations (Gaussian, inferno)...")
t_render_fit = @elapsed begin
    result_fit = render(smld_filtered;
        strategy=GaussianRender(),
        zoom=20,
        colormap=:inferno,
        filename=joinpath(output_dir, "fitted_gaussian.png")
    )
end
log_and_print("  Saved fitted_gaussian.png ($(round(t_render_fit*1000, digits=1))ms)")
log_and_print("  Image size: $(size(result_fit.image))")
log_and_print("")

# ============================================================================
# Summary
# ============================================================================
log_and_print("="^70)
log_and_print("WORKFLOW COMPLETE")
log_and_print("="^70)
log_and_print("Total time: $(round(t_sim + t_img + t_detect + t_fit, digits=1))s")
log_and_print("")
log_and_print("Output files in: $output_dir/")
log_and_print("  - step_by_step_log.txt")
log_and_print("  - raw_frame_*.png")
log_and_print("  - roi_examples.png")
log_and_print("  - ground_truth_gaussian.png")
log_and_print("  - fitted_gaussian.png")

close(io)
println("\n✓ Log saved to: $log_file")
