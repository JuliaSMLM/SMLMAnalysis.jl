"""
    basic_workflow.jl

Demonstrates the basic simulate → detect → fit workflow using SMLMAnalysis.

This example shows:
1. Setting up a simulation with realistic parameters
2. Running the complete workflow in one function call
3. Accessing and inspecting results
4. Viewing the recorded workflow steps
5. Rendering localizations to PNG images

To run in VSCode:
- Open this file and click "Run Code" or use Julia REPL
- Output images will be saved to examples/output/
"""

# Activate the examples environment (allows running from VSCode)
import Pkg
Pkg.activate(@__DIR__)

using SMLMAnalysis

# Setup camera (from SMLMRender octamer example)
println("Setting up camera...")
camera = IdealCamera(128, 128, 0.1)  # 128×128 pixels, 100nm pixel size
println("Camera: 128×128 pixels, 100nm pixel size\n")

# Setup simulation parameters (optimized for localization)
println("Setting up simulation...")
sim_params = StaticSMLMParams(
    density = 2.0,          # 2 patterns per μm² (16 total emitters/μm²)
    σ_psf = 0.13,           # 130nm PSF width
    nframes = 500,          # Full acquisition
    framerate = 20.0,       # 20 fps
    ndatasets = 2,          # 2 datasets
    ndims = 2               # 2D simulation
)

pattern = Nmer2D(n=8, d=0.15)  # Octamer, 150nm diameter

fluor = GenericFluor(
    photons = 20000.0,  # 20k photons/sec → ~1000 photons/frame @ 20fps
    k_off = 20.0,       # Switch off at 20 Hz
    k_on = 0.06         # Lower on-rate for sparse labeling
)

println("Simulation params:")
println("  Pattern: 8-mer (150nm diameter)")
println("  Density: $(sim_params.density) patterns/μm²")
println("  PSF σ: $(sim_params.σ_psf) μm")
println("  Frames: $(sim_params.nframes)\n")

# Run complete workflow
println("Running complete workflow...")
println("="^60)

detect_params = (
    boxsize = 11,
    overlap = 2.0,
    minval = 5.0,    # Lower threshold for sparse detections
    use_gpu = false  # Set to true if CUDA is available
)

fit_params = (
    iterations = 20,
    device = :cpu     # Or :gpu, :auto
)

result = simulate_detect_fit_workflow(
    sim_params,
    camera;
    pattern = pattern,
    molecule = fluor,
    detect_params = detect_params,
    fit_params = fit_params
)

println("="^60)
println("\nWorkflow complete!\n")

# Save intermediate images for inspection
output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)

println("Saving intermediate images...")

# Save example raw camera frames
using FileIO, ImageCore, Statistics

# Write statistics to file
stats_file = joinpath(output_dir, "workflow_statistics.txt")
open(stats_file, "w") do io
    println(io, "="^70)
    println(io, "SMLM WORKFLOW STATISTICS")
    println(io, "="^70)
    println(io)

    # Simulation parameters
    println(io, "SIMULATION PARAMETERS:")
    println(io, "  Pattern: 8-mer, diameter 150nm")
    println(io, "  Density: $(sim_params.density) patterns/μm²")
    println(io, "  Total emitters: $(sim_params.density * 8) emitters/μm²")
    println(io, "  PSF σ: $(sim_params.σ_psf) μm")
    println(io, "  Frames: $(sim_params.nframes)")
    println(io, "  k_on: $(fluor.γ), k_off: $(fluor.q)")
    println(io)

    # Camera info
    npix_x = length(camera.pixel_edges_x) - 1
    npix_y = length(camera.pixel_edges_y) - 1
    pix_size = camera.pixel_edges_x[2] - camera.pixel_edges_x[1]
    fov_x = npix_x * pix_size
    fov_y = npix_y * pix_size
    fov_area = fov_x * fov_y
    println(io, "CAMERA:")
    println(io, "  Size: $(npix_x)×$(npix_y) pixels")
    println(io, "  Pixel size: $(pix_size*1000) nm")
    println(io, "  FOV: $(round(fov_x, digits=2))×$(round(fov_y, digits=2)) μm")
    println(io, "  FOV area: $(round(fov_area, digits=2)) μm²")
    println(io)

    # Ground truth stats
    n_gt = length(result.smld_ground_truth.emitters)
    println(io, "GROUND TRUTH:")
    println(io, "  Total emitters: $n_gt")
    println(io, "  Emitters/frame: $(round(n_gt/sim_params.nframes, digits=1))")
    println(io, "  Emitters/μm²/frame: $(round(n_gt/sim_params.nframes/fov_area, digits=2))")
    if n_gt > 0
        gt_photons = [e.photons for e in result.smld_ground_truth.emitters]
        println(io, "  Mean photons: $(round(mean(gt_photons), digits=1))")
        println(io, "  Photon range: $(round(minimum(gt_photons), digits=1)) - $(round(maximum(gt_photons), digits=1))")
    end
    println(io)

    # Image statistics
    println(io, "RAW IMAGES:")
    println(io, "  Stack size: $(size(result.images))")
    img_stats = [extrema(result.images[:,:,i]) for i in 1:min(10, size(result.images, 3))]
    mean_min = mean([s[1] for s in img_stats])
    mean_max = mean([s[2] for s in img_stats])
    println(io, "  Mean intensity range (first 10 frames): $(round(mean_min, digits=1)) - $(round(mean_max, digits=1))")
    println(io)

    # Detection stats
    n_detected = length(result.detections)
    println(io, "DETECTION:")
    println(io, "  Threshold (minval): $(detect_params.minval)")
    println(io, "  Detections: $n_detected")
    println(io, "  Detection rate: $(round(n_detected/n_gt*100, digits=2))%")
    println(io, "  Detections/frame: $(round(n_detected/sim_params.nframes, digits=1))")
    println(io)

    # Fitting stats
    if n_detected > 0
        n_fitted = length(result.smld_fitted.emitters)
        println(io, "FITTING:")
        println(io, "  Fitted: $n_fitted")
        println(io, "  Fit success rate: $(round(n_fitted/n_detected*100, digits=1))%")
        fit_photons = [e.photons for e in result.smld_fitted.emitters]
        fit_bg = [e.bg for e in result.smld_fitted.emitters]
        println(io, "  Mean photons: $(round(mean(fit_photons), digits=1))")
        println(io, "  Mean background: $(round(mean(fit_bg), digits=1))")
        x_errors = [e.σ_x for e in result.smld_fitted.emitters]
        y_errors = [e.σ_y for e in result.smld_fitted.emitters]
        println(io, "  Mean σ_x: $(round(mean(x_errors)*1000, digits=1)) nm")
        println(io, "  Mean σ_y: $(round(mean(y_errors)*1000, digits=1)) nm")
    else
        println(io, "FITTING:")
        println(io, "  No fits (no detections)")
    end
    println(io)
    println(io, "="^70)
end
println("  Saved workflow_statistics.txt")
println()
for i in 1:min(3, size(result.images, 3))
    frame = result.images[:, :, i]
    # Normalize to 0-1 (handle case where min==max)
    fmin, fmax = extrema(frame)
    frame_norm = if fmax > fmin
        (frame .- fmin) ./ (fmax - fmin)
    else
        zeros(Float32, size(frame))
    end
    save(joinpath(output_dir, "raw_frame_$(lpad(i,3,'0')).png"), Gray.(frame_norm))
    println("  Saved raw_frame_$(lpad(i,3,'0')).png")
end

# Save example ROI boxes if any were detected
if length(result.detections) > 0
    n_rois_to_save = min(9, length(result.detections))
    roi_montage = zeros(Float32, result.detections.roi_size * 3, result.detections.roi_size * 3)
    for i in 1:n_rois_to_save
        row_idx = div(i-1, 3)
        col_idx = mod(i-1, 3)
        roi_data = result.detections.data[:, :, i]
        roi_norm = (roi_data .- minimum(roi_data)) ./ (maximum(roi_data) - minimum(roi_data) + 1e-10)

        r_start = row_idx * result.detections.roi_size + 1
        r_end = (row_idx + 1) * result.detections.roi_size
        c_start = col_idx * result.detections.roi_size + 1
        c_end = (col_idx + 1) * result.detections.roi_size

        roi_montage[r_start:r_end, c_start:c_end] = roi_norm
    end
    save(joinpath(output_dir, "roi_examples.png"), Gray.(roi_montage))
    println("  Saved roi_examples.png ($(n_rois_to_save) ROIs)")
end

println()

# Inspect results
println("="^60)
println("RESULTS SUMMARY")
println("="^60)

println("\nGround Truth:")
println("  Emitters: $(length(result.smld_ground_truth.emitters))")
println("  Frames: $(result.smld_ground_truth.n_frames)")

println("\nDetection:")
println("  Detections: $(length(result.detections))")
println("  Box size: $(result.detections.roi_size)×$(result.detections.roi_size)")

println("\nFitting:")
println("  Fitted: $(length(result.smld_fitted.emitters))")
println("  Mean photons: $(round(mean([e.photons for e in result.smld_fitted.emitters]), digits=1))")
println("  Mean background: $(round(mean([e.bg for e in result.smld_fitted.emitters]), digits=1))")

# Calculate precision
x_errors = [e.σ_x for e in result.smld_fitted.emitters]
y_errors = [e.σ_y for e in result.smld_fitted.emitters]
println("  Mean σ_x: $(round(mean(x_errors)*1000, digits=1)) nm")
println("  Mean σ_y: $(round(mean(y_errors)*1000, digits=1)) nm")

# Show workflow
println("\n")
println("="^60)
println("WORKFLOW PROVENANCE")
println("="^60)
println(result.workflow)

println("\n")
println("="^60)
println("RENDERING LOCALIZATIONS")
println("="^60)

# Create output directory if it doesn't exist
output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)

# Render ground truth with SMLMRender
println("\nRendering ground truth...")
gt_histogram = joinpath(output_dir, "ground_truth_histogram.png")
gt_gaussian = joinpath(output_dir, "ground_truth_gaussian.png")

render_result_gt = render(result.smld_ground_truth;
    strategy=GaussianRender(), zoom=20, colormap=:inferno, filename=gt_gaussian)
println("  Saved: $gt_gaussian")
println("  Image size: $(size(render_result_gt.image))")

# Render fitted localizations
println("\nRendering fitted localizations...")
fitted_gaussian = joinpath(output_dir, "fitted_gaussian.png")

render_result_fit = render(result.smld_fitted;
    strategy=GaussianRender(), zoom=20, colormap=:inferno, filename=fitted_gaussian)
println("  Saved: $fitted_gaussian")
println("  Image size: $(size(render_result_fit.image))")

println("\n")
println("="^60)
println("Example complete!")
println("="^60)
println("\nResults available in `result` variable:")
println("  result.smld_fitted          - Fitted localizations")
println("  result.smld_ground_truth    - Ground truth from simulation")
println("  result.images               - Generated camera images")
println("  result.detections           - Detection results (ROIBatch)")
println("  result.workflow             - Complete workflow record")
println("\nRendered images saved to: $output_dir/")
println("  - ground_truth_histogram.png  (ground truth, histogram render)")
println("  - ground_truth_gaussian.png   (ground truth, Gaussian render)")
println("  - fitted_histogram.png        (fitted locs, histogram render)")
println("  - fitted_gaussian.png         (fitted locs, Gaussian render)")
