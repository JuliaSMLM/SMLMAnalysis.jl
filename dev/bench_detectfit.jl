# DetectFit Performance Benchmark
#
# Profiles the detectfit step with hexabody dataset to identify bottlenecks.
# Run on descent for 10gbps connection to NAS.
#
# Usage: julia --project=. dev/bench_detectfit.jl

using SMLMAnalysis
using Statistics
using Printf

println("="^60)
println("DETECTFIT PERFORMANCE BENCHMARK")
println("="^60)
println("Host: ", gethostname())
println()

# =============================================================================
# Data Path - Hexabody dataset
# =============================================================================
h5file = "/mnt/nas/cellpath/Genmab/Data/20250603_A431_SaturatingIgG10min+C1q/A431_IgG1-2F8-AF647_5ugml_10min+C1q/Cell_01/Label_01/Data_2025-6-4-17-24-11.h5"

# Check file exists
if !isfile(h5file)
    error("Data file not found: $h5file")
end

info = load_lidkelab_h5_info(h5file)
println("Data file: $(basename(h5file))")
println("  Total frames: $(info.n_frames)")
println("  Datasets (blocks): $(info.n_blocks)")
println("  Frames per dataset: $(info.frames_per_block[1])")
println("  Image size: $(info.width)×$(info.height)")
println("  File size: $(round(info.file_size_gb, digits=2)) GB")
println()

# =============================================================================
# Camera Setup
# =============================================================================
pixel_size = 0.0978f0  # microns
cal = load_lidkelab_h5_calibration_for_scmos(h5file)
camera = SCMOSCamera(info.width, info.height, pixel_size, cal.readnoise;
    offset = cal.offset, gain = cal.gain, qe = 0.82f0)
println("Camera: $(info.width)×$(info.height), $(pixel_size*1000)nm pixels")
println()

# =============================================================================
# Benchmark Parameters
# =============================================================================
# How many datasets to benchmark (use all for full benchmark, fewer for quick test)
n_datasets_bench = min(5, info.n_blocks)  # Benchmark first 5 datasets

# Detection parameters
boxsize = 9
min_photons = 500.0
psf_sigma = 0.130
backend = :auto

# Fit parameters
psf_model = :variable
iterations = 20

println("Benchmark config:")
println("  Datasets to process: $n_datasets_bench of $(info.n_blocks)")
println("  Detection: boxsize=$boxsize, min_photons=$min_photons, psf_sigma=$psf_sigma, backend=$backend")
println("  Fitting: model=$psf_model, iterations=$iterations")
println()

# =============================================================================
# Setup fitter and timing arrays
# =============================================================================
psf = if psf_model == :fixed
    GaussianXYNB(Float32(psf_sigma))
elseif psf_model == :variable
    GaussianXYNBS()
elseif psf_model == :anisotropic
    GaussianXYNBSXSY()
end
fitter = GaussMLEConfig(psf_model=psf, iterations=iterations)

# Timing storage
load_times = Float64[]
detect_times = Float64[]
fit_times = Float64[]
rois_per_dataset = Int[]
fits_per_dataset = Int[]

# =============================================================================
# Warmup (compile everything with first dataset)
# =============================================================================
println("--- WARMUP (compiling) ---")
frames_per_ds = info.frames_per_block[1]

t_warmup = @elapsed begin
    warmup_images = load_lidkelab_h5_block(h5file, 1)
    (warmup_roi, _) = getboxes(warmup_images, camera;
        boxsize=boxsize, overlap=2.0, min_photons=min_photons, psf_sigma=psf_sigma, backend=backend)
    (warmup_smld, _) = fit(warmup_roi, fitter)
end
println("Warmup complete: $(round(t_warmup, digits=2))s")
println()

# Clear warmup data
warmup_images = nothing
warmup_roi = nothing
warmup_smld = nothing
GC.gc()

# =============================================================================
# Benchmark Loop - Using BLOCK-BASED loading (efficient)
# =============================================================================
println("--- BENCHMARK ($n_datasets_bench datasets, BLOCK-BASED loading) ---")
println()

for ds in 1:n_datasets_bench
    print("Dataset $ds/$n_datasets_bench (block $ds)... ")

    # 1. Load - USE BLOCK-BASED LOADING (efficient)
    t_load = @elapsed begin
        images = load_lidkelab_h5_block(h5file, ds)
    end
    push!(load_times, t_load)

    # 2. Detect
    t_detect = @elapsed begin
        (roi_batch, _) = getboxes(images, camera;
            boxsize=boxsize, overlap=2.0, min_photons=min_photons, psf_sigma=psf_sigma, backend=backend)
    end
    push!(detect_times, t_detect)
    n_rois = length(roi_batch)
    push!(rois_per_dataset, n_rois)

    # 3. Fit
    t_fit = @elapsed begin
        (smld_ds, _) = fit(roi_batch, fitter)
    end
    push!(fit_times, t_fit)
    n_fits = length(smld_ds.emitters)
    push!(fits_per_dataset, n_fits)

    total_ds = t_load + t_detect + t_fit
    println(@sprintf("%.2fs (load=%.2f, detect=%.2f, fit=%.2f) → %d ROIs, %d fits",
                     total_ds, t_load, t_detect, t_fit, n_rois, n_fits))

    # Free memory
    images = nothing
    roi_batch = nothing
    smld_ds = nothing
    GC.gc()
end

# =============================================================================
# Results Summary
# =============================================================================
println()
println("="^60)
println("BENCHMARK RESULTS")
println("="^60)

total_load = sum(load_times)
total_detect = sum(detect_times)
total_fit = sum(fit_times)
total_time = total_load + total_detect + total_fit

total_rois = sum(rois_per_dataset)
total_fits = sum(fits_per_dataset)
total_frames = n_datasets_bench * frames_per_ds

println()
println("TOTAL TIME: $(round(total_time, digits=2))s for $n_datasets_bench datasets")
println("  - Load:   $(round(total_load, digits=2))s ($(round(100*total_load/total_time, digits=1))%)")
println("  - Detect: $(round(total_detect, digits=2))s ($(round(100*total_detect/total_time, digits=1))%)")
println("  - Fit:    $(round(total_fit, digits=2))s ($(round(100*total_fit/total_time, digits=1))%)")
println()

println("PER-DATASET AVERAGES:")
println("  - Load:   $(round(mean(load_times), digits=3))s ± $(round(std(load_times), digits=3))s")
println("  - Detect: $(round(mean(detect_times), digits=3))s ± $(round(std(detect_times), digits=3))s")
println("  - Fit:    $(round(mean(fit_times), digits=3))s ± $(round(std(fit_times), digits=3))s")
println()

println("THROUGHPUT:")
println("  - Frames/s: $(round(total_frames/total_time, digits=1))")
println("  - ROIs/s: $(round(total_rois/total_time, digits=1))")
println("  - Fits/s: $(round(total_fits/total_time, digits=1))")
println("  - Fit throughput: $(round(total_fits/total_fit/1000, digits=1))k fits/s (fitting only)")
println()

println("DATA SUMMARY:")
println("  - Total frames processed: $total_frames")
println("  - Total ROIs detected: $total_rois")
println("  - Total fits: $total_fits")
println("  - ROIs per frame: $(round(total_rois/total_frames, digits=1))")
println("  - Fit rate: $(round(100*total_fits/total_rois, digits=1))%")
println()

# Extrapolate to full dataset
if n_datasets_bench < info.n_blocks
    println("EXTRAPOLATED (full $(info.n_blocks) datasets):")
    scale = info.n_blocks / n_datasets_bench
    println("  - Total time: $(round(total_time * scale / 60, digits=1)) minutes")
    println("  - Total fits: ~$(round(Int, total_fits * scale / 1000))k localizations")
end
println()

# =============================================================================
# Detailed Breakdown
# =============================================================================
println("="^60)
println("DETAILED BREAKDOWN (per dataset)")
println("="^60)
println()
println("Dataset │ Load (s) │ Detect (s) │ Fit (s) │ ROIs │ Fits │ Fits/s")
println("────────┼──────────┼────────────┼─────────┼──────┼──────┼───────")
for ds in 1:n_datasets_bench
    fits_per_sec = fits_per_dataset[ds] / fit_times[ds]
    println(@sprintf("   %2d   │  %6.3f  │   %6.3f   │ %6.3f  │ %5d│ %5d│ %6.0f",
                     ds, load_times[ds], detect_times[ds], fit_times[ds],
                     rois_per_dataset[ds], fits_per_dataset[ds], fits_per_sec))
end
println()

# =============================================================================
# Bottleneck Analysis
# =============================================================================
println("="^60)
println("BOTTLENECK ANALYSIS")
println("="^60)
println()

bottleneck = if total_load > total_detect && total_load > total_fit
    "DATA LOADING"
elseif total_detect > total_fit
    "DETECTION"
else
    "FITTING"
end

println("Primary bottleneck: $bottleneck")
println()

if bottleneck == "DATA LOADING"
    load_rate = (info.width * info.height * frames_per_ds * sizeof(UInt16)) / mean(load_times) / 1e9
    println("Load analysis:")
    println("  - Data rate: $(round(load_rate, digits=2)) GB/s")
    println("  - Expected max for 10gbps: ~1.2 GB/s")
    if load_rate < 0.8
        println("  → Network may be limiting factor")
    else
        println("  → Loading is near network maximum")
    end
elseif bottleneck == "DETECTION"
    println("Detection analysis:")
    println("  - Frames per dataset: $frames_per_ds")
    println("  - ROIs per frame: $(round(mean(rois_per_dataset)/frames_per_ds, digits=1))")
    println("  - Detection time per frame: $(round(mean(detect_times)*1000/frames_per_ds, digits=2))ms")
    if backend != :cpu
        println("  → GPU detection; may benefit from batching or kernel optimization")
    else
        println("  → CPU detection; consider enabling GPU")
    end
else
    println("Fitting analysis:")
    println("  - Fit time per ROI: $(round(mean(fit_times)*1000000/mean(rois_per_dataset), digits=2))μs")
    println("  - Iterations: $iterations")
    if backend != :cpu
        println("  → GPU fitting; check batch size and occupancy")
    else
        println("  → CPU fitting; enable GPU for speedup")
    end
end
println()
