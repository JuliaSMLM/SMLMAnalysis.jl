# Gatta Beads Analysis - Full Pipeline with DetectFitConfig
#
# Uses DetectFitConfig for combined detect+fit step.
# This is a CONTINUOUS dataset (single acquisition), so n_datasets=1.
#
# Pipeline: detectfit → filter → frameconnect → drift → isolated → render

using SMLMAnalysis

println("="^60)
println("GATTA BEADS ANALYSIS")
println("="^60)

# =============================================================================
# Data Path
# =============================================================================
h5file = joinpath(@__DIR__, "..", "data", "gatta_ruler", "2025-10-23",
    "20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5")

# Check file info
info = load_smart_h5_info(h5file)
println("Data: $(info.nframes) frames, $(info.width)×$(info.height), $(round(info.file_size_gb, digits=2)) GB")

# =============================================================================
# Camera Setup - ORCA-Fusion (C14440-20UP)
# =============================================================================
# 78nm pixels, 0.7 e-/ADU readnoise equivalent
# SCMOSCamera(nx, ny, ...) where nx=cols, ny=rows
# After smart_h5_to_array: data is (rows, cols, frames) = (height, width, frames)
camera = SCMOSCamera(info.width, info.height, 0.078f0, 0.7f0;
    offset = 100.0f0, gain = 0.24f0, qe = 1.0f0)
println("Camera: $(info.width)×$(info.height) (nx×ny), 78nm pixels")

# =============================================================================
# Analysis Setup
# =============================================================================
outdir = joinpath(@__DIR__, "output", "gatta_beads")

# Create Analysis without data (data comes from DetectFitConfig)
a = Analysis(camera; outdir, verbose=Verbosity.DETAILED, checkpoint=true)
println("Output: $outdir")
println()

# =============================================================================
# DetectFit - Single continuous acquisition
# =============================================================================
println("--- DETECTFIT ---")
run_step!(a, DetectFitConfig(
    path = h5file,
    n_datasets = 1,  # Continuous acquisition
    # Detection
    boxsize = 11,
    min_photons = 1000.0,
    psf_sigma = 0.135,
    backend = :auto,
    # Fitting - variable sigma
    psf_model = :variable,
    iterations = 20,
    # Filter preview (for fit_quality plot)
    filter_min_photons = 500.0,
    filter_max_precision = 0.007,  # 7nm
    filter_min_pvalue = 1e-6
))

# =============================================================================
# Filter
# =============================================================================
println("\n--- FILTER ---")
run_step!(a, FilterConfig(
    photons = (500.0, Inf),
    precision = (0.0, 0.007),  # max 7nm
    pvalue = (1e-6, 1.0),
    psf_sigma = :auto  # mode ± 10%
))

# =============================================================================
# Frame Connection + Uncertainty Calibration
# =============================================================================
println("\n--- FRAMECONNECT ---")
run_step!(a, FrameConnectConfig(
    maxframegap = 5,
    nsigmadev = 5.0,
    calibrate = true
))

# =============================================================================
# Drift Correction (with 4 chunks for 20k frames)
# =============================================================================
println("\n--- DRIFTCORRECT ---")
run_step!(a, DriftCorrectConfig(
    degree = 3,
    continuous = true,  # Continuous acquisition
    n_chunks = 4,       # Split into 4 chunks (5000 frames each)
    quality = :iterative
))

# =============================================================================
# Isolated Emitter Filter
# =============================================================================
println("\n--- ISOLATED ---")
run_step!(a, IsolatedConfig(
    n_sigma = 2.0,
    min_neighbors = :auto
))

# =============================================================================
# Render
# =============================================================================
println("\n--- RENDER ---")
run_step!(a, RenderConfig(
    renders = [
        RenderSpec(strategy=:gaussian, zoom=20, colormap=:inferno, clip_percentile=0.999),
        RenderSpec(strategy=:histogram, zoom=10, colormap=:turbo, color_by=:absolute_frame),
    ]
))

# =============================================================================
# Summary
# =============================================================================
println("\n" * "="^60)
println("ANALYSIS COMPLETE")
println("="^60)
println(a)

# Quick stats
if a.smld !== nothing && length(a.smld.emitters) > 0
    using Statistics
    emitters = a.smld.emitters
    photons = [e.photons for e in emitters]
    σ_x = [e.σ_x for e in emitters]
    σ_y = [e.σ_y for e in emitters]

    println("\nFinal localizations:")
    println("  Count: $(length(emitters))")
    println("  Photons: median=$(round(median(photons), digits=0))")
    println("  Precision: σ_x=$(round(median(σ_x)*1000, digits=1))nm, σ_y=$(round(median(σ_y)*1000, digits=1))nm")
end

println("\nOutput directory: $outdir")
