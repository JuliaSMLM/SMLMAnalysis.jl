# Hexabody dSTORM Analysis - Multi-dataset Pipeline
#
# Lidke Lab H5 format with registered acquisitions (TYPE 2).
# Each H5 file contains one cell with multiple datasets (blocks) separated by
# stage registration.
#
# Data: A431 cells + IgG (2F8 wild-type, E345R, E430G, RGY mutants) + C1q
# 20 datasets × 5000 frames = 100k frames per cell
#
# Pipeline: detectfit → filter → frameconnect → drift → isolated → render

using SMLMAnalysis

println("="^60)
println("HEXABODY dSTORM ANALYSIS")
println("="^60)

# =============================================================================
# Data Path - Single cell file (wild-type 2F8)
# =============================================================================
h5file = "/mnt/nas/cellpath/Genmab/Data/20250603_A431_SaturatingIgG10min+C1q/A431_IgG1-2F8-AF647_5ugml_10min+C1q/Cell_01/Label_01/Data_2025-6-4-17-24-11.h5"

# Check file info
info = load_lidkelab_h5_info(h5file)
println("Data: $(info.n_frames) frames in $(info.n_blocks) datasets")
println("Frames per dataset: $(info.frames_per_block[1])")
println("Image size: $(info.width)×$(info.height)")
println("File size: $(round(info.file_size_gb, digits=2)) GB")

# =============================================================================
# Camera Setup - Hamamatsu sCMOS with per-pixel calibration
# =============================================================================
# 97.8nm pixels (from cellpath config)
pixel_size = 0.0978f0  # microns

# Load per-pixel calibration from H5 file
cal = load_lidkelab_h5_calibration_for_scmos(h5file)
using Statistics
println("Calibration: gain=$(round(median(cal.gain), digits=3)), offset=$(round(median(cal.offset), digits=1))")

# Create camera with per-pixel calibration
# SCMOSCamera(nx, ny, pixel_size, readnoise; offset, gain, qe)
camera = SCMOSCamera(info.width, info.height, pixel_size, cal.readnoise;
    offset = cal.offset, gain = cal.gain, qe = 0.82f0)
println("Camera: $(info.width)×$(info.height), $(pixel_size*1000)nm pixels (per-pixel calibration)")

# =============================================================================
# Analysis Setup
# =============================================================================
outdir = joinpath(@__DIR__, "output", "hexabody_dstorm")

# Create Analysis without data (data comes from DetectFitConfig)
a = Analysis(camera; outdir, verbose=Verbosity.DETAILED, checkpoint=true)  # Disk checkpoints after detectfit/frameconnect/drift only
println("Output: $outdir")
println()

# =============================================================================
# DetectFit - Multi-dataset registered acquisition
# =============================================================================
println("--- DETECTFIT ---")
run_step!(a, DetectFitConfig(
    path = h5file,
    h5_format = :mic,
    n_datasets = info.n_blocks,  # 20 datasets
    # Detection
    boxsize = 9,
    min_photons = 500.0,
    psf_sigma = 0.130,  # ~1.3 pixels at 97.8nm
    backend = :auto,
    # Fitting - variable sigma for dSTORM
    psf_model = :variable,
    iterations = 20,
    # Filter preview (for fit_quality plot)
    filter_min_photons = 500.0,
    filter_max_precision = 0.007,  # 7nm
    filter_min_pvalue = 1e-6
))

# =============================================================================
# Filter - Strict filtering with per-pixel calibration
# =============================================================================
println("\n--- FILTER ---")
run_step!(a, FilterConfig(
    photons = (500.0, Inf),    # 500 photon minimum
    precision = (0.0, 0.007),  # max 7nm precision (tight filter)
    pvalue = (1e-6, 1.0),      # p-value filter (per-pixel cal should work)
    psf_sigma = nothing        # No PSF filter
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
# Drift Correction - Registered mode (TYPE 2)
# =============================================================================
println("\n--- DRIFTCORRECT ---")
run_step!(a, DriftCorrectConfig(
    degree = 2,
    continuous = false,  # TYPE 2: registered acquisitions
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
run_step!(a, RenderConfig(zoom=20, colormap=:inferno))

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
