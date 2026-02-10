# Hexabody dSTORM Analysis - Multi-dataset Pipeline
#
# Lidke Lab H5 format with registered acquisitions (TYPE 2).
# Each H5 file contains one cell with multiple datasets (blocks) separated by
# stage registration.
#
# Data: A431 cells + IgG (2F8 wild-type, E345R, E430G, RGY mutants) + C1q
# 20 datasets x 5000 frames = 100k frames per cell
#
# Pipeline: detectfit -> filter -> frameconnect -> drift -> densityfilter -> render

using SMLMAnalysis

println("="^60)
println("HEXABODY dSTORM ANALYSIS")
println("="^60)

# =============================================================================
# Data Path - Single cell file (wild-type 2F8)
# =============================================================================
h5file = "/mnt/nas/cellpath/Genmab/Data/20250603_A431_SaturatingIgG10min+C1q/A431_IgG1-2F8-AF647_5ugml_10min+C1q/Cell_01/Label_01/Data_2025-6-4-17-24-11.h5"

# Check file info
info = load_mic_h5_info(h5file)
println("Data: $(info.n_frames) frames in $(info.n_blocks) datasets")
println("Frames per dataset: $(info.frames_per_block[1])")
println("Image size: $(info.width)x$(info.height)")
println("File size: $(round(info.file_size_gb, digits=2)) GB")

# =============================================================================
# Camera Setup - Hamamatsu sCMOS with per-pixel calibration
# =============================================================================
# 97.8nm pixels (from cellpath config)
pixel_size = 0.0978f0  # microns

# Load per-pixel calibration from H5 file
cal = load_mic_h5_calibration_for_scmos(h5file)
using Statistics
println("Calibration: gain=$(round(median(cal.gain), digits=3)), offset=$(round(median(cal.offset), digits=1))")

# Create camera with per-pixel calibration
# SCMOSCamera(nx, ny, pixel_size, readnoise; offset, gain, qe)
camera = SCMOSCamera(info.width, info.height, pixel_size, cal.readnoise;
    offset = cal.offset, gain = cal.gain, qe = 0.82f0)
println("Camera: $(info.width)x$(info.height), $(pixel_size*1000)nm pixels (per-pixel calibration)")

# =============================================================================
# Analysis Setup - File-based (MIC format auto-detects blocks as datasets)
# =============================================================================
outdir = joinpath(@__DIR__, "output", "hexabody_dstorm")

config = AnalysisConfig(
    camera = camera,
    steps = [
        DetectFitConfig(
            path = h5file,
            h5_format = :mic,
            boxsize = 9,
            min_photons = 500.0,
            psf_sigma = 0.130,
            backend = :auto,
            psf_model = :variable,
            iterations = 20,
            filter_min_photons = 500.0,
            filter_max_precision = 0.007,
            filter_min_pvalue = 1e-6
        ),
        FilterConfig(
            photons = (500.0, Inf),
            precision = (0.0, 0.007),
            pvalue = (1e-6, 1.0)
        ),
        FrameConnectConfig(
            max_frame_gap = 5,
            max_sigma_dist = 5.0,
            calibrate = true
        ),
        DriftCorrectConfig(
            degree = 2,
            continuous = false,
            quality = :iterative
        ),
        DensityFilterConfig(
            n_sigma = 2.0,
            min_neighbors = :auto
        ),
        RenderConfig(zoom=20, colormap=:inferno),
        RenderConfig(strategy=HistogramRender(), zoom=10, colormap=:turbo, color_by=:absolute_frame),
        RenderConfig(strategy=CircleRender(), zoom=50, colormap=:turbo, color_by=:absolute_frame),
    ],
    outdir = outdir,
    verbose = Verbosity.DETAILED,
)

println("Output: $outdir")
println()

# File-based: analyze(config) loads data from DetectFitConfig.path
(result, analysis_info) = analyze(config)

# =============================================================================
# Summary
# =============================================================================
println("\n" * "="^60)
println("ANALYSIS COMPLETE")
println("="^60)

# Quick stats
if length(result.smld.emitters) > 0
    emitters = result.smld.emitters
    photons = [e.photons for e in emitters]
    σ_x = [e.σ_x for e in emitters]
    σ_y = [e.σ_y for e in emitters]

    println("\nFinal localizations:")
    println("  Count: $(length(emitters))")
    println("  Photons: median=$(round(median(photons), digits=0))")
    println("  Precision: σ_x=$(round(median(σ_x)*1000, digits=1))nm, σ_y=$(round(median(σ_y)*1000, digits=1))nm")
end

println("\nOutput directory: $outdir")
