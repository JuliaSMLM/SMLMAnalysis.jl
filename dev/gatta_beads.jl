# Gatta Beads Analysis - Full Pipeline (right-half ROI)
#
# Uses ROI to crop to the right half of the FOV.
# This is a CONTINUOUS dataset (single acquisition), so 1 image stack.
#
# Pipeline: detectfit -> filter -> frameconnect -> drift -> densityfilter -> render

using SMLMAnalysis

println("="^60)
println("GATTA BEADS ANALYSIS (right-half ROI)")
println("="^60)

# =============================================================================
# Data Path
# =============================================================================
h5file = joinpath(@__DIR__, "..", "data", "gatta_ruler", "2025-10-23",
    "20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5")

# Check file info
info = load_smart_h5_info(h5file)
println("Data: $(info.nframes) frames, $(info.width)x$(info.height), $(round(info.file_size_gb, digits=2)) GB")

# =============================================================================
# Load images
# =============================================================================
println("Loading images...")
t_load = @elapsed begin
    images, _ = smart_h5_to_array(h5file)
end
println("  Loaded $(size(images)) in $(round(t_load, digits=1))s")

# =============================================================================
# Camera Setup - ORCA-Fusion (C14440-20UP)
# =============================================================================
# 78nm pixels, 0.7 e-/ADU readnoise equivalent
# SCMOSCamera(nx, ny, ...) where nx=cols, ny=rows
# After smart_h5_to_array: data is (rows, cols, frames) = (height, width, frames)
camera = SCMOSCamera(info.width, info.height, 0.078f0, 0.7f0;
    offset = 100.0f0, gain = 0.24f0, qe = 1.0f0)
println("Camera: $(info.width)x$(info.height) (nxxny), 78nm pixels")

# =============================================================================
# Analysis Setup - ROI: right half of FOV
# =============================================================================
# x = columns (width), y = rows (height)
roi = (x = info.width÷2+1:info.width, y = 1:info.height)
println("ROI: x=$(roi.x), y=$(roi.y) (right half)")

outdir = joinpath(@__DIR__, "output", "gatta_beads")

# Single continuous acquisition = 1-element Vector
image_stacks = [images]

config = AnalysisConfig(
    camera = camera,
    steps = [
        DetectFitConfig(
            boxer = BoxerConfig(boxsize=11, min_photons=1000.0, psf_sigma=0.135),
            fitter = GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20),
        ),
        FilterConfig(
            photons = (500.0, Inf),
            precision = (0.0, 0.007),
            pvalue = (1e-6, 1.0),
            psf_sigma = :auto
        ),
        FrameConnectConfig(
            max_frame_gap = 5,
            max_sigma_dist = 5.0,
            calibration = CalibrationConfig(),
        ),
        DriftConfig(
            degree = 3,
            dataset_mode = :continuous,
            n_chunks = 4,
            quality = :iterative
        ),
        DensityFilterConfig(
            n_sigma = 2.0,
            min_neighbors = :auto
        ),
        RenderConfig(zoom=20, colormap=:inferno, clip_percentile=0.999),
        RenderConfig(strategy=HistogramRender(), zoom=10, colormap=:turbo, color_by=:absolute_frame),
        RenderConfig(strategy=CircleRender(), zoom=50, colormap=:turbo, color_by=:absolute_frame),
    ],
    roi = roi,
    outdir = outdir,
    verbose = Verbosity.DETAILED,
)

println("Output: $outdir")
println()

(result, analysis_info) = analyze(image_stacks, config)

# =============================================================================
# Summary
# =============================================================================
println("\n" * "="^60)
println("ANALYSIS COMPLETE")
println("="^60)

# Quick stats
if length(result.smld.emitters) > 0
    using Statistics
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
