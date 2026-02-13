# Hexabody dSTORM Analysis - Single Cell (Parameterized)
#
# Usage: julia -t auto --project=dev dev/genmab_hexabody_cell.jl <h5file> <outname>
#
# Standard registered-acquisition pipeline with degree=2 singlepass drift correction.

import Pkg
Pkg.activate(@__DIR__)

using SMLMAnalysis
using Statistics

length(ARGS) >= 2 || error("Usage: julia ... genmab_hexabody_cell.jl <h5file> <outname>")
h5file = ARGS[1]
outname = ARGS[2]

isfile(h5file) || error("File not found: $h5file")

println("="^60)
println("HEXABODY dSTORM - $outname")
println("="^60)

# Data info
info = load_mic_h5_info(h5file)
println("Data: $(info.n_frames) frames in $(info.n_blocks) datasets")
println("Image size: $(info.width)x$(info.height)")
println("File size: $(round(info.file_size_gb, digits=2)) GB")

# Camera setup - Hamamatsu sCMOS with per-pixel calibration
pixel_size = 0.0978f0
cal = load_mic_h5_calibration_for_scmos(h5file)
camera = SCMOSCamera(info.width, info.height, pixel_size, cal.readnoise;
    offset = cal.offset, gain = cal.gain, qe = 0.82f0)

# Analysis
outdir = joinpath(@__DIR__, "output", outname)

config = AnalysisConfig(
    camera = camera,
    steps = [
        DetectFitConfig(
            path = h5file,
            h5_format = :mic,
            boxer = BoxerConfig(boxsize=9, min_photons=500.0, psf_sigma=0.130),
            fitter = GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20),
        ),
        FilterConfig(
            photons = (500.0, Inf),
            precision = (0.0, 0.007),
            pvalue = (1e-6, 1.0)
        ),
        FrameConnectConfig(
            max_frame_gap = 5,
            max_sigma_dist = 5.0,
            calibration = CalibrationConfig(),
        ),
        DriftConfig(
            degree = 2,
            dataset_mode = :registered,
            quality = :iterative,
            auto_roi = false
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

println("Output: $outdir\n")

(result, analysis_info) = analyze(config)

# Summary
println("\n" * "="^60)
println("COMPLETE: $outname")
println("="^60)
if length(result.smld.emitters) > 0
    emitters = result.smld.emitters
    println("  Localizations: $(length(emitters))")
    println("  Photons: median=$(round(median([e.photons for e in emitters]), digits=0))")
    println("  Precision: $(round(median([e.σ_x for e in emitters])*1000, digits=1))nm")
end
