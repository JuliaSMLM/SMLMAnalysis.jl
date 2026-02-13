"""
    analysis_example.jl

SMLM analysis using the one-liner analyze() function.

This example demonstrates:
1. Generating realistic simulated SMLM data (4 datasets x 2000 frames)
2. Running complete pipeline with AnalysisConfig (primary interface)

The analyze() approach is best for:
- Quick analysis with sensible defaults
- Production workflows with standard parameters
- Reproducible configs (save/share AnalysisConfig)

For interactive parameter tuning, see stepwise_example.jl
"""

import Pkg
Pkg.activate(@__DIR__)

using SMLMAnalysis

include("generate_data.jl")

# ============================================================================
# Configuration
# ============================================================================

const OUTPUT_DIR = joinpath(@__DIR__, "output", "analysis_example")

# ============================================================================
# Load or generate simulated data
# ============================================================================

println("="^60)
println("Loading SMLM data")
println("="^60)

data = load_or_generate("single_target")

image_stacks = data["image_stacks"]  # Vector{SubArray}, dataset boundaries from data structure
camera = IdealCamera(data["camera_nx"], data["camera_ny"], data["camera_pixelsize"])
psf_sigma = data["psf_sigma"]

println("Camera: $(data["camera_nx"])x$(data["camera_ny"]) pixels, $(data["camera_pixelsize"]*1000)nm/pixel")
println("Data: $(length(image_stacks)) datasets x $(data["n_frames"]) frames/dataset")
println()

# ============================================================================
# Run Analysis with AnalysisConfig (primary interface)
# ============================================================================

println("="^60)
println("Running analysis with AnalysisConfig")
println("="^60)

mkpath(OUTPUT_DIR)

config = AnalysisConfig(
    camera = camera,
    steps = [
        DetectFitConfig(
            boxer = BoxerConfig(boxsize=7, min_photons=500.0, psf_sigma=psf_sigma, backend=:cpu),
            fitter = GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20, backend=:cpu),
        ),
        FilterConfig(
            photons = (500.0, Inf),
            precision = (0.0, 0.007),
            pvalue = (1e-3, 1.0)
        ),
        FrameConnectConfig(max_frame_gap = 5),
        CalibrationConfig(),
        DriftConfig(degree = 2),
        DensityFilterConfig(n_sigma=2.0, min_neighbors=:auto),
        RenderConfig(zoom=20, colormap=:inferno),
        RenderConfig(strategy=HistogramRender(), zoom=10, colormap=:turbo, color_by=:absolute_frame),
        RenderConfig(strategy=CircleRender(), zoom=50, colormap=:turbo, color_by=:absolute_frame),
    ],
    outdir = OUTPUT_DIR,
    verbose = Verbosity.STANDARD
)

# analyze() returns (AnalysisResult, AnalysisInfo) tuple
(result, info) = analyze(image_stacks, config)

# ============================================================================
# Summary
# ============================================================================

println()
println("="^60)
println("Analysis Complete")
println("="^60)
println()
println("Results:")
println("  Final localizations: $(length(result.smld.emitters))")
println("  Datasets: $(result.smld.n_datasets)")
println("  Frames per dataset: $(result.smld.n_frames)")
println("  Total time: $(round(info.elapsed_s, digits=2))s")
println()
println("Step info available:")
for (name, _) in info.steps
    println("  info.steps[:$name]")
end
println()
println("Output saved to: $OUTPUT_DIR")
