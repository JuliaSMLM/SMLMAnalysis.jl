"""
    nanoruler_example.jl

Super-resolution analysis of simulated GATTAquant-style nanorulers.

This example demonstrates:
1. Simulating SMLM data of `Nanoruler2D` structures (SMLMSim) — n collinear,
   evenly-spaced marks per ruler, randomly placed/oriented in the field
   (4 datasets x 2000 frames)
2. Running the complete pipeline with AnalysisConfig
3. Rendering the reconstruction so the sub-diffraction marks resolve

Nanorulers are a standard SMLM calibration structure: a fixed, known spacing
between fluorophore marks (here 40 nm, 3 marks).

Acquisition is tuned to a realistic dSTORM-like blinking regime (see
generate_data.jl for the full target→rate derivation): ~2-frame blinks, ~1000
photons per saturated frame (→ ~4 nm precision, σ_psf/√photons), and ~5 blinks
per emitter. Rulers are sparse (2/μm², ~700 nm apart vs ~80 nm long) so each one
is cleanly isolated and the 40 nm marks resolve into a distinct 3-mark structure.

Analogous to analysis_example.jl (single-emitter Nmer2D), but with the
nanoruler pattern. For interactive parameter tuning, see stepwise_example.jl.
"""

import Pkg
Pkg.activate(@__DIR__)

using SMLMAnalysis

include("generate_data.jl")

# ============================================================================
# Configuration
# ============================================================================

const OUTPUT_DIR = joinpath(@__DIR__, "output", "nanoruler_example")

# ============================================================================
# Load or generate simulated data
# ============================================================================

println("="^60)
println("Loading nanoruler SMLM data")
println("="^60)

data = load_or_generate("nanoruler")

image_stacks = data["image_stacks"]  # Vector{SubArray}, dataset boundaries from data structure
camera = IdealCamera(data["camera_nx"], data["camera_ny"], data["camera_pixelsize"])
psf_sigma = data["psf_sigma"]

println("Camera: $(data["camera_nx"])x$(data["camera_ny"]) pixels, $(data["camera_pixelsize"]*1000)nm/pixel")
println("Data: $(length(image_stacks)) datasets x $(data["n_frames"]) frames/dataset")
println("Ruler: $(data["ruler_marks"]) marks, $(data["ruler_spacing"]*1000)nm spacing")
println()

# ============================================================================
# Run Analysis with AnalysisConfig
# ============================================================================

println("="^60)
println("Running analysis with AnalysisConfig")
println("="^60)

rm(OUTPUT_DIR; force=true, recursive=true)
mkpath(OUTPUT_DIR)

config = AnalysisConfig(
    camera = camera,
    steps = [
        DetectFitConfig(
            boxer = BoxerConfig(boxsize=7, min_photons=500.0, psf_sigma=psf_sigma),
            fitter = GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20),
        ),
        FilterConfig(
            photons = (500.0, Inf),
            precision = (0.0, 0.007),
            pvalue = (1e-3, 1.0)
        ),
        FrameConnectConfig(max_frame_gap = 5, calibration=CalibrationConfig()),
        DriftConfig(degree = 2),
        # NOTE: no DensityFilterConfig here. Nanorulers are sparse, isolated
        # structures (a few marks per ruler, rulers far apart), so a density
        # filter — which removes localizations with few neighbors — would erode
        # the very structure we want to see. It is appropriate for dense
        # clustered samples (see analysis_example.jl), not nanorulers.
        # High zoom so the 40 nm marks separate in the reconstruction.
        RenderConfig(zoom=40, colormap=:inferno, scalebar=true),
        RenderConfig(strategy=CircleRender(), zoom=80, colormap=:turbo,
                     color_by=:absolute_frame, scalebar=true),
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
println("Output saved to: $OUTPUT_DIR")
