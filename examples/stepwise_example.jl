"""
    stepwise_example.jl

SMLM analysis using step-by-step pure function calls.

This example demonstrates:
1. Generating realistic simulated SMLM data (4 datasets x 2000 frames)
2. Calling each pipeline step as a pure function
3. Threading SMLD state through the pipeline manually
4. Saving/loading intermediate SMLDs for iteration

Pipeline steps:
1. detectfit - combined detection and fitting
2. filter_step - quality filtering
3. frameconnect_step - link localizations across frames
4. driftcorrect_step - correct sample drift
5. render_step - generate super-resolution image

The stepwise approach is best for:
- Interactive parameter tuning
- Exploring different analysis strategies
- Development and debugging

For production workflows, see analysis_example.jl
"""

import Pkg
Pkg.activate(@__DIR__)

using SMLMAnalysis
using Statistics

include("generate_data.jl")

# ============================================================================
# Configuration
# ============================================================================

const OUTPUT_DIR = joinpath(@__DIR__, "output", "stepwise_example")

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

mkpath(OUTPUT_DIR)

# ============================================================================
# Step 1: Detection + Fitting
# ============================================================================

println("="^60)
println("Step 1: Detection + Fitting")
println("="^60)

(smld, df_info) = detectfit(image_stacks, camera, DetectFitConfig(
    boxsize = 7,
    min_photons = 500.0,
    psf_sigma = psf_sigma,
    backend = :cpu,
    psf_model = :variable,
    iterations = 20
); outdir=OUTPUT_DIR, step_number=1, verbose=Verbosity.STANDARD)

smld_raw = df_info.smld_raw

println("  Fitted emitters: $(length(smld.emitters))")
println("  Mean photons: $(round(mean([e.photons for e in smld.emitters]), digits=0))")
println("  Mean precision: $(round(mean([e.σ_x for e in smld.emitters])*1000, digits=1)) nm")
println()

# Save checkpoint after expensive detectfit
save_smld(joinpath(OUTPUT_DIR, "after_detectfit.h5"), smld)

# ============================================================================
# Step 2: Filtering
# ============================================================================

println("="^60)
println("Step 2: Filtering")
println("="^60)

n_before = length(smld.emitters)

(smld, f_info) = filter_step(smld, FilterConfig(
    photons = (500.0, Inf),
    precision = (0.0, 0.015),
    pvalue = (1e-3, 1.0)
); smld_raw=smld_raw, outdir=OUTPUT_DIR, step_number=2, verbose=Verbosity.STANDARD)

n_after = length(smld.emitters)
println("  Kept: $n_after / $n_before ($(round(100*n_after/n_before, digits=1))%)")
println()

# ============================================================================
# Step 3: Frame Connection
# ============================================================================

println("="^60)
println("Step 3: Frame Connection")
println("="^60)

(smld, fc_info) = frameconnect_step(smld, FrameConnectConfig(
    max_frame_gap = 5,
    max_sigma_dist = 5.0
); outdir=OUTPUT_DIR, step_number=3, verbose=Verbosity.STANDARD)

smld_connected = fc_info.smld_connected
n_tracks = length(unique([e.track_id for e in smld_connected.emitters]))
println("  Tracks: $n_tracks")
println()

# ============================================================================
# Step 4: Drift Correction
# ============================================================================

println("="^60)
println("Step 4: Drift Correction")
println("="^60)

(smld, dc_info) = driftcorrect_step(smld, DriftCorrectConfig(
    degree = 2,
    continuous = false
); outdir=OUTPUT_DIR, step_number=4, verbose=Verbosity.STANDARD)

drift_model = dc_info.drift_model
println()

# ============================================================================
# Step 5: Render
# ============================================================================

println("="^60)
println("Step 5: Render")
println("="^60)

(_, r_info) = render_step(smld, RenderConfig(zoom=20, colormap=:inferno);
    outdir=OUTPUT_DIR, step_number=5, verbose=Verbosity.STANDARD)

println()

# ============================================================================
# Build AnalysisConfig from the steps we ran (for reproducibility)
# ============================================================================

println("="^60)
println("Building AnalysisConfig for reproducibility")
println("="^60)

config = AnalysisConfig(
    camera = camera,
    steps = [
        DetectFitConfig(boxsize=7, min_photons=500.0, psf_sigma=psf_sigma,
                        backend=:cpu, psf_model=:variable, iterations=20),
        FilterConfig(photons=(500.0, Inf), precision=(0.0, 0.015), pvalue=(1e-3, 1.0)),
        FrameConnectConfig(max_frame_gap=5, max_sigma_dist=5.0),
        DriftCorrectConfig(degree=2, continuous=false),
        RenderConfig(zoom=20, colormap=:inferno),
    ],
    outdir = OUTPUT_DIR,
)
println("  Steps: $(length(config.steps))")
for (i, step) in enumerate(config.steps)
    println("    $i. $(step_name(step))")
end
println()

# ============================================================================
# Summary & Iteration Guide
# ============================================================================

println("="^60)
println("Analysis Complete")
println("="^60)
println()
println("Final results:")
println("  Localizations: $(length(smld.emitters))")
println()
println("Output: $OUTPUT_DIR")
println()
println("-"^60)
println("Iteration Guide")
println("-"^60)
println("""
To iterate on parameters (in REPL):

  # Load saved SMLD from after detectfit
  smld = load_smld("$(joinpath(OUTPUT_DIR, "after_detectfit.h5"))")

  # Try different filter settings
  (smld2, _) = filter_step(smld, FilterConfig(photons=(300.0, Inf), precision=(0.0, 0.020)))
  (smld2, _) = frameconnect_step(smld2, FrameConnectConfig(max_frame_gap=5))
  (smld2, _) = driftcorrect_step(smld2, DriftCorrectConfig(degree=2))
  (_, _) = render_step(smld2, RenderConfig(zoom=20, colormap=:inferno))

  # Or use the config for a full re-run:
  (result, info) = analyze(new_image_stacks, config)
""")
