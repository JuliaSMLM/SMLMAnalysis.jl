"""
    stepwise_example.jl

SMLM analysis using step-by-step workflow with checkpointing.

This example demonstrates:
1. Generating realistic simulated SMLM data (4 datasets x 2000 frames)
2. Creating an Analysis object with checkpointing enabled
3. Running each step individually with run_step!()
4. Using reset!() to iterate on parameters

Pipeline steps:
1. DetectFitConfig - combined detection and fitting
2. FilterConfig - quality filtering
3. FrameConnectConfig - link localizations across frames
4. DriftCorrectConfig - correct sample drift
5. RenderConfig - generate super-resolution image

The stepwise approach is best for:
- Interactive parameter tuning
- Exploring different analysis strategies
- Development and debugging

For production workflows, see recipe_example.jl
"""

import Pkg
Pkg.activate(@__DIR__)

using SMLMAnalysis
using MicroscopePSFs
using Statistics

# ============================================================================
# Configuration
# ============================================================================

const OUTPUT_DIR = joinpath(@__DIR__, "output", "stepwise_example")

# Simulation parameters
const N_FRAMES = 2000      # Frames per dataset
const N_DATASETS = 4       # Number of datasets
const PIXEL_SIZE = 0.1     # 100 nm pixels
const PSF_SIGMA = 0.13     # 130 nm PSF width

# ============================================================================
# Helper: Generate images per dataset (workaround for SMLMSim bug)
# ============================================================================

"""
Generate images for a specific dataset by filtering emitters first.
Workaround for gen_images not respecting dataset parameter.
"""
function gen_images_for_dataset(smld, psf, dataset::Int; kwargs...)
    emitters_d = filter(e -> e.dataset == dataset, smld.emitters)
    smld_d = BasicSMLD(emitters_d, smld.camera, smld.n_frames, 1, smld.metadata)
    (images, _) = gen_images(smld_d, psf; dataset=1, kwargs...)
    images
end

# ============================================================================
# Step 1: Generate Simulated Data
# ============================================================================

println("="^60)
println("Generating simulated SMLM data")
println("="^60)

camera = IdealCamera(256, 128, PIXEL_SIZE)  # width=256, height=128 (rectangular to catch permutation bugs)

sim_params = StaticSMLMConfig(
    density = 2.0,
    σ_psf = PSF_SIGMA,
    nframes = N_FRAMES,
    ndatasets = N_DATASETS
)
pattern = Nmer2D(n=8, d=0.15)
fluor = GenericFluor(photons=50000.0, k_off=20.0, k_on=0.02)

println("Simulating $(N_DATASETS) datasets x $(N_FRAMES) frames/dataset...")
t_sim = @elapsed begin
    (_, sim_info) = simulate(sim_params; pattern=pattern, molecule=fluor, camera=camera)
    smld_truth = sim_info.smld_true
    smld_model = sim_info.smld_model
end
println("  Ground truth emitters: $(length(smld_truth.emitters))")
println("  Emitter appearances: $(length(smld_model.emitters)) ($(round(t_sim, digits=1))s)")

# Generate images per dataset (workaround for SMLMSim bug)
psf = MicroscopePSFs.GaussianPSF(PSF_SIGMA)
println("Generating images per dataset...")
t_img = @elapsed begin
    image_stacks = [gen_images_for_dataset(smld_model, psf, d; bg=10.0, poisson_noise=true)
                    for d in 1:N_DATASETS]
end
println("  $(N_DATASETS) stacks of $(size(image_stacks[1])) ($(round(t_img, digits=1))s)")

# Concatenate for Analysis (until we implement Vector{Array} dispatch)
images = cat(image_stacks...; dims=3)
println("  Combined: $(size(images))")
println()

# ============================================================================
# Step 2: Create Analysis Object
# ============================================================================

println("="^60)
println("Creating Analysis with checkpointing")
println("="^60)

mkpath(OUTPUT_DIR)

a = Analysis(images, camera;
    n_datasets = N_DATASETS,
    outdir = OUTPUT_DIR,
    verbose = Verbosity.STANDARD,
    checkpoint = true
)

println("Analysis created:")
println("  Datasets: $(a.n_datasets)")
println("  Frames per dataset: $(a.n_frames_per_dataset)")
println("  Output: $OUTPUT_DIR")
println()

# ============================================================================
# Step 3: Detection + Fitting
# ============================================================================

println("="^60)
println("Step 1: Detection + Fitting")
println("="^60)

run_step!(a, DetectFitConfig(
    boxsize = 7,
    min_photons = 500.0,
    psf_sigma = PSF_SIGMA,
    backend = :cpu,
    psf_model = :variable,
    iterations = 20
))

println("  Fitted emitters: $(length(a.smld.emitters))")
println("  Mean photons: $(round(mean([e.photons for e in a.smld.emitters]), digits=0))")
println("  Mean precision: $(round(mean([e.σ_x for e in a.smld.emitters])*1000, digits=1)) nm")
println()

# ============================================================================
# Step 4: Filtering
# ============================================================================

println("="^60)
println("Step 2: Filtering")
println("="^60)

n_before = length(a.smld.emitters)

run_step!(a, FilterConfig(
    photons = (500.0, Inf),
    precision = (0.0, 0.015),
    pvalue = (1e-3, 1.0)
))

n_after = length(a.smld.emitters)
println("  Kept: $n_after / $n_before ($(round(100*n_after/n_before, digits=1))%)")
println()

# ============================================================================
# Step 5: Frame Connection
# ============================================================================

println("="^60)
println("Step 3: Frame Connection")
println("="^60)

run_step!(a, FrameConnectConfig(
    maxframegap = 5,
    nsigmadev = 5.0
))

n_tracks = length(unique([e.track_id for e in a.smld_connected.emitters]))
println("  Tracks: $n_tracks")
println()

# ============================================================================
# Step 6: Drift Correction
# ============================================================================

println("="^60)
println("Step 4: Drift Correction")
println("="^60)

run_step!(a, DriftCorrectConfig(
    degree = 2,
    continuous = false
))

println("  Drift model datasets: $(a.drift_model.ndatasets)")
println()

# ============================================================================
# Step 7: Render
# ============================================================================

println("="^60)
println("Step 5: Render")
println("="^60)

run_step!(a, RenderConfig(zoom = 20))

println()

# ============================================================================
# Step 8: BaGoL (Bayesian Grouping)
# ============================================================================

println("="^60)
println("Step 6: BaGoL (Bayesian Grouping of Localizations)")
println("="^60)

run_step!(a, BaGoLConfig(
    n_iterations = 10000,
    burn_in = 2000,
    α = :auto,           # Estimate from data
    learn_α = true       # Update during MCMC
    # render_zoom = 50 (default)
))

if a.bagol_result !== nothing && a.bagol_smld !== nothing
    println("  Grouped emitters: $(a.bagol_result.n_emitters)")
    compression = length(a.smld.emitters) / max(1, a.bagol_result.n_emitters)
    println("  Compression: $(round(compression, digits=1))x")

    # Compare against ground truth
    n_truth = length(smld_truth.emitters)
    println()
    println("  Ground truth comparison:")
    println("    True emitters: $n_truth")
    println("    BaGoL found: $(a.bagol_result.n_emitters)")
    println("    Difference: $(a.bagol_result.n_emitters - n_truth)")

    # Render truth vs BaGoL overlay with ellipses
    # Green = ground truth positions, Red = BaGoL result
    render([smld_truth, a.bagol_smld];
        colors = [:green, :red],
        strategy = EllipseRender(),
        zoom = 50,
        filename = joinpath(OUTPUT_DIR, "06_bagol", "truth_vs_bagol.png")
    )
    println("    Truth vs BaGoL overlay saved (green=truth, red=BaGoL)")
end
println()

# ============================================================================
# Extract reproducible config from step history
# ============================================================================

println("="^60)
println("Extracting AnalysisConfig")
println("="^60)

config = get_config(a)
println("  Steps: $(length(config.steps))")
for (i, step) in enumerate(config.steps)
    println("    $i. $(step_name(step))")
end

# Extract AnalysisInfo from step records (tuple-pattern)
info = get_analysis_info(a)
println("  Total time: $(round(info.elapsed_s, digits=2))s")
println()

# ============================================================================
# Summary & Iteration Guide
# ============================================================================

println("="^60)
println("Analysis Complete")
println("="^60)
println()
println("Final results:")
println("  Localizations: $(length(a.smld.emitters))")
println("  Checkpoints: $(sort(collect(keys(a.checkpoints))))")
println()
println("Output: $OUTPUT_DIR")
println()
println("-"^60)
println("Iteration Guide")
println("-"^60)
println("""
To iterate on parameters (in REPL):

  # Reset to after detectfit (step 1) and try different filter
  reset!(a, 1)
  run_step!(a, FilterConfig(photons=(300.0, Inf), precision=(0.0, 0.020)))
  run_step!(a, FrameConnectConfig(maxframegap=5))
  run_step!(a, DriftCorrectConfig(degree=2))
  run_step!(a, RenderConfig(zoom=20))

  # Extract config for reproducibility
  config = get_config(a)
  (result, info) = analyze(new_images, config)

  # Resume from disk after closing Julia:
  a = resume_analysis(OUTPUT_DIR; images=images)
  reset!(a, 1)
""")
