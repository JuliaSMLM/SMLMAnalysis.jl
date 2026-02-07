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
using MicroscopePSFs

# ============================================================================
# Configuration
# ============================================================================

const OUTPUT_DIR = joinpath(@__DIR__, "output", "analysis_example")

# Simulation parameters - realistic SMLM acquisition
const N_FRAMES = 2000      # Frames per dataset
const N_DATASETS = 4       # Number of datasets (total: 8000 frames)
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

# Camera: 256x128 pixels @ 100nm = 25.6 x 12.8 um FOV (rectangular to catch permutation bugs)
camera = IdealCamera(256, 128, PIXEL_SIZE)
println("Camera: 256x128 pixels, $(PIXEL_SIZE*1000)nm/pixel")

# Simulation: octamer pattern (8-mer @ 150nm diameter)
sim_params = StaticSMLMConfig(
    density = 2.0,          # 2 patterns/um^2
    σ_psf = PSF_SIGMA,
    nframes = N_FRAMES,     # Frames per dataset
    ndatasets = N_DATASETS
)
pattern = Nmer2D(n=8, d=0.05)  # Octamer, 150nm diameter

# Fluorophore: realistic blinking kinetics
fluor = GenericFluor(
    photons = 50000.0,      # 50k photons/sec -> ~2500 photons/frame @ 20fps
    k_off = 20.0,           # 20 Hz off-rate (50ms on-time)
    k_on = 0.02             # Sparse activation
)

println("Simulation: $(N_DATASETS) datasets x $(N_FRAMES) frames/dataset")
println("Pattern: 8-mer @ 50nm, density = $(sim_params.density) patterns/um^2")

# Run simulation
t_sim = @elapsed begin
    (_, sim_info) = simulate(sim_params; pattern=pattern, molecule=fluor, camera=camera)
    smld_model = sim_info.smld_model
end
println("Simulated $(length(smld_model.emitters)) emitter appearances ($(round(t_sim, digits=1))s)")

# Generate camera images per dataset (workaround for SMLMSim bug)
psf = MicroscopePSFs.GaussianPSF(PSF_SIGMA)
println("Generating images per dataset...")
t_img = @elapsed begin
    image_stacks = [gen_images_for_dataset(smld_model, psf, d; bg=20.0, poisson_noise=true)
                    for d in 1:N_DATASETS]
end
println("  $(N_DATASETS) stacks of $(size(image_stacks[1])) ($(round(t_img, digits=1))s)")

# Concatenate for Analysis (until we implement Vector{Array} dispatch)
images = cat(image_stacks...; dims=3)
println("  Combined: $(size(images))")
println()

# ============================================================================
# Step 2: Run Analysis with AnalysisConfig (primary interface)
# ============================================================================

println("="^60)
println("Running analysis with AnalysisConfig")
println("="^60)

mkpath(OUTPUT_DIR)

config = AnalysisConfig(
    camera = camera,
    steps = [
        DetectFitConfig(
            boxsize = 7,
            min_photons = 500.0,
            psf_sigma = PSF_SIGMA,
            backend = :cpu,
            psf_model = :variable,
            iterations = 20,
            n_datasets = N_DATASETS
        ),
        FilterConfig(
            photons = (500.0, Inf),
            precision = (0.0, 0.007),
            pvalue = (1e-3, 1.0)
        ),
        FrameConnectConfig(maxframegap = 5),
        DriftCorrectConfig(degree = 2),
        DensityFilterConfig(n_sigma=2.0, min_neighbors=:auto),
        RenderConfig(zoom=20, colormap=:inferno),
        RenderConfig(strategy=HistogramRender(), zoom=10, colormap=:turbo, color_by=:absolute_frame),
        RenderConfig(strategy=CircleRender(), zoom=50, colormap=:turbo, color_by=:absolute_frame),
    ],
    outdir = OUTPUT_DIR,
    verbose = Verbosity.STANDARD
)

# analyze() returns (Analysis, AnalysisInfo) tuple
(result, info) = analyze(images, config)

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
