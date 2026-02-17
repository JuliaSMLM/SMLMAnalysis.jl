"""
    stepwise_example.jl

SMLM analysis using step-by-step analyze() calls.

This example demonstrates:
1. Generating realistic simulated SMLM data (4 datasets x 2000 frames)
2. Calling analyze() with individual step configs
3. Threading SMLD state through the pipeline manually
4. Saving/loading intermediate SMLDs for iteration

Pipeline steps (all via analyze() dispatch):
1. analyze(images, DetectFitConfig(...)) - detection and fitting
2. analyze(smld, FilterConfig(...)) - quality filtering
3. analyze(smld, FrameConnectConfig(...)) - link across frames + calibrate uncertainties
4. analyze(smld, DriftConfig(...)) - correct sample drift
5. analyze(smld, RenderConfig(...)) - super-resolution image

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

rm(OUTPUT_DIR; force=true, recursive=true)
mkpath(OUTPUT_DIR)

# ============================================================================
# Step 1: Detection + Fitting
# ============================================================================

println("="^60)
println("Step 1: Detection + Fitting")
println("="^60)

(smld, df_info) = analyze(image_stacks, DetectFitConfig(
    camera = camera,
    boxer = BoxerConfig(boxsize=7, min_photons=500.0, psf_sigma=psf_sigma, backend=:auto),
    fitter = GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20, backend=:auto),
); outdir=OUTPUT_DIR, step_number=1, verbose=Verbosity.STANDARD)

smld_raw = smld  # Raw SMLD from detectfit, before filtering

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

(smld, f_info) = analyze(smld, FilterConfig(
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

(smld, fc_info) = analyze(smld, FrameConnectConfig(
    max_frame_gap = 5,
    max_sigma_dist = 5.0,
    calibration = CalibrationConfig(clamp_k_to_one=true)
); outdir=OUTPUT_DIR, step_number=3, verbose=Verbosity.STANDARD)

smld_connected = fc_info.info.connected  # FrameConnectInfo.connected

n_tracks = length(unique([e.track_id for e in smld_connected.emitters]))
println("  Tracks: $n_tracks")

# Calibration diagnostics from FrameConnectInfo
cal = fc_info.info.calibration
if cal !== nothing && cal.calibration_applied
    println("  Calibration: k=$(round(cal.k_scale, digits=2)), sigma_motion=$(round(cal.sigma_motion_nm, digits=1))nm")
end
println()

# ============================================================================
# Step 4: Drift Correction
# ============================================================================

println("="^60)
println("Step 4: Drift Correction")
println("="^60)

(smld, dc_info) = analyze(smld, DriftConfig(
    degree = 2,
    dataset_mode = :registered
); outdir=OUTPUT_DIR, step_number=4, verbose=Verbosity.STANDARD)

drift_model = dc_info.info.model  # DriftInfo.model
println()

# ============================================================================
# Step 5-7: Render (three strategies)
# ============================================================================

println("="^60)
println("Step 5: Render (Gaussian)")
println("="^60)

(_, r_info) = analyze(smld, RenderConfig(zoom=20, colormap=:inferno, scalebar=true);
    outdir=OUTPUT_DIR, step_number=5, verbose=Verbosity.STANDARD)

println()

println("="^60)
println("Step 6: Render (Histogram)")
println("="^60)

(_, _) = analyze(smld, RenderConfig(strategy=HistogramRender(), zoom=10, colormap=:turbo, color_by=:absolute_frame, clip_percentile=nothing, scalebar=true);
    outdir=OUTPUT_DIR, step_number=6, verbose=Verbosity.STANDARD)

println()

println("="^60)
println("Step 7: Render (Circle)")
println("="^60)

(_, _) = analyze(smld, RenderConfig(strategy=CircleRender(), zoom=50, colormap=:turbo, color_by=:absolute_frame, scalebar=true);
    outdir=OUTPUT_DIR, step_number=7, verbose=Verbosity.STANDARD)

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
        DetectFitConfig(
            boxer=BoxerConfig(boxsize=7, min_photons=500.0, psf_sigma=psf_sigma, backend=:auto),
            fitter=GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20, backend=:auto)),
        FilterConfig(photons=(500.0, Inf), precision=(0.0, 0.015), pvalue=(1e-3, 1.0)),
        FrameConnectConfig(max_frame_gap=5, max_sigma_dist=5.0,
            calibration=CalibrationConfig(clamp_k_to_one=true)),
        DriftConfig(degree=2, dataset_mode=:registered),
        RenderConfig(zoom=20, colormap=:inferno, scalebar=true),
        RenderConfig(strategy=HistogramRender(), zoom=10, colormap=:turbo, color_by=:absolute_frame, clip_percentile=nothing, scalebar=true),
        RenderConfig(strategy=CircleRender(), zoom=50, colormap=:turbo, color_by=:absolute_frame, scalebar=true),
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
  (smld2, _) = analyze(smld, FilterConfig(photons=(300.0, Inf), precision=(0.0, 0.020)))
  (smld2, _) = analyze(smld2, FrameConnectConfig(max_frame_gap=5))
  (smld2, _) = analyze(smld2, DriftConfig(degree=2))
  (_, _) = analyze(smld2, RenderConfig(zoom=20, colormap=:inferno, scalebar=true))

  # Or use the config for a full re-run:
  (result, info) = analyze(new_image_stacks, config)
""")
