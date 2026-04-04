"""
    bagol_example.jl

BaGoL (Bayesian Grouping of Localizations) demonstration with simulated DNA-PAINT data.

Simulates a hexameric pattern (6-mer, 40nm diameter) with DNA-PAINT-like kinetics
that produce ~10-20 localizations per binding site. BaGoL groups these scattered
localizations into precise emitter positions.

Pipeline:
1. Simulate DNA-PAINT data (6-mer clusters, 4000 frames)
2. DetectFit — detect and fit single-molecule events
3. Filter — quality filtering
4. FrameConnect — link blinking events + calibrate uncertainties
5-7. Render — 3 pre-BaGoL renders (Gaussian, Histogram, Circle)
8. BaGoL — Bayesian grouping → MAP-N emitters + diagnostics
         BaGoL step folder gets: posterior, NN distances, count distribution,
         acceptance rates, partition circles, overlay (locs + MAP-N)
9. Render — post-BaGoL Gaussian MAP-N (precision-weighted)

Requires: julia -t auto --project=. bagol_example.jl
"""

import Pkg
Pkg.activate(@__DIR__)

using SMLMAnalysis
using MicroscopePSFs
using Statistics

# ============================================================================
# Simulation: DNA-PAINT hexamer clusters
# ============================================================================

println("="^60)
println("Generating DNA-PAINT simulation")
println("="^60)

camera = IdealCamera(128, 128, 0.1)  # 128x128 px, 100nm/px
psf_sigma = 0.130  # 130nm PSF width

sim_params = StaticSMLMConfig(
    density = 1.0,       # 1 cluster/μm² — sparse enough for clean BaGoL
    σ_psf = psf_sigma,
    nframes = 4000,      # Long acquisition for ~10-20 locs/site
    ndatasets = 1,
)

# Hexamer: 6 binding sites on a 40nm-diameter circle
pattern = Nmer2D(n=6, d=0.040)

# DNA-PAINT kinetics: short binding events, moderate rebinding rate
fluor = GenericFluor(photons=50000.0, k_off=40.0, k_on=0.08)

t_sim = @elapsed begin
    (_, sim_info) = simulate(sim_params; pattern=pattern, molecule=fluor, camera=camera)
    smld_model = sim_info.smld_model

    psf = MicroscopePSFs.GaussianPSF(psf_sigma)
    (images, _) = gen_images(smld_model, psf; bg=20.0, poisson_noise=true)
    image_stacks = [images]  # Single dataset
end

n_true = sim_info.n_emitters
n_locs = sim_info.n_localizations
n_patterns = sim_info.n_patterns
println("  Clusters: $n_patterns ($(n_patterns * 6) binding sites)")
println("  True localizations: $n_locs (~$(round(n_locs / max(1, n_true), digits=1)) per site)")
println("  Image: $(size(images)) ($(round(t_sim, digits=1))s)")
println()

# ============================================================================
# Pipeline
# ============================================================================

const OUTPUT_DIR = joinpath(@__DIR__, "output", "bagol_example")
rm(OUTPUT_DIR; force=true, recursive=true)
mkpath(OUTPUT_DIR)

# --- Step 1: DetectFit ---
println("="^60)
println("Step 1: Detection + Fitting")
println("="^60)

(smld, df_info) = analyze(image_stacks, DetectFitConfig(
    camera = camera,
    boxer = BoxerConfig(boxsize=7, min_photons=300.0, psf_sigma=psf_sigma, backend=:auto),
    fitter = GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20, backend=:auto),
); outdir=OUTPUT_DIR, step_number=1, verbose=Verbosity.STANDARD)

println("  Fitted: $(length(smld.emitters)) localizations")
println("  Mean photons: $(round(mean([e.photons for e in smld.emitters]), digits=0))")
println("  Mean precision: $(round(mean([e.σ_x for e in smld.emitters])*1000, digits=1)) nm")
println()

# --- Step 2: Filter ---
println("="^60)
println("Step 2: Quality Filtering")
println("="^60)

n_before = length(smld.emitters)

(smld, f_info) = analyze(smld, FilterConfig(
    photons = (300.0, Inf),
    precision = (0.0, 0.020),
    pvalue = (1e-3, 1.0),
); outdir=OUTPUT_DIR, step_number=2, verbose=Verbosity.STANDARD)

println("  Kept: $(length(smld.emitters)) / $n_before")
println()

# --- Step 3: FrameConnect + Calibration ---
println("="^60)
println("Step 3: Frame Connection + Calibration")
println("="^60)

(smld, fc_info) = analyze(smld, FrameConnectConfig(
    max_frame_gap = 5,
    max_sigma_dist = 5.0,
    calibration = CalibrationConfig(clamp_k_to_one=true),
); outdir=OUTPUT_DIR, step_number=3, verbose=Verbosity.STANDARD)

cal = fc_info.info.calibration
if cal !== nothing && cal.calibration_applied
    println("  Calibration: k=$(round(cal.k_scale, digits=2)), σ_motion=$(round(cal.sigma_motion_nm, digits=1))nm")
end
println("  Localizations after connection: $(length(smld.emitters))")
println()

# --- Steps 5-7: Pre-BaGoL Renders ---
println("="^60)
println("Steps 5-7: Pre-BaGoL Renders")
println("="^60)

(_, _) = analyze(smld, RenderConfig(zoom=20, colormap=:inferno, scalebar=true);
    outdir=OUTPUT_DIR, step_number=5, verbose=Verbosity.STANDARD)

(_, _) = analyze(smld, RenderConfig(strategy=HistogramRender(), zoom=10, colormap=:turbo,
    color_by=:absolute_frame, clip_percentile=nothing, scalebar=true);
    outdir=OUTPUT_DIR, step_number=6, verbose=Verbosity.STANDARD)

(_, _) = analyze(smld, RenderConfig(strategy=CircleRender(), zoom=50, colormap=:turbo,
    color_by=:absolute_frame, scalebar=true);
    outdir=OUTPUT_DIR, step_number=7, verbose=Verbosity.STANDARD)
println()

# --- Step 8: BaGoL ---
println("="^60)
println("Step 8: BaGoL Grouping")
println("="^60)

(smld_bagol, bagol_step_info) = analyze(smld, BaGoLConfig(
    μ = 10.0,                    # Expected ~10 locs per binding site
    shape = 2.0,                 # DNA-PAINT: peaked count distribution
    learn_distribution = true,   # Let MCMC refine μ and shape
    n_iterations = 10000,
    burn_in = 2000,
    partition_sigma = 3.0,
    posterior_pixel_size = 0.002, # 2nm posterior image
); outdir=OUTPUT_DIR, step_number=8, verbose=Verbosity.STANDARD)

bagol_info = bagol_step_info.info
println("  Input: $(bagol_info.n_locs_in) localizations")
println("  Output: $(bagol_info.n_emitters) emitters")
println("  Compression: $(bagol_info.compression)x")
println("  Final μ: $(round(bagol_info.final_μ, digits=1)) locs/emitter")
println("  Final shape: $(round(bagol_info.final_shape, digits=2))")
println()

# --- Step 9: Post-BaGoL Gaussian Render ---
println("="^60)
println("Step 9: Post-BaGoL Render (Gaussian MAP-N)")
println("="^60)

(_, _) = analyze(smld_bagol, RenderConfig(zoom=50, colormap=:inferno, scalebar=true);
    outdir=OUTPUT_DIR, step_number=9, verbose=Verbosity.STANDARD)
println()

# ============================================================================
# Summary
# ============================================================================

println("="^60)
println("BaGoL Analysis Complete")
println("="^60)
println()

# Precision comparison
σ_pre = mean(sqrt(e.σ_x^2 + e.σ_y^2) for e in smld.emitters) * 1000
σ_post = if !isempty(smld_bagol.emitters)
    mean(sqrt(e.σ_x^2 + e.σ_y^2) for e in smld_bagol.emitters) * 1000
else
    NaN
end

println("  Pre-BaGoL:  $(length(smld.emitters)) localizations, σ=$(round(σ_pre, digits=1)) nm")
println("  Post-BaGoL: $(bagol_info.n_emitters) emitters, σ=$(round(σ_post, digits=1)) nm")
if isfinite(σ_post) && σ_post > 0
    println("  Precision improvement: $(round(σ_pre / σ_post, digits=1))x")
end
println()
println("  Output: $OUTPUT_DIR")
println()
println("  Key outputs:")
println("    05_render/          — Pre-BaGoL Gaussian")
println("    06_render/          — Pre-BaGoL Histogram")
println("    07_render/          — Pre-BaGoL Circle")
println("    08_bagol/           — BaGoL diagnostics, partition circles, overlay")
println("    09_render/          — Post-BaGoL Gaussian (MAP-N)")
