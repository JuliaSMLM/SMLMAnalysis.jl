"""
    plot_refinement_workflow.jl

Iterative workflow for refining analysis step plots using simulated data.

This workflow uses PERSISTENT CHECKPOINTS - you can close Julia and resume later!

Run this script section-by-section in the REPL to iterate on plot parameters:
1. Run once: Simulation setup (generates & saves data)
2. Run once: Analysis creation (with checkpoint=true for persistence)
3. Run repeatedly: Individual steps (tweak params, re-run)
4. Use reset!(a, N) to go back to checkpoint N and try different parameters
5. Use resume_analysis(OUTPUT_DIR) to continue from a previous session

Output structure:
    output/plot_refinement/
    ├── checkpoints/         # Persistent checkpoints (JLD2 files)
    │   ├── step_001_detect.jld2
    │   ├── step_002_fit.jld2
    │   └── ...
    ├── simulation/          # Saved images + ground truth
    ├── 01_detect/           # Detection outputs
    ├── 02_fit/              # Fit quality plots
    ├── 03_filter/           # Filter stats
    ├── 04_frameconnect/     # Track/calibration plots
    ├── 05_driftcorrect/     # Drift trajectory
    ├── 06_isolated/         # Neighbor histogram
    └── 07_render/           # Final renders
"""

# ============================================================================
# SECTION 1: Setup & Activate Environment
# ============================================================================
import Pkg
Pkg.activate(@__DIR__)

using SMLMAnalysis
using MicroscopePSFs
using JLD2
using Statistics

# Output directory
const OUTPUT_DIR = joinpath(@__DIR__, "output", "plot_refinement")
const SIM_DIR = joinpath(OUTPUT_DIR, "simulation")

# ============================================================================
# SECTION 2: Simulation (Run Once)
# ============================================================================
# This section generates simulated SMLM data and saves it to disk.
# Skip this section on subsequent runs by loading from disk instead.

println("="^70)
println("SECTION 2: Simulation Setup")
println("="^70)

# Check if simulation already exists
sim_file = joinpath(SIM_DIR, "sim_data.jld2")
if isfile(sim_file)
    println("Simulation data found at: $sim_file")
    println("To regenerate, delete this file and re-run.")
else
    println("Generating new simulation data...")
    mkpath(SIM_DIR)

    # Camera: 128x128 @ 100nm pixels = 12.8 x 12.8 um FOV
    camera = IdealCamera(128, 128, 0.1)
    println("  Camera: 128x128 pixels, 100nm/pixel, FOV = 12.8 x 12.8 um")

    # Fluorophore: realistic emission, sparse activation
    fluor = GenericFluor(
        photons = 50000.0,      # 50k photons/sec → ~1000 photons/frame @ 20fps
        k_off = 20.0,           # 20 Hz off-rate (50ms on-time)
        k_on = 0.02             # ~8 blinks over 400s acquisition (4 datasets × 100s)
    )
    println("  Fluorophore: 50k photons/s, k_off=20Hz, k_on=0.02Hz")

    # === MIXED SIMULATION: Hexamers + Random Background ===
    # This creates bimodal neighbor distribution for testing isolated filter

    # 1. Hexamer clusters: 6-mer @ 50nm, density 2 patterns/μm²
    hex_params = StaticSMLMParams(
        density = 2.0,          # 2 patterns/um^2 × 6 emitters = 12 emitters/μm²
        σ_psf = 0.13,
        nframes = 2000,
        ndatasets = 4
    )
    hex_pattern = Nmer2D(n=6, d=0.05)
    println("  Hexamers: 6-mer @ 50nm, density=2 patterns/μm²")

    # 2. Random background: monomers at density 6/μm² (isolated emitters)
    rand_params = StaticSMLMParams(
        density = 6.0,          # 6 emitters/μm² (using n=1 monomer) - higher to create distinct population
        σ_psf = 0.13,
        nframes = 2000,
        ndatasets = 4
    )
    rand_pattern = Nmer2D(n=1, d=0.0)  # Monomer = random point
    println("  Random: monomers, density=6 emitters/μm²")

    # Run both simulations
    t_sim = @elapsed begin
        hex_true, hex_model, hex_noisy = simulate(hex_params;
            pattern=hex_pattern, molecule=fluor, camera=camera)
        rand_true, rand_model, rand_noisy = simulate(rand_params;
            pattern=rand_pattern, molecule=fluor, camera=camera)
    end
    println("  Simulation time: $(round(t_sim, digits=2))s")
    println("  Hexamer emitters: $(length(hex_model.emitters))")
    println("  Random emitters: $(length(rand_model.emitters))")

    # Merge SMLDs (offset track_ids to avoid collision)
    function merge_smlds(smld1, smld2)
        offset = isempty(smld1.emitters) ? 0 : maximum(e.track_id for e in smld1.emitters)
        # Copy emitters from smld2 with offset track_id
        # Emitter2DFit fields: x, y, photons, bg, σ_x, σ_y, σ_photons, σ_bg, frame, dataset, track_id, id
        emitters2 = map(smld2.emitters) do e
            Emitter2DFit(e.x, e.y, e.photons, e.bg, e.σ_x, e.σ_y, e.σ_photons, e.σ_bg,
                         e.frame, e.dataset, e.track_id + offset, e.id)
        end
        BasicSMLD(vcat(smld1.emitters, emitters2), smld1.camera,
                  smld1.n_frames, smld1.n_datasets, smld1.metadata)
    end

    smld_true = merge_smlds(hex_true, rand_true)
    smld_model = merge_smlds(hex_model, rand_model)
    smld_noisy = merge_smlds(hex_noisy, rand_noisy)
    sim_params = hex_params  # Use hex params for metadata

    println("  Total emitters for imaging: $(length(smld_model.emitters))")

    # Generate camera images for ALL datasets
    # Use smld_model (true positions + photon counts), NOT smld_noisy (already has coordinate noise)
    # Using smld_noisy would double-apply noise and break CRLB validation
    psf = MicroscopePSFs.GaussianPSF(sim_params.σ_psf)
    t_img = @elapsed begin
        # gen_images defaults to dataset=1, need to call for each dataset and concatenate
        image_stacks = [gen_images(smld_model, psf;
            dataset=d,
            bg=10.0,
            poisson_noise=true,
            camera_noise=false
        ) for d in 1:sim_params.ndatasets]
        images = cat(image_stacks...; dims=3)
    end
    println("  Image generation time: $(round(t_img, digits=2))s")
    println("  Image stack: $(size(images)) ($(sim_params.ndatasets) datasets x $(sim_params.nframes) frames)")

    # Save simulation data
    @save sim_file images camera smld_true smld_model smld_noisy sim_params
    println("  Saved to: $sim_file")
end
println()

# ============================================================================
# SECTION 3: Load Simulation & Create Analysis (or Resume)
# ============================================================================
println("="^70)
println("SECTION 3: Load Simulation & Create Analysis")
println("="^70)

# Load simulation data
@load sim_file images camera smld_true smld_model smld_noisy sim_params
n_pixels_x = length(camera.pixel_edges_x) - 1
n_pixels_y = length(camera.pixel_edges_y) - 1
pixel_size = camera.pixel_edges_x[2] - camera.pixel_edges_x[1]
println("Loaded simulation data:")
println("  Images: $(size(images))")
println("  Camera: $(n_pixels_x)x$(n_pixels_y), $(pixel_size) um/pixel")
println("  Emitters in images (smld_model): $(length(smld_model.emitters))")
println()

# Check for existing checkpoints - offer to resume
checkpoint_dir = joinpath(OUTPUT_DIR, "checkpoints")
if isdir(checkpoint_dir) && !isempty(readdir(checkpoint_dir))
    println("Found existing checkpoints in $checkpoint_dir")
    println("To resume from previous session, run:")
    println("  a = resume_analysis(OUTPUT_DIR; images=images)")
    println("To start fresh, delete the checkpoints/ directory first.")
    println()
end

# Create Analysis object with DETAILED verbosity and CHECKPOINTING ENABLED
a = Analysis(images, camera;
    n_datasets = sim_params.ndatasets,  # Tell Analysis about multi-dataset structure
    outdir = OUTPUT_DIR,
    verbose = Verbosity.DETAILED,
    checkpoint = true  # Enable persistent checkpoints!
)
println("Analysis created with DETAILED verbosity")
println("Multi-dataset: $(a.n_datasets) datasets x $(a.n_frames_per_dataset) frames")
println("Checkpointing: ENABLED (saves to $checkpoint_dir)")
println("Output directory: $OUTPUT_DIR")
println()

# ============================================================================
# SECTION 4: Detection Step
# ============================================================================
# Key params to tweak:
#   - boxsize: ROI size (7, 9, 11)
#   - min_photons: detection threshold (lower = more detections)
#   - psf_sigma: expected PSF width
#
# Plots generated (DETAILED):
#   - detection_overlay.png: frames with detection boxes
#   - rois_per_frame.png: temporal detection rate

println("="^70)
println("SECTION 4: Detection")
println("="^70)

run_step!(a, DetectConfig(
    boxsize = 7,              # Compact ROIs for 130nm PSF @ 100nm pixels
    min_photons = 500.0,      # Lower threshold catches dimmer emitters
    psf_sigma = 0.13,         # Match simulation PSF
    use_gpu = false           # CPU for portability
))
println()

# ============================================================================
# SECTION 5: Fitting Step
# ============================================================================
# Key params to tweak:
#   - psf_model: :fixed (fastest), :variable (fit sigma), :anisotropic (fit sx/sy)
#   - iterations: MLE iterations (20 usually sufficient)
#
# Plots generated (DETAILED):
#   - fit_quality.png: photon/bg/precision/pvalue distributions
#   - fit_overlay.png: boxes colored by fit status
#   - precision_vs_photons.png: diagnostic scatter

println("="^70)
println("SECTION 5: Fitting")
println("="^70)

run_step!(a, FitConfig(
    psf_model = :variable,     # :fixed (fastest), :variable (fit sigma), :anisotropic (fit sx/sy)
    iterations = 20,
    device = nothing           # Auto-select GPU if available
))
println()

# ============================================================================
# SECTION 6: Filtering Step
# ============================================================================
# All filters use (min, max) tuples. Use Inf/-Inf for unbounded.
# Key params to tweak:
#   - photons: (min, max) - reject outside range (high can be doubles)
#   - precision: (min, max) - in microns (high = poor localization)
#   - pvalue: (min, max) - reject poor fit quality (low pvalue = bad fit)
#   - psf_sigma: :auto (mode ± 10%) or (min, max) tuple in microns
#
# Plots generated (DETAILED):
#   - detailed_stats.md: per-criterion rejection breakdown

println("="^70)
println("SECTION 6: Filtering")
println("="^70)

run_step!(a, FilterConfig(
    photons = (500.0, Inf),    # min 500, no max
    precision = (0.0, 0.015),  # max 15nm precision
    pvalue = (1e-3, 1.0),      # min 0.001 pvalue
    psf_sigma = :auto          # mode ± 10% (or specify (min, max) tuple in microns)
))
println()

# ============================================================================
# SECTION 7: Frame Connection Step
# ============================================================================
# Key params to tweak:
#   - maxframegap: max frames between appearances in same track
#   - nsigmadev: position matching tolerance (in sigma)
#   - calibrate: run uncertainty calibration from track statistics
#
# Plots generated (STANDARD+):
#   - track_histogram.png: localizations per track
#   - chi2_histogram.png: uncertainty validation
#   - uncertainty_calibration.png: variance calibration plot

println("="^70)
println("SECTION 7: Frame Connection")
println("="^70)

run_step!(a, FrameConnectConfig(
    maxframegap = 5,
    nsigmadev = 5.0,
    calibrate = true
))
println()

# ============================================================================
# SECTION 8: Drift Correction Step
# ============================================================================
# Key params to tweak:
#   - degree: polynomial degree (2=quadratic, 3=cubic)
#   - intramodel: "Polynomial" or "LegendrePoly"
#   - cost_fun: "Kdtree" (fast) or "Entropy" (robust)
#
# Plots generated (STANDARD+):
#   - drift_trajectory.png: X/Y drift over frames
# Plots generated (DETAILED):
#   - per_dataset.md: per-dataset max drift table

println("="^70)
println("SECTION 8: Drift Correction")
println("="^70)

run_step!(a, DriftCorrectConfig(
    degree = 2,                # Quadratic fit
    intramodel = "LegendrePoly",  # Better numerical stability (fixed frame number bug)
    cost_fun = "Kdtree",
    continuous = true          # TYPE 1: continuous acquisition (no registration between datasets)
))
println()

# ============================================================================
# SECTION 9: Isolated Emitter Filter Step
# ============================================================================
# Key params to tweak:
#   - n_sigma: neighbor search radius in multiples of combined sigma
#   - min_neighbors: :auto (triangle method) or explicit integer
#
# Plots generated (STANDARD+):
#   - neighbor_histogram.png: neighbor count distribution with threshold

println("="^70)
println("SECTION 9: Isolated Filter")
println("="^70)

run_step!(a, IsolatedConfig(
    n_sigma = 2.0,
    min_neighbors = :auto
))
println()

# ============================================================================
# SECTION 10: Render Step
# ============================================================================
# Default renders (at STANDARD verbosity):
#   1. Gaussian inferno 20x (density)
#   2. Histogram viridis 10x (time-encoded)
#   3. Circles viridis 50x (time-encoded)
#
# Custom renders via RenderSpec:
#   run_step!(a, RenderConfig(renders=[
#       RenderSpec(strategy=:gaussian, zoom=30, colormap=:turbo),
#       RenderSpec(strategy=:histogram, zoom=20, colormap=:inferno, color_by=:photons),
#   ]))

println("="^70)
println("SECTION 10: Render")
println("="^70)

run_step!(a, RenderConfig())  # Uses DEFAULT_RENDERS
println()

# ============================================================================
# SECTION 11: Ground Truth Comparison
# ============================================================================
# Render the ground truth simulation data for comparison

println("="^70)
println("SECTION 11: Ground Truth Comparison")
println("="^70)

gt_dir = joinpath(OUTPUT_DIR, "ground_truth")
mkpath(gt_dir)

# Render ground truth
println("Rendering ground truth (Gaussian, inferno)...")
render(smld_noisy;
    strategy = GaussianRender(),
    zoom = 20,
    colormap = :inferno,
    filename = joinpath(gt_dir, "ground_truth_gaussian.png")
)

render(smld_noisy;
    strategy = GaussianRender(),
    zoom = 20,
    colormap = :viridis,
    filename = joinpath(gt_dir, "ground_truth_viridis.png")
)
println("Ground truth renders saved to: $gt_dir")
println()

# ============================================================================
# SECTION 12: Summary & Iteration Guide
# ============================================================================
println("="^70)
println("WORKFLOW COMPLETE")
println("="^70)
println()
println("Analysis summary:")
println(a)
println()

println("="^70)
println("ITERATION GUIDE")
println("="^70)
println("""
To iterate on a specific step (SAME SESSION):

1. Reset to a checkpoint:
   reset!(a, 2)  # Go back to after fit (step 2)

2. Re-run with different parameters:
   run_step!(a, FilterConfig(min_photons=300))  # Try looser filter

3. Continue pipeline from there:
   run_step!(a, FrameConnectConfig(...))
   run_step!(a, DriftCorrectConfig(...))
   # etc.

Available checkpoints (memory): $(sort(collect(keys(a.checkpoints))))

RESUME FROM PREVIOUS SESSION:

Checkpoints are persisted to disk! To continue after closing Julia:

1. Load simulation data:
   @load sim_file images camera smld_true smld_model smld_noisy sim_params

2. Resume analysis:
   a = resume_analysis(OUTPUT_DIR; images=images)

3. Reset to any step and continue:
   reset!(a, 2)  # Loads from disk if not in memory
   run_step!(a, FilterConfig(min_photons=300))

Checkpoint files: $(checkpoint_dir)/

Tips:
- Checkpoints are saved after every step (not just expensive ones)
- reset!() works across sessions - loads from disk automatically
- To start completely fresh, delete the checkpoints/ directory
- Compare renders in output/plot_refinement/07_render/ vs ground_truth/
- Check 01_detect/rois_per_frame.png for detection stability
- Check 02_fit/precision_vs_photons.png for fit quality trends
- Check 04_frameconnect/uncertainty_calibration.png for k and sigma_motion
""")
