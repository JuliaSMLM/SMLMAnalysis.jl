# Ground-truth validation: double-emitter rate (p₂) estimation
# Simulates tight dimers (Nmer2D, n=2, d=50nm) and compares
# estimated p₂ from mixture decomposition against true double-on rate.

using SMLMAnalysis
using SMLMAnalysis: MicroscopePSFs
using Statistics

outdir = joinpath(@__DIR__, "output", "test_mixture_gt")
rm(outdir; force=true, recursive=true)
mkpath(outdir)

# --- Simulation parameters ---
pixel_size = 0.1  # μm
npixels = 128
σ_psf = 0.130     # μm
nframes = 10000
density = 0.8     # patterns per μm²

camera = IdealCamera(npixels, npixels, pixel_size)

# Higher duty cycle for statistical power: k_on=0.1 → d ≈ 0.5%
# With ~100 dimers × 10000 frames, expect ~25 double-on events
fluor = GenericFluor(photons=50000.0, k_off=20.0, k_on=0.1)

sim_params = StaticSMLMConfig(
    density=density,
    σ_psf=σ_psf,
    nframes=nframes,
    minphotons=500,
    ndatasets=1,
    framerate=50.0,
)

pattern = Nmer2D(n=2, d=0.05)  # 50nm dimer separation

println("Simulating dimers...")
t_sim = @elapsed (smld_noisy, sim_info) = simulate(sim_params;
    pattern=pattern,
    molecule=fluor,
    camera=camera,
)

smld_model = sim_info.smld_model
println("  $(sim_info.n_emitters) fluorophore positions ($(sim_info.n_patterns) dimers)")
println("  $(sim_info.n_localizations) localizations in $nframes frames ($(round(t_sim, digits=1))s)")
println("  smld_model: $(length(smld_model.emitters)) blink events")

# --- Ground-truth p₂ from smld_model ---
# Group by (dataset, frame, id) — id identifies the pattern instance
# A group with 2+ unique track_ids = both dimer partners ON simultaneously
println("\nComputing ground-truth p₂...")

events = smld_model.emitters
# Build groups: (dataset, frame, id) → Set of track_ids
groups = Dict{Tuple{Int,Int,Int}, Set{Int}}()
for e in events
    key = (e.dataset, e.frame, e.id)
    if !haskey(groups, key)
        groups[key] = Set{Int}()
    end
    push!(groups[key], e.track_id)
end

n_double = count(tids -> length(tids) >= 2, values(groups))
n_single = count(tids -> length(tids) == 1, values(groups))
n_total_events = n_single + n_double
p2_true = n_total_events > 0 ? n_double / n_total_events : 0.0

println("  Single-on events: $n_single")
println("  Double-on events: $n_double")
println("  Total events: $n_total_events")
println("  p₂_true = $(round(100 * p2_true, digits=3))%")

# --- Generate camera images ---
println("\nGenerating camera images...")
t_gen = @elapsed (images, img_info) = gen_images(smld_model, MicroscopePSFs.GaussianPSF(σ_psf);
    bg=5.0, poisson_noise=true)
println("  Image size: $(size(images)) ($(round(t_gen, digits=1))s)")

# --- Run pipeline: detectfit → filter → intensity filter ---
println("\nRunning analysis pipeline...")

config = AnalysisConfig(
    camera = camera,
    steps = [
        DetectFitConfig(
            boxer = BoxerConfig(boxsize=9, min_photons=500.0, psf_sigma=σ_psf),
            fitter = GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20),
        ),
        FilterConfig(
            photons = (500.0, Inf),
            precision = (0.0, 0.010),
            pvalue = (1e-6, 1.0),
        ),
        IntensityFilterConfig(
            field_mode = :uniform,
            estimate_p2 = true,
            p2_tail_threshold = 1.0,
            p2_n_bins = 200,
        ),
    ],
    outdir = outdir,
    verbose = Verbosity.DETAILED,
)

t_pipe = @elapsed (result, analysis_info) = analyze([images], config)

# --- Extract p₂ estimate ---
if_info = analysis_info.step_infos[end].info
p2_est = if_info.p2_estimate

println("\n" * "="^60)
println("RESULTS: Double-Emitter Rate Estimation")
println("="^60)
println("  Ground truth p₂:  $(round(100 * p2_true, digits=3))%  ($n_double / $n_total_events)")
println("  Estimated p₂:     $(p2_est !== nothing ? "$(round(100 * p2_est, digits=3))%" : "N/A")")
if p2_est !== nothing && p2_true > 0
    ratio = p2_est / p2_true
    println("  Ratio (est/true): $(round(ratio, digits=2))")
end
println("  Localizations:    $(length(result.smld.emitters))")
println("  Pipeline time:    $(round(t_pipe, digits=1))s")
println("="^60)

# Save summary
open(joinpath(outdir, "gt_comparison.md"), "w") do io
    println(io, "# Ground-Truth p₂ Comparison\n")
    println(io, "## Simulation")
    println(io, "- Patterns: $(sim_info.n_patterns) dimers (Nmer2D, d=50nm)")
    println(io, "- Frames: $nframes")
    println(io, "- Duty cycle: k_on=$(fluor.q[2,1]), k_off=$(fluor.q[1,2])")
    println(io, "- Blink events: $(length(smld_model.emitters))")
    println(io, "")
    println(io, "## Ground Truth")
    println(io, "- Single-on events: $n_single")
    println(io, "- Double-on events: $n_double")
    println(io, "- **p₂_true = $(round(100 * p2_true, digits=3))%**")
    println(io, "")
    println(io, "## Estimation")
    println(io, "- Localizations after filter: $(length(result.smld.emitters))")
    if p2_est !== nothing
        println(io, "- **p₂_estimated = $(round(100 * p2_est, digits=3))%**")
        p2_true > 0 && println(io, "- Ratio (est/true): $(round(p2_est / p2_true, digits=2))")
    else
        println(io, "- p₂ estimation failed")
    end
end

println("\nOutput: $outdir")
