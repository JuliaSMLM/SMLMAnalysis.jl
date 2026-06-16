#
# dev/deepfit_sim.jl — stages 2-5 of the deepfit sim chain (train -> sim -> infer -> render + compare).
# Run in the dev/ env on the GPU; consumes psf.h5 from dev/learn_psf_sim.jl.
#
#   DEEPFIT_EPOCHS=40 julia -t auto --project=. dev/deepfit_sim.jl        # full run
#   DEEPFIT_MODEL=<path/model.jld2> julia -t auto --project=. dev/deepfit_sim.jl   # skip training, reuse a model
#
using SMLMData
using SMLMAnalysis
using SMLMDeepFit
using SMLMSim
using MicroscopePSFs
using SMLMRender
using CairoMakie
using NearestNeighbors
using Statistics
using Printf
using Random
import Optimisers

include(joinpath(@__DIR__, "deepfit_sim_common.jl"))

const TRAIN_EPOCHS = parse(Int, get(ENV, "DEEPFIT_EPOCHS", "40"))

Random.seed!(1)
mkpath(OUTDIR)

# ============================================================
# Stage 2 — deepfit_training
# ============================================================
function make_decode(psf_path)
    p = SMLMDeepFit.Decode(; sz = 40, ρ = 1.0f0, photons = 500.0f0, bg = 5.0f0,
                           n_train = 500, n_test = 100, psffile = psf_path,
                           pixelsize = Float32(PIXEL_SIZE))   # explicit camera pixel size (no silent 0.1 guess)
    p.netconfig = SMLMDeepFit.UNetConfig(n_features = 16, depth = 2)   # depth=2 ⇒ FOV÷4
    p
end

function train_decode(p)
    args = SMLMDeepFit.TrainConfig(; traintype = p, use_reactant = true, epochs = TRAIN_EPOCHS,
                                   batchsize = 10, optimiser = Optimisers.Adam(5.0f-4),  # 1e-3 diverged late
                                   infotime = 5, checktime = 0, tblogger = false,
                                   savepath = joinpath(OUTDIR, "02_deepfit_training"))
    (result, sinfo) = analyze(args; outdir = OUTDIR, step_number = 2,
                              verbose = SMLMAnalysis.Verbosity.PROGRESS)
    plot_training_curves(sinfo.info)
    @info "stage 2 (deepfit_training) done" model_path = result.model_path epochs = TRAIN_EPOCHS
    result.model_path
end

# Loss/accuracy curves from TrainInfo (SMLMDeepFit's own LossPlot/AccuracyPlot come out empty).
function plot_training_curves(ti)
    if ti === nothing || isempty(ti.train_losses)
        @warn "no loss history in TrainInfo (train_losses empty) — flag to @smlmdeepfit"
        return
    end
    try
        dir = joinpath(OUTDIR, "02_deepfit_training"); ep = 1:length(ti.train_losses)
        fig = Figure(size = (940, 380))
        ax1 = Axis(fig[1, 1], xlabel = "logged step", ylabel = "loss", title = "DECODE training loss")
        lines!(ax1, ep, ti.train_losses, label = "train")
        length(ti.test_losses) == length(ep) && lines!(ax1, ep, ti.test_losses, label = "test")
        axislegend(ax1)
        ax2 = Axis(fig[1, 2], xlabel = "logged step", ylabel = "accuracy / efficiency", title = "DECODE accuracy")
        lines!(ax2, ep, ti.train_accuracies, label = "train")
        length(ti.test_accuracies) == length(ep) && lines!(ax2, ep, ti.test_accuracies, label = "test")
        axislegend(ax2)
        CairoMakie.save(joinpath(dir, "loss_accuracy.png"), fig)
        @info "training curves" final_train_loss = ti.train_losses[end] best_epoch = ti.best_epoch n = length(ep)
    catch e
        @warn "training curve plot failed" exception = e
    end
end

# ============================================================
# Stage 3 — known 3D structure -> movie + ground truth
# ============================================================
function simulate_structure(psf_path; n_frames = 200)
    psf = MicroscopePSFs.load_psf(psf_path)
    cam = IdealCamera(FOV, FOV, PIXEL_SIZE)
    sim = StaticSMLMConfig(density = 1.0, σ_psf = 0.13, nframes = n_frames, ndatasets = 1,
                           ndims = 3, zrange = [-0.5, 0.5])
    (_, si) = simulate(sim;
        pattern  = Nmer3D(n = 8, d = 0.05),
        molecule = GenericFluor(photons = 5.0e3, k_off = 20.0, k_on = 0.04),
        camera   = cam)
    (movie, _) = gen_images(si.smld_model, psf; dataset = 1, bg = 10.0, poisson_noise = true)
    cam, movie, si.smld_model
end

# ============================================================
# Stage 4 — deepfit_inference
# ============================================================
function run_inference(movie, p, model_path, cam)
    movie_f = Float32.(movie)
    movie_f ./= maximum(movie_f)                 # match training scale ≈[0,1]; inference does NOT normalize
    cfg = SMLMDeepFit.DeepFitConfig(; model_path = model_path, traintype = p, camera = cam, batchsize = 8,
                                    mintotalp = parse(Float32, get(ENV, "DEEPFIT_MINTOTALP", "0.8")))  # sweep optimum: prec=1.0 @ recall plateau
    (smld, _) = analyze(movie_f, cfg; outdir = OUTDIR, step_number = 4,
                        verbose = SMLMAnalysis.Verbosity.PROGRESS)
    smld
end

# ============================================================
# Stage 5 — compare to ground truth + render
# ============================================================
xy(s) = (Float64[e.x for e in s.emitters], Float64[e.y for e in s.emitters])

"""Lateral NN match within `radius` µm (pooled over frames). Recall = GT found; precision = inferred TP."""
function metrics_at(inferred, gt, radius)
    isempty(inferred.emitters) && return (radius = radius, recall = 0.0, precision = 0.0, rmse_xy = NaN, n_tp = 0)
    gx, gy = xy(gt); ix, iy = xy(inferred)
    G = permutedims(hcat(gx, gy)); I = permutedims(hcat(ix, iy))     # 2×N point clouds
    _, di = knn(KDTree(G), I, 1); tp = [d[1] <= radius for d in di]
    sq = [d[1]^2 for (m, d) in zip(tp, di) if m]
    _, dg = knn(KDTree(I), G, 1); n_found = count(d[1] <= radius for d in dg)
    (radius = radius, recall = n_found / length(gx), precision = count(tp) / length(ix),
     rmse_xy = isempty(sq) ? NaN : sqrt(mean(sq)), n_tp = count(tp))
end

"""Direct alignment check: GT vs inferred positions on one axis."""
function scatter_overlay(inferred, gt, dir)
    gx, gy = xy(gt); ix, iy = xy(inferred)
    fig = Figure(size = (560, 560))
    ax = Axis(fig[1, 1], xlabel = "x (µm)", ylabel = "y (µm)",
              title = "positions: GT (gray) vs inferred (red)", aspect = DataAspect())
    scatter!(ax, gx, gy, color = (:gray, 0.6), markersize = 10, label = "GT")
    scatter!(ax, ix, iy, color = (:red, 0.7), markersize = 5, label = "inferred")
    axislegend(ax)
    CairoMakie.save(joinpath(dir, "positions_overlay.png"), fig)
end

function render_compare(inferred, gt, dir)
    # GT-vs-inferred comparison = positions_overlay.png (scatter). GT model emitters have σ=0,
    # which blanks every σ-based SR render (Gaussian/Circle), so render only the inferred SR here.
    if isempty(inferred.emitters)
        @warn "inferred SMLD empty — skipping render"
        return
    end
    (img_i, _) = render(inferred; strategy = GaussianRender(), zoom = 20)
    save_image(joinpath(dir, "inferred_sr.png"), img_i)
end

"""deepfit_inference detectfit-style overlay: inferred localizations (red) on sample movie frames."""
function inference_overlay(movie, inferred, dir; n_show = 6)
    isempty(inferred.emitters) && return
    H, W, T = size(movie)
    fr_with = sort(unique(Int[e.frame for e in inferred.emitters if 1 <= e.frame <= T]))
    isempty(fr_with) && return
    sel = fr_with[round.(Int, range(1, length(fr_with), length = min(n_show, length(fr_with))))]
    nc = 3; nr = cld(length(sel), nc)
    fig = Figure(size = (nc * 250, nr * 250 + 30))
    Label(fig[0, 1:nc], "deepfit_inference: localizations (red ×) on movie frames", fontsize = 12)
    for (i, fr) in enumerate(sel)
        r = div(i - 1, nc) + 1; c = mod(i - 1, nc) + 1
        ax = Axis(fig[r, c], title = "frame $fr", aspect = DataAspect(), yreversed = true)
        heatmap!(ax, movie[:, :, fr]', colormap = :grays)
        ex = Float64[e.x / PIXEL_SIZE for e in inferred.emitters if e.frame == fr]
        ey = Float64[e.y / PIXEL_SIZE for e in inferred.emitters if e.frame == fr]
        scatter!(ax, ex, ey, color = :red, markersize = 9, marker = :xcross)
        hidedecorations!(ax)
    end
    CairoMakie.save(joinpath(dir, "inference_overlay.png"), fig)
end

"""Validation summary: reconstruction metrics at several match radii."""
function write_validation(ms, inferred, gt, dir)
    open(joinpath(dir, "validation.md"), "w") do io
        println(io, "# deepfit_inference validation (vs simulated ground truth)\n")
        println(io, "- inferred localizations: $(length(inferred.emitters))")
        println(io, "- ground-truth emitters:  $(length(gt.emitters))\n")
        println(io, "| match radius | recall | precision | RMSE_xy (nm) | TP |")
        println(io, "|---|---|---|---|---|")
        for m in ms
            println(io, "| $(round(Int, 1000*m.radius)) nm | $(round(m.recall, digits=3)) | $(round(m.precision, digits=3)) | $(round(1000*m.rmse_xy, digits=1)) | $(m.n_tp) |")
        end
    end
end

function main()
    isfile(PSF_PATH) || error("psf.h5 missing — run `julia --project=. dev/learn_psf_sim.jl` first")
    p = make_decode(PSF_PATH)
    model_path = get(ENV, "DEEPFIT_MODEL", "")
    if isempty(model_path)
        @info "=== stage 2 deepfit_training ==="; model_path = train_decode(p)
    else
        @info "=== skip training — reusing model ===" model_path
    end
    @info "=== stage 3 structure sim ==="; cam, movie, gt = simulate_structure(PSF_PATH)
    @info "movie" size = size(movie) n_gt = length(gt.emitters)
    @info "=== stage 4 deepfit_inference ==="; smld = run_inference(movie, p, model_path, cam)
    @info "=== stage 5 compare + render ==="
    dir = joinpath(OUTDIR, "05_compare"); mkpath(dir)
    infdir = joinpath(OUTDIR, "04_deepfit_inference"); mkpath(infdir)
    # inference_overlay.png + inferred_sr.png are now produced by SMLMAnalysisDeepFitExt during
    # analyze(movie, DeepFitConfig) — the step owns its production diagnostics (like detectfit).
    # The harness keeps only the GT-based validation, which can't live in a step (no ground truth
    # at inference time): the GT-vs-inferred scatter, the multi-radius metrics, and validation.md.
    scatter_overlay(smld, gt, dir)
    ms = [metrics_at(smld, gt, r) for r in (0.05, 0.1, 0.2)]
    write_validation(ms, smld, gt, infdir)
    println("\nRESULT: ", length(smld.emitters), " inferred vs ", length(gt.emitters), " GT")
    for m in ms
        @printf "  r=%3dnm: recall=%.2f precision=%.2f rmse_xy=%.1fnm (%d TP)\n" round(Int, 1000 * m.radius) m.recall m.precision 1000 * m.rmse_xy m.n_tp
    end
    println("  → ", dir)
    (smld, gt, ms)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
