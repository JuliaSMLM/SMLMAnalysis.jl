#
# dev/deepfit_threshold_sweep.jl — map DECODE detection precision/recall vs mintotalp.
#
# Reuses one trained model + one simulated movie; re-runs only deepfit() inference per
# threshold (network forward is cheap after the first Reactant compile). Produces the
# precision-recall tradeoff curve so we pick mintotalp on evidence, not a guess.
#
#   DEEPFIT_MODEL=<path/model.jld2> julia -t auto --project=. dev/deepfit_threshold_sweep.jl
#
include(joinpath(@__DIR__, "deepfit_sim.jl"))   # helpers only; main() is guarded, not run

function sweep()
    model_path = get(ENV, "DEEPFIT_MODEL", "")
    isempty(model_path) && error("set DEEPFIT_MODEL=<model.jld2>")
    p = make_decode(PSF_PATH)
    cam, movie, gt = simulate_structure(PSF_PATH)
    movie_f = Float32.(movie); movie_f ./= maximum(movie_f)
    @info "sweep setup" model_path n_gt = length(gt.emitters) movie = size(movie)

    thresholds = Float32[0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
    rows = NamedTuple[]
    for mt in thresholds
        cfg = SMLMDeepFit.DeepFitConfig(; model_path = model_path, traintype = p, camera = cam,
                                        batchsize = 8, mintotalp = mt)
        (smld, _) = analyze(movie_f, cfg; outdir = OUTDIR, step_number = 4,
                            verbose = SMLMAnalysis.Verbosity.SILENT)
        m1 = metrics_at(smld, gt, 0.1); m2 = metrics_at(smld, gt, 0.2)
        push!(rows, (mt = mt, n = length(smld.emitters),
                     r100 = m1.recall, p100 = m1.precision, rmse100 = m1.rmse_xy,
                     r200 = m2.recall, p200 = m2.precision, rmse200 = m2.rmse_xy))
        @printf "mintotalp=%.2f  n=%3d  | 100nm r=%.2f p=%.2f rmse=%.0f | 200nm r=%.2f p=%.2f rmse=%.0f\n" mt length(smld.emitters) m1.recall m1.precision 1000*m1.rmse_xy m2.recall m2.precision 1000*m2.rmse_xy
    end

    dir = joinpath(OUTDIR, "04_deepfit_inference"); mkpath(dir)
    open(joinpath(dir, "threshold_sweep.md"), "w") do io
        println(io, "# deepfit_inference — detection threshold sweep\n")
        println(io, "GT emitters: $(length(gt.emitters)) | model: `$(basename(dirname(model_path)))`\n")
        println(io, "| mintotalp | n_inferred | recall@100 | prec@100 | RMSE@100 | recall@200 | prec@200 | RMSE@200 |")
        println(io, "|---|---|---|---|---|---|---|---|")
        for r in rows
            @printf io "| %.2f | %d | %.2f | %.2f | %.0f | %.2f | %.2f | %.0f |\n" r.mt r.n r.r100 r.p100 1000*r.rmse100 r.r200 r.p200 1000*r.rmse200
        end
    end

    mts = [r.mt for r in rows]
    fig = Figure(size = (760, 340))
    ax = Axis(fig[1, 1], xlabel = "mintotalp", ylabel = "recall / precision",
              title = "DECODE detection tradeoff (200nm match)")
    lines!(ax, mts, [r.r200 for r in rows], color = :seagreen, label = "recall")
    scatter!(ax, mts, [r.r200 for r in rows], color = :seagreen)
    lines!(ax, mts, [r.p200 for r in rows], color = :crimson, label = "precision")
    scatter!(ax, mts, [r.p200 for r in rows], color = :crimson)
    axislegend(ax, position = :rc)
    ax2 = Axis(fig[1, 2], xlabel = "mintotalp", ylabel = "n inferred",
               title = "detections vs threshold")
    lines!(ax2, mts, [Float64(r.n) for r in rows], color = :steelblue)
    scatter!(ax2, mts, [Float64(r.n) for r in rows], color = :steelblue)
    hlines!(ax2, [Float64(length(gt.emitters))], color = (:black, 0.5), linestyle = :dash)
    CairoMakie.save(joinpath(dir, "threshold_sweep.png"), fig)
    @info "threshold sweep written" dir
    rows
end

sweep()
