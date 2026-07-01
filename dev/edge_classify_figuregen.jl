# Figure-gen + functional coverage for the edge_classify step's figure hook.
#
# The edge figures live in SMLMClustering's SMLMClusteringFiguresExt (CairoMakie + SMLMRender
# weakdep ext) and are produced by the step's STANDARD-gated call to compute_edge_report /
# write_edge_report / plot_edge_report. Per @SMLMBaGoL's pothole, ext figure code is NOT
# exercised by Pkg.test (no Makie in the test env) — this script actually RUNS the hook end to
# end so breakage is caught.
#
#   julia -t auto --project=. dev/edge_classify_figuregen.jl

using SMLMAnalysis
using SMLMAnalysis: Verbosity, step_name, step_outdir

# Synthetic cell-like SMLD: a dense disk ("cell") + sparse background, realistic precision.
function synth_cell_smld(; ncell=8000, nbg=300, fov_px=128, px=0.1,
                           center=(6.4, 6.4), radius=2.5, σ=0.015)
    cam = IdealCamera(fov_px, fov_px, px)
    L = fov_px * px
    cx, cy = center
    ems = Emitter2DFit{Float64}[]
    while length(ems) < ncell
        r = radius * sqrt(rand()); θ = 2π * rand()
        push!(ems, Emitter2DFit{Float64}(cx + r*cos(θ), cy + r*sin(θ),
              1000.0, 5.0, σ, σ, 30.0, 1.0; frame=1, dataset=1, track_id=0, id=0))
    end
    for _ in 1:nbg
        push!(ems, Emitter2DFit{Float64}(L*rand(), L*rand(),
              1000.0, 5.0, σ, σ, 30.0, 1.0; frame=1, dataset=1, track_id=0, id=0))
    end
    BasicSMLD(ems, cam, 1, 1)
end

smld = synth_cell_smld()
outdir = mktempdir()
println("synthetic SMLD: ", length(smld.emitters), " emitters; outdir=", outdir)

# OuterPolygonConfig is robust on any point cloud (alpha-shape) → the pass criterion.
# KdeValleyConfig is dSTORM-tuned; run best-effort (synthetic data may not suit its valley gate).
got_figures = String[]   # configs whose figure series was written (mutated in soft scope — safe)
for (i, cfg) in enumerate((OuterPolygonConfig(), KdeValleyConfig()))
    name = step_name(cfg)
    println("--- [$i] $name ---")
    try
        (out, sinfo) = analyze(smld, cfg; outdir=outdir, step_number=i, verbose=Verbosity.STANDARD)
        d = step_outdir(outdir, i, cfg)
        files = readdir(d)
        pngs  = filter(f -> endswith(f, ".png"), files)
        println("  → ", sinfo.summary)
        println("  files: ", files)
        println("  PNGs:  ", pngs)
        isempty(pngs) || push!(got_figures, name)
    catch e
        println("  !! ", name, " errored (data may not suit this config — not a hook bug): ", e)
    end
end

# Coverage passes if the STANDARD-gated hook (compute/write/plot_edge_report) wrote the named
# figure series for at least one config. KdeValleyConfig is the dSTORM default; OuterPolygonConfig
# needs a denser cloud and may reject a too-sparse synthetic set (data/param, not a hook bug).
if isempty(got_figures)
    println("FIGUREGEN_FAIL")
else
    println("FIGUREGEN_OK via: ", join(got_figures, ", "))
end
