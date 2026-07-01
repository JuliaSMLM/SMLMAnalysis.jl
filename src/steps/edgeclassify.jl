"""
Edge-classification step for SMLMAnalysis.

Wraps SMLMClustering's `classify_emitters` as a pipeline step. Each localization
is labeled `:interior` / `:membrane` / `:outside` relative to one or more cell
masks carved from the localization cloud (multi-cell masks + multi-scale adaptive
alpha). The classifier is **state-modifying but emitter-preserving**: the returned
SMLD wraps the *same* emitters, with only the cell-mask **geometry** mirrored into
`metadata`:

- `metadata["edge_cells"]`          — the `MultiCellMask` (`Vector{CellPolygon}`, largest-first)
- `metadata["edge_outer_polygon"]`  — the dominant cell's outer ring (back-compat / Hopkins)

The per-emitter class is **not** mirrored into `metadata` — a per-emitter side-list
desyncs the moment a downstream step subsets emitters. It lives only in the step's
`EdgeClassifyInfo` (`info.class`, with the `in_cell` / `interior_mask` accessors);
consume it at the classify point. Only the geometry travels downstream — e.g.
`HopkinsConfig` reads `edge_cells` / `edge_outer_polygon` as its observation region.

Upstream owns the config types (`OuterPolygonConfig`, `KdeValleyConfig`, both
`<: AbstractEdgeClassifyConfig`), re-exported by SMLMAnalysis as `const` aliases.
This file adds only the `analyze()` dispatch that threads them through the
pipeline, so edge classification composes like any other step. The method is
selected entirely by the concrete config type. See the SMLMClustering
documentation for the algorithm details and per-config parameters.
"""

# Output folder/log label = the STEP's role (`edge_classify`), NOT the upstream
# gate-method name (`kde_valley` / `outer_polygon`), which is internal-mechanism
# jargon at the pipeline level — it doesn't read as "edge classification" at a
# glance. The gate method still rides in `config.toml` + the StepInfo summary
# (`:method`) for provenance, so nothing is lost.
step_name(cfg::SMLMClustering.AbstractEdgeClassifyConfig) = "edge_classify"

_step_summary(info::SMLMClustering.EdgeClassifyInfo) = Dict{Symbol,Any}(
    :method     => SMLMClustering.method_name(info.config),
    :n_emitters => info.n_emitters,
    :n_interior => info.n_interior,
    :n_membrane => info.n_membrane,
    :n_outside  => info.n_outside,
    :n_cells    => length(info.cells),
)

"""
    analyze(smld, cfg::AbstractEdgeClassifyConfig; kwargs...) -> (smld_out, StepInfo)

Classify each localization as `:interior` / `:membrane` / `:outside` against the
cell mask(s) carved by SMLMClustering. The method is selected by the concrete type
of `cfg` (`KdeValleyConfig` — the dSTORM-tuned default — or `OuterPolygonConfig`).
The emitters are **not modified** and are shared into the returned SMLD; only the
cell-mask geometry is mirrored into `metadata` (`edge_cells`, `edge_outer_polygon`)
for downstream steps. The authoritative per-emitter class and the full geometry
live in the step's `EdgeClassifyInfo` (`info.class` / `interior_mask(info)`).
"""
function analyze(smld::BasicSMLD, cfg::SMLMClustering.AbstractEdgeClassifyConfig;
                 outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD,
                 checkpoint::Int=Checkpoint.EXPENSIVE, kwargs...)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)
    v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg))" n_locs=length(smld.emitters)

    t = @elapsed ((smld_out, info) = SMLMClustering.classify_emitters(smld, cfg))

    if dir !== nothing && v >= Verbosity.STANDARD
        mkpath(dir)
        _save_config!(dir, cfg)
        _save_info!(dir, info)
        # Edge report: diagnostics (core) + named figure series (SMLMClusteringFiguresExt,
        # active since SMLMAnalysis hard-deps CairoMakie + SMLMRender). Owner-owned figures,
        # called here — same pattern as bagol_step → compute/write/plot_report.
        report = SMLMClustering.compute_edge_report(smld_out, info)
        SMLMClustering.write_edge_report(report; output_dir=dir)
        SMLMClustering.plot_edge_report(report; output_dir=dir, prefix="edge")
    end
    if dir !== nothing && checkpoint >= Checkpoint.ALL
        _save_step_smld(dir, smld_out; filename="smld_edgeclassified.jld2")
    end

    v >= Verbosity.PROGRESS && @info "  → $(info.n_interior) interior / $(info.n_membrane) membrane / $(info.n_outside) outside, $(length(info.cells)) cell(s)"
    (smld_out, StepInfo(step_number, cfg, t, _step_summary(info); info=info))
end
