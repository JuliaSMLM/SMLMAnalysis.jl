"""
Clustering step for SMLMAnalysis.

Wraps SMLMClustering's two verbs as pipeline steps:

- `cluster` (label backends: DBSCAN / HDBSCAN / Voronoi / Hierarchical) — assigns
  a cluster id to every localization (`emitter.id`: `0` = noise, `1..K` = cluster).
- `cluster_statistics` (read-only backends: Hopkins / Voronoi density) — computes
  a spatial-statistic summary and returns the SMLD unchanged.

Upstream owns the config types (`DBSCANConfig`, `VoronoiConfig`, `HopkinsConfig`,
…), already re-exported by SMLMAnalysis as `const` aliases. This file adds only
the `analyze()` dispatch that threads those configs through the pipeline, so they
compose like any other step — e.g. `RenderConfig` after `DBSCANConfig` renders the
labeled SMLD. The backend is selected entirely by the concrete config type.
"""

# Cluster-label summary (cluster() backends: DBSCAN / HDBSCAN / Voronoi / Hierarchical)
_step_summary(info::SMLMClustering.ClusterInfo) = Dict{Symbol,Any}(
    :algorithm   => info.algorithm,
    :n_locs_in   => info.n_locs_in,
    :n_clustered => info.n_clustered,
    :n_noise     => info.n_noise,
    :n_clusters  => info.n_clusters,
)

# Spatial-statistic summary (cluster_statistics() backends: Hopkins / Voronoi density)
_step_summary(info::SMLMClustering.ClusterStatisticsInfo) = Dict{Symbol,Any}(
    :algorithm      => info.algorithm,
    :n_locs_in      => info.n_locs_in,
    :statistic_name => info.statistic_name,
    :statistic      => info.statistic,
)

"""
    analyze(smld, cfg::AbstractClusterConfig; kwargs...) -> (smld_out, StepInfo)

Label each localization with a cluster id via SMLMClustering. The backend is
selected by the concrete type of `cfg` (`DBSCANConfig`, `HDBSCANConfig`,
`VoronoiConfig`, `HierarchicalConfig`). The input SMLD is **not modified** — a
deep-copied, labeled SMLD is returned and threaded onward (`emitter.id`: `0` =
noise, `1..K` = cluster). See the SMLMClustering documentation for the algorithm
details and per-backend configuration.
"""
function analyze(smld::BasicSMLD, cfg::SMLMClustering.AbstractClusterConfig;
                 outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD,
                 checkpoint::Int=Checkpoint.EXPENSIVE, kwargs...)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)
    v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg))" n_locs=length(smld.emitters)

    t = @elapsed ((smld_out, info) = SMLMClustering.cluster(smld, cfg))

    if dir !== nothing && v >= Verbosity.STANDARD
        mkpath(dir)
        _save_config!(dir, cfg)
        _save_info!(dir, info)
    end
    if dir !== nothing && checkpoint >= Checkpoint.ALL
        _save_step_smld(dir, smld_out; filename="smld_clustered.jld2")
    end

    v >= Verbosity.PROGRESS && @info "  → $(info.n_clusters) clusters, $(info.n_clustered)/$(info.n_locs_in) clustered ($(info.n_noise) noise)"
    (smld_out, StepInfo(step_number, cfg, t, _step_summary(info); info=info))
end

"""
    analyze(smld, cfg::AbstractStatisticsConfig; kwargs...) -> (smld, StepInfo)

Compute a read-only spatial statistic via SMLMClustering — e.g. Hopkins
clustering tendency (`HopkinsConfig`) or Voronoi density
(`VoronoiDensityConfig`). The SMLD is returned **unchanged**; the result scalar
and any per-emitter/per-dataset vectors live in the step's `ClusterStatisticsInfo`
(`info.statistic`, `info.extras`).
"""
function analyze(smld::BasicSMLD, cfg::SMLMClustering.AbstractStatisticsConfig;
                 outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD,
                 checkpoint::Int=Checkpoint.EXPENSIVE, kwargs...)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)
    v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg))" n_locs=length(smld.emitters)

    t = @elapsed ((smld_out, info) = SMLMClustering.cluster_statistics(smld, cfg))

    if dir !== nothing && v >= Verbosity.STANDARD
        mkpath(dir)
        _save_config!(dir, cfg)
        _save_info!(dir, info)
    end

    v >= Verbosity.PROGRESS && @info "  → $(info.statistic_name) = $(round(info.statistic, digits=4))"
    (smld_out, StepInfo(step_number, cfg, t, _step_summary(info); info=info))
end
