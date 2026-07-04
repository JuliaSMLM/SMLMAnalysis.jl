# Step Configs & Info

```@meta
CurrentModule = SMLMAnalysis
```

Each pipeline step is driven by a config struct and returns a typed info struct.
See [The Pipeline Model](@ref) for how `analyze` dispatches on these types, and
[Pipeline Steps: Overview](@ref) for the per-step reference pages.

## Step Configs

Native to SMLMAnalysis:

```@docs
DetectFitConfig
FilterConfig
IntensityFilterConfig
DensityFilterConfig
```

Re-exported from upstream packages (the owning package documents the algorithm;
SMLMAnalysis dispatches `analyze` on these types — see [The Pipeline Model](@ref)):

```@docs
FrameConnectConfig
CalibrationConfig
DriftConfig
RenderConfig
BaGoLConfig
```

!!! note "CalibrationConfig is a sub-config"
    `CalibrationConfig` configures uncertainty calibration *inside*
    `FrameConnectConfig` (via its `calibration=` field) — it is **not** a separate
    pipeline step. `CalibrationResult` holds its output.

### Clustering (re-exported)

The clustering verbs `cluster` / `cluster_statistics` and their config types
(`DBSCANConfig`, `HDBSCANConfig`, `HierarchicalConfig`, `VoronoiConfig`,
`HopkinsConfig`, `VoronoiDensityConfig`) are re-exported from
[SMLMClustering](https://github.com/JuliaSMLM/SMLMClustering.jl) and dispatched as
pipeline steps — see [Clustering](@ref clustering-step). Their full API and
algorithm reference lives in the SMLMClustering manual; SMLMAnalysis adds only the
`analyze()` dispatch (it does not re-document the upstream API here).

## Step Info Types

Each step's `analyze()` returns a [`StepInfo`](@ref) whose `.info` field holds the
step's typed info struct.

```@docs
DetectFitInfo
FilterInfo
IntensityFilterInfo
DensityFilterInfo
BaGoLInfo
```

!!! note "Upstream info types"
    Some steps return their upstream packages' info structs:
    `SMLMFrameConnection.FrameConnectInfo`, `SMLMDriftCorrection.DriftInfo`,
    `SMLMRender.RenderInfo`, and the clustering `ClusterInfo` /
    `ClusterStatisticsInfo` — documented in those packages. The connected tracks
    and drift model are also surfaced on [`AnalysisResult`](@ref) as
    `smld_connected` and `drift_model`.
