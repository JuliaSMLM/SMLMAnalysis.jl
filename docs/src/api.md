# [API Reference](@id API-Reference)

```@meta
CurrentModule = SMLMAnalysis
```

```@docs
SMLMAnalysis
```

## Core Functions

```@docs
analyze
```

## Types

```@docs
AnalysisConfig
AnalysisResult
AnalysisInfo
StepInfo
DataSource
Checkpoint
```

### Verbosity Levels

`Verbosity` is a module with integer constants controlling output detail:

| Level | Constant | Output |
|-------|----------|--------|
| 0 | `Verbosity.SILENT` | Errors only |
| 1 | `Verbosity.PROGRESS` | Step names, counts, timing |
| 2 | `Verbosity.STANDARD` | + stats.md, basic figures |
| 3 | `Verbosity.DETAILED` | + diagnostic plots, per-filter breakdowns |
| 4 | `Verbosity.DEBUG` | + MP4 animations, frame-by-frame analysis |

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
CompositeRenderInfo
CrossAlignInfo
CrossCorrInfo
```

!!! note "Upstream info types"
    Some steps return their upstream packages' info structs:
    `SMLMFrameConnection.FrameConnectInfo`, `SMLMDriftCorrection.DriftInfo`,
    `SMLMRender.RenderInfo`, and the clustering `ClusterInfo` /
    `ClusterStatisticsInfo` — documented in those packages. The connected tracks
    and drift model are also surfaced on [`AnalysisResult`](@ref) as
    `smld_connected` and `drift_model`.

## Multi-Target

```@docs
MultiTargetConfig
MultiTargetResult
MultiTargetInfo
AbstractMultiTargetStep
CompositeRenderConfig
CrossAlignConfig
CrossCorrConfig
```

## I/O

```@docs
save_smld
load_smld
smld_info
save_pipeline_state
load_pipeline_state
load_smart_h5
load_smart_h5_info
load_smart_h5_frame
smart_h5_to_array
load_mic_h5
load_mic_h5_info
load_mic_h5_block
load_mic_h5_calibration
load_mic_h5_calibration_for_scmos
build_camera_from_mic_h5
```

## Utilities

```@docs
crop_camera
crop_images
step_name
step_outdir
n_datasets
n_frames_per_dataset
```

## Internals

```@autodocs
Modules = [SMLMAnalysis]
Public = false
```
