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
StepRecord
DataSource
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

```@docs
DetectFitConfig
FilterConfig
FrameConnectConfig
DriftCorrectConfig
DensityFilterConfig
```

!!! note "RenderConfig"
    `RenderConfig` is re-exported from SMLMRender.jl and used directly as a step config.
    See [SMLMRender documentation](https://github.com/JuliaSMLM/SMLMRender.jl) for details.

## Multi-Target

```@docs
MultiTargetConfig
MultiTargetResult
MultiTargetInfo
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
```

## Calibration

```@docs
analyze_frameconnect_drift
apply_uncertainty_calibration
recombine_tracks
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
