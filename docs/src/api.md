# [API Reference](@id API-Reference)

```@meta
CurrentModule = SMLMAnalysis
```

## Core Functions

```@docs
analyze
reset!
checkpoint!
debug!
get_analysis_info
get_config
```

## Types

```@docs
Analysis
AnalysisConfig
AnalysisInfo
StepRecord
DataSource
```

## Step Configs

```@docs
DetectFitConfig
FilterConfig
FrameConnectConfig
DriftCorrectConfig
DensityFilterConfig
BaGoLConfig
```

!!! note "RenderConfig"
    `RenderConfig` is re-exported from SMLMRender.jl and used directly as a step config.
    See [SMLMRender documentation](https://github.com/JuliaSMLM/SMLMRender.jl) for details.

## I/O

```@docs
save_smld
load_smld
smld_info
resume_analysis
load_smart_h5
load_smart_h5_info
load_lidkelab_h5
load_lidkelab_h5_info
load_lidkelab_h5_block
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
```
