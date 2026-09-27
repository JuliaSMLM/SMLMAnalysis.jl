# Multi-Target & I/O

```@meta
CurrentModule = SMLMAnalysis
```

## Multi-Target

Multi-channel (multi-color) analysis types. The orchestrator config
[`MultiTargetConfig`](@ref) is documented on the [Multi-Channel Analysis](@ref)
page; its result, info, and step types are below. See [Composite Render](@ref) /
[Cross-Alignment](@ref) / [Cross-Correlation](@ref) for the individual steps.

```@docs
SMLMAnalysis.MultiTargetResult
SMLMAnalysis.MultiTargetInfo
SMLMAnalysis.AbstractMultiTargetStep
CompositeRenderConfig
CrossAlignConfig
CrossCorrConfig
SMLMAnalysis.CompositeRenderInfo
SMLMAnalysis.CrossAlignInfo
SMLMAnalysis.CrossCorrInfo
```

`MultiTargetResult`/`MultiTargetInfo`/`AbstractMultiTargetStep` and the
`*Info` result structs are not exported — reach them as `SMLMAnalysis.Name`.

## I/O

Save and load localization results, and import microscope `.h5` data. See
[I/O & Resume](@ref) for the workflow.

```@docs
save_smld
load_smld
smld_info
load_smart_h5
load_smart_h5_info
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
SMLMAnalysis.step_name
SMLMAnalysis.step_outdir
```
