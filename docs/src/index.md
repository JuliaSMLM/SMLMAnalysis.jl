```@meta
CurrentModule = SMLMAnalysis
```

# SMLMAnalysis.jl

High-level integration package for the [JuliaSMLM](https://github.com/JuliaSMLM) ecosystem. Orchestrates detection, fitting, filtering, frame connection, drift correction, density filtering, rendering, and Bayesian grouping into reproducible SMLM analysis pipelines.

## Two Workflow Patterns

### Config-driven (production)

```julia
config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(boxsize=9, psf_model=:variable),
        FilterConfig(photons=(500.0, Inf)),
        DriftCorrectConfig(degree=2),
        RenderConfig(zoom=20, colormap=:inferno),
    ],
    outdir = "output/",
)
(result, info) = analyze(image_stacks, config)
```

### Step-by-step (exploration)

```julia
(smld, df_info) = analyze(image_stacks, DetectFitConfig(camera=camera, boxsize=9))
smld_raw = df_info.smld_raw

(smld, _) = analyze(smld, FilterConfig(photons=(500.0, Inf)); smld_raw=smld_raw)
(smld, _) = analyze(smld, FrameConnectConfig(max_frame_gap=5))
(smld, _) = analyze(smld, DriftCorrectConfig(degree=2))
(img, _)  = analyze(smld, RenderConfig(zoom=20, colormap=:inferno))

# Save intermediate state for later resume
save_smld("after_detectfit.h5", smld)
```

## Documentation

- **[Tutorial](@ref)** -- Full pipeline walkthrough with figures at every step
- **[Guide](@ref)** -- Concepts: multi-dataset architecture, drift modes, I/O
- **[API Reference](@ref)** -- Complete function and type documentation

## Installation

```julia
using Pkg
Pkg.add("SMLMAnalysis")
```

## Ecosystem

```
SMLMData (core types)
    |
    +-- SMLMBoxer (ROI detection)
    +-- GaussMLE (GPU-accelerated MLE fitting)
    +-- SMLMFrameConnection (linking across frames)
    +-- SMLMDriftCorrection (entropy-based drift correction)
    +-- SMLMRender (super-resolution rendering)
    +-- SMLMSim (simulation + image generation)
    +-- MicroscopePSFs (PSF models)
    +-- SMLMBaGoL (Bayesian grouping)
    |
    +-- SMLMAnalysis (integrates all)
```

```@index
```
