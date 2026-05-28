```@meta
CurrentModule = SMLMAnalysis
```

# SMLMAnalysis.jl

High-level integration package for the [JuliaSMLM](https://github.com/JuliaSMLM) ecosystem. Orchestrates detection, fitting, filtering, frame connection, drift correction, density filtering, rendering, and Bayesian grouping into reproducible SMLM analysis pipelines.

## Two Workflow Patterns

### Config-driven (production)

```julia
# Steps are composable — reorder, repeat, or omit freely after DetectFitConfig
config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(
            boxer=BoxerConfig(boxsize=9, psf_sigma=0.130),
            fitter=GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20)),
        FilterConfig(photons=(500.0, Inf)),
        DriftConfig(degree=2),
        RenderConfig(zoom=20, colormap=:inferno),
    ],
    outdir = "output/",
)
(result, info) = analyze(image_stacks, config)
```

### Step-by-step (exploration)

```julia
(smld, df_info) = analyze(image_stacks, DetectFitConfig(
    camera=camera, boxer=BoxerConfig(boxsize=9, psf_sigma=0.130)))

(smld, _) = analyze(smld, FilterConfig(photons=(500.0, Inf)))
(smld, _) = analyze(smld, FrameConnectConfig(max_frame_gap=5,
    calibration=CalibrationConfig(clamp_k_to_one=true)))
(smld, _) = analyze(smld, DriftConfig(degree=2))
(img, _)  = analyze(smld, RenderConfig(zoom=20, colormap=:inferno))

# Save intermediate state for later resume
save_smld("after_detectfit.h5", smld)
```

## Documentation

- **[Tutorial](@ref)** -- Full pipeline walkthrough with figures at every step
- **[Guide](@ref)** -- Concepts: pipeline architecture, extending with new steps, multi-dataset, drift modes, I/O
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
    +-- SMLMFrameConnection (linking + uncertainty calibration)
    +-- SMLMDriftCorrection (entropy-based drift correction)
    +-- SMLMRender (super-resolution rendering)
    +-- SMLMSim (simulation + image generation)
    +-- MicroscopePSFs (PSF models)
    +-- SMLMBaGoL (Bayesian grouping)
    +-- SMLMClustering (DBSCAN / Hierarchical / Voronoi)
    |
    +-- SMLMAnalysis (integrates all)
```

```@index
```
