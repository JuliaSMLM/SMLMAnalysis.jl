# SMLMAnalysis.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSMLM.github.io/SMLMAnalysis.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/SMLMAnalysis.jl/dev/)
[![Build Status](https://github.com/JuliaSMLM/SMLMAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSMLM/SMLMAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSMLM/SMLMAnalysis.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSMLM/SMLMAnalysis.jl)

High-level integration package for the [JuliaSMLM](https://github.com/JuliaSMLM) ecosystem. Orchestrates detection, fitting, filtering, frame connection, drift correction, density filtering, rendering, and Bayesian grouping into reproducible SMLM analysis pipelines.

## Philosophy

- **Functional pipeline** -- pure step functions returning `(result, info)` tuples
- **Typed configs** for every step -- reproducible, serializable, shareable
- **Dataset boundaries from data structure** -- `Vector{Array}` encodes multi-dataset boundaries
- **All coordinates in microns** -- consistent across the ecosystem

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

## Installation

```julia
using Pkg
Pkg.add("SMLMAnalysis")
```

## Quick Start

### One-liner with AnalysisConfig

```julia
using SMLMAnalysis

config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(boxsize=9, psf_model=:variable),
        FilterConfig(photons=(500.0, Inf), precision=(0.0, 0.007)),
        FrameConnectConfig(max_frame_gap=5),
        DriftCorrectConfig(degree=2),
        DensityFilterConfig(n_sigma=2.0),
        RenderConfig(zoom=20, colormap=:inferno),
    ],
    outdir = "output/",
)

(result, info) = analyze(image_stacks, config)
result.smld                # Final SMLD
info.steps[:detectfit]     # Per-step info from upstream packages
```

### Step-by-step with pure functions

```julia
(smld, df_info) = detectfit(image_stacks, camera, DetectFitConfig(boxsize=9, psf_model=:variable))
smld_raw = df_info.smld_raw

(smld, _) = filter_step(smld, FilterConfig(photons=(500.0, Inf)); smld_raw=smld_raw)
(smld, _) = frameconnect_step(smld, FrameConnectConfig(max_frame_gap=5))
(smld, _) = driftcorrect_step(smld, DriftCorrectConfig(degree=2))
(img, _)  = render_step(smld, RenderConfig(zoom=20, colormap=:inferno))

# Save intermediate state for resume
save_smld("after_detectfit.h5", smld)
smld = load_smld("after_detectfit.h5")
```

### Rendered output

![Super-resolution render](docs/src/assets/render_gaussian.png)

## Pipeline Steps

| Step | Config | Description |
|------|--------|-------------|
| Detection + Fitting | `DetectFitConfig` | ROI detection and GPU-accelerated MLE fitting |
| Filtering | `FilterConfig` | Filter by photons, precision, p-value, PSF width |
| Frame Connection | `FrameConnectConfig` | Link localizations across frames, uncertainty calibration |
| Drift Correction | `DriftCorrectConfig` | Entropy-based drift correction (continuous or registered) |
| Density Filter | `DensityFilterConfig` | Remove isolated localizations by neighbor count |
| Render | `RenderConfig` | Gaussian, histogram, circle, or ellipse rendering |
| BaGoL | `BaGoLConfig` | Bayesian grouping of localizations |

## Related Packages

| Package | Description |
|---------|-------------|
| [SMLMData.jl](https://github.com/JuliaSMLM/SMLMData.jl) | Core types: Emitter, Camera, BasicSMLD |
| [SMLMBoxer.jl](https://github.com/JuliaSMLM/SMLMBoxer.jl) | ROI detection from raw images |
| [GaussMLE.jl](https://github.com/JuliaSMLM/GaussMLE.jl) | GPU-accelerated MLE fitting |
| [SMLMFrameConnection.jl](https://github.com/JuliaSMLM/SMLMFrameConnection.jl) | Linking localizations across frames |
| [SMLMDriftCorrection.jl](https://github.com/JuliaSMLM/SMLMDriftCorrection.jl) | Entropy-based drift correction |
| [SMLMRender.jl](https://github.com/JuliaSMLM/SMLMRender.jl) | Super-resolution image rendering |
| [SMLMSim.jl](https://github.com/JuliaSMLM/SMLMSim.jl) | SMLM data simulation |
| [MicroscopePSFs.jl](https://github.com/JuliaSMLM/MicroscopePSFs.jl) | PSF models |
| [SMLMBaGoL.jl](https://github.com/JuliaSMLM/SMLMBaGoL.jl) | Bayesian grouping of localizations |

## License

MIT License. See [LICENSE](LICENSE) for details.
