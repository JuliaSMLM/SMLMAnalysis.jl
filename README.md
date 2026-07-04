# SMLMAnalysis

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSMLM.github.io/SMLMAnalysis.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/SMLMAnalysis.jl/dev/)
[![Build Status](https://github.com/JuliaSMLM/SMLMAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSMLM/SMLMAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSMLM/SMLMAnalysis.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSMLM/SMLMAnalysis.jl)

SMLM analysis pipeline for the [JuliaSMLM](https://github.com/JuliaSMLM) ecosystem: detection, fitting, filtering, frame connection, drift correction, and super-resolution rendering with provenance tracking.

## Installation

```julia
using Pkg
Pkg.add("SMLMAnalysis")
```

## Quick Start

```julia
using SMLMAnalysis

cam = IdealCamera(512, 512, 0.1)  # 512x512, 100nm pixels

config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(boxer=BoxerConfig(boxsize=9, psf_sigma=0.130),
            fitter=GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20)),
        FilterConfig(photons=(500.0, Inf), precision=(0.0, 0.007)),
        FrameConnectConfig(max_frame_gap=5),
        DriftConfig(degree=2, dataset_mode=:registered),
        RenderConfig(zoom=20, colormap=:inferno),
    ],
    outdir = "output/",
)
(result, info) = analyze(image_stacks, config)
```

![Super-resolution render](docs/src/assets/render_gaussian.png)

## Pipeline Architecture

The pipeline folds `analyze()` over a steps vector, with Julia's method dispatch routing each call by `(state_type, config_type)`:

```
for (i, cfg) in enumerate(steps)
    (state, step_info) = analyze(state, cfg)
end
```

Each step defines `analyze(input, ::MyConfig)` and returns `(result, StepInfo)`. Adding a step means defining a config type, implementing `analyze()` for it, and exporting -- it works immediately in both config-driven and step-by-step workflows.

Steps from upstream packages (`DriftConfig`, `FrameConnectConfig`, `RenderConfig`) are re-exported and dispatch directly to their source implementations.

## Step-by-Step Usage

Each step works standalone with the same `analyze()` dispatch:

```julia
cam = IdealCamera(512, 512, 0.1)

(smld, _) = analyze(image_stacks, DetectFitConfig(
    camera=cam, boxer=BoxerConfig(boxsize=9, psf_sigma=0.130)))
(smld, _) = analyze(smld, FilterConfig(photons=(500.0, Inf)))
(smld, _) = analyze(smld, FrameConnectConfig(max_frame_gap=5))
(smld, _) = analyze(smld, DriftConfig(degree=2))
(img, _)  = analyze(smld, RenderConfig(zoom=20, colormap=:inferno))

# Save/load intermediate results
save_smld("checkpoint.h5", smld)
smld = load_smld("checkpoint.h5")
```

## Composability

After `DetectFitConfig` (which produces localizations from images), steps can be used in any combination, order, or repetition:

```julia
# Minimal: detect and render
steps = [DetectFitConfig(boxer=BoxerConfig(boxsize=9)), RenderConfig(zoom=20)]

# Multiple renders at different scales
steps = [DetectFitConfig(...), DriftConfig(degree=2),
         RenderConfig(zoom=10, colormap=:viridis),
         RenderConfig(zoom=20, colormap=:inferno)]

# Repeated filtering: coarse before connection, tight after
steps = [DetectFitConfig(...),
         FilterConfig(photons=(500.0, Inf)),
         FrameConnectConfig(max_frame_gap=5),
         FilterConfig(precision=(0.0, 0.005)),
         DriftConfig(degree=2),
         RenderConfig(zoom=20)]
```

## Multi-Dataset Workflows

Dataset boundaries are encoded in the data structure:

```julia
# Vector of arrays = multiple datasets
(result, info) = analyze([dataset1, dataset2, dataset3], config)

# Single array = one dataset
(result, info) = analyze(single_stack, config)

# File-based: MIC format auto-detects blocks as datasets
config = AnalysisConfig(
    camera = cam,
    steps = [DetectFitConfig(path="data.h5", h5_format=:mic), ...],
)
(result, info) = analyze(config)
```

## Multi-Target (Multi-Color)

Each channel runs its own `AnalysisConfig` pipeline, then cross-channel steps (alignment, composite rendering) run via `AbstractMultiTargetStep` dispatch:

```julia
mt = MultiTargetConfig(
    labels = [:IgG, :C1q],
    colors = [:cyan, :magenta],
    steps = [
        CompositeRenderConfig(zoom=20.0, strategy=GaussianRender()),
        CrossAlignConfig(method=:entropy),
        CompositeRenderConfig(zoom=20.0, strategy=GaussianRender()),
    ],
    outdir = "output/cell1/",
)

(result, info) = analyze([
    (images_647, config_647),
    (images_568, config_568),
], mt)

result[:IgG].smld    # Per-channel access
result.smlds         # All SMLDs
```

## Output and Provenance

Every `analyze()` call returns `(result, info)`:

| Field | Description |
|-------|-------------|
| `result.smld` | Final BasicSMLD with corrected localizations |
| `result.smld_connected` | Connected SMLD with track info |
| `result.drift_model` | Fitted drift model |
| `info.elapsed_s` | Total wall time |
| `info.steps[:detectfit]` | Typed info struct from upstream package |
| `info.step_infos` | Vector of StepInfos with per-step timing, config, and summary |

When `outdir` is set, each step writes to `outdir/01_detectfit/`, `outdir/02_filter/`, etc. with saved configs, summary stats, and diagnostic plots.

## JuliaSMLM Ecosystem

```
SMLMData (core types: Emitter, Camera, BasicSMLD)
    +-- SMLMBoxer (ROI detection)
    +-- GaussMLE (GPU-accelerated MLE fitting)
    +-- SMLMFrameConnection (linking across frames)
    +-- SMLMDriftCorrection (entropy-based drift correction)
    +-- SMLMRender (super-resolution rendering)
    +-- SMLMSim (simulation + image generation)
    +-- MicroscopePSFs (PSF models)
    +-- SMLMAnalysis (integrates all)
```

All packages share [SMLMData.jl](https://github.com/JuliaSMLM/SMLMData.jl) types. Coordinates are in microns throughout.

## Documentation

- [Stable docs](https://JuliaSMLM.github.io/SMLMAnalysis.jl/stable/) - Full guide, configuration reference, and API
- [API Overview](api_overview.md) - LLM-parseable API reference

## License

MIT License - see [LICENSE](LICENSE) file.
