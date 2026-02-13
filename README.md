# SMLMAnalysis

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSMLM.github.io/SMLMAnalysis.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/SMLMAnalysis.jl/dev/)
[![Build Status](https://github.com/JuliaSMLM/SMLMAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSMLM/SMLMAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSMLM/SMLMAnalysis.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSMLM/SMLMAnalysis.jl)

Complete SMLM (Single Molecule Localization Microscopy) analysis pipeline: detection, fitting, filtering, frame connection, drift correction, density filtering, and super-resolution rendering. Orchestrates the [JuliaSMLM](https://github.com/JuliaSMLM) ecosystem into reproducible workflows with provenance tracking. All coordinates are in microns.

## Installation

```julia
using Pkg
Pkg.add("SMLMAnalysis")
```

## Quick Start

### Pipeline with AnalysisConfig

```julia
using SMLMAnalysis

# Define camera (512x512 pixels, 100nm pixel size)
cam = IdealCamera(512, 512, 0.1)

# Configure and run the full pipeline
config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(
            boxer=BoxerConfig(boxsize=9, psf_sigma=0.130),
            fitter=GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20)),
        FilterConfig(photons=(500.0, Inf), precision=(0.0, 0.007)),
        FrameConnectConfig(max_frame_gap=5),
        CalibrationConfig(),
        DriftConfig(degree=2),
        RenderConfig(zoom=20, colormap=:inferno),
    ],
    outdir = "output/",
)
(result, info) = analyze(image_stacks, config)

result.smld               # Final BasicSMLD with corrected localizations
result.drift_model        # Drift model (if DriftConfig was used)
info.steps[:detectfit]    # Per-step info from upstream packages
info.elapsed_s            # Total wall time
```

### Step-by-step with analyze() dispatch

Each step uses `analyze()` with a typed config:

```julia
using SMLMAnalysis

cam = IdealCamera(512, 512, 0.1)

# 1. Detect ROIs and fit localizations (GPU-accelerated)
(smld, df_info) = analyze(image_stacks, DetectFitConfig(
    camera=cam, boxer=BoxerConfig(boxsize=9, psf_sigma=0.130)))

# 2. Filter by quality
(smld, _) = analyze(smld, FilterConfig(photons=(500.0, Inf)))

# 3. Link localizations across frames
(smld, fc_info) = analyze(smld, FrameConnectConfig(max_frame_gap=5))

# 4. Calibrate localization uncertainties
(smld, cal_info) = analyze(smld, CalibrationConfig())

# 5. Correct sample drift
(smld, dc_info) = analyze(smld, DriftConfig(degree=2))

# 6. Render super-resolution image
(img, _) = analyze(smld, RenderConfig(zoom=20, colormap=:inferno))

# Save/load intermediate results for parameter iteration
save_smld("checkpoint.h5", smld)
smld = load_smld("checkpoint.h5")
```

## Composable Pipeline

The `steps` vector is composable — after `DetectFitConfig` (which must be first since it produces localizations from raw images), you can use any combination, order, or repetition of steps.

**Minimal pipeline** — detect+fit and render, skipping all intermediate processing:

```julia
config = AnalysisConfig(
    camera = cam,
    steps = [DetectFitConfig(boxer=BoxerConfig(boxsize=9)), RenderConfig(zoom=20)],
)
```

**Multiple renders** — different zooms or colormaps in one pipeline:

```julia
steps = [
    DetectFitConfig(boxer=BoxerConfig(boxsize=9)),
    FilterConfig(photons=(500.0, Inf)),
    DriftConfig(degree=2),
    RenderConfig(zoom=10, colormap=:viridis),
    RenderConfig(zoom=20, colormap=:inferno),
]
```

**Repeated filtering** — filter, process, filter again:

```julia
steps = [
    DetectFitConfig(boxer=BoxerConfig(boxsize=9)),
    FilterConfig(photons=(500.0, Inf)),           # coarse filter
    FrameConnectConfig(max_frame_gap=5),
    CalibrationConfig(),
    FilterConfig(precision=(0.0, 0.005)),          # tighter filter after connection
    DriftConfig(degree=2),
    RenderConfig(zoom=20),
]
```

**Config provenance**: `DetectFitConfig` and `FilterConfig` are defined in SMLMAnalysis. `FrameConnectConfig` and `DriftConfig` are re-exported from upstream packages (SMLMFrameConnection, SMLMDriftCorrection). `RenderConfig` is re-exported from SMLMRender.

**Extensibility**: Define a new `struct YourConfig <: AbstractSMLMConfig`, implement `analyze(smld, YourConfig)`, and add it to the steps vector.

### Rendered output

![Super-resolution render](docs/src/assets/render_gaussian.png)

## Configuration Reference

### DetectFitConfig

Combined ROI detection (SMLMBoxer) and MLE fitting (GaussMLE). Embeds native upstream configs for full access to all detection/fitting parameters.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `boxer` | `BoxerConfig(boxsize=11, psf_sigma=0.135)` | Detection config (SMLMBoxer) |
| `fitter` | `GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20)` | Fitting config (GaussMLE) |
| `camera` | `nothing` | Camera model (required for step-by-step; injected by AnalysisConfig) |
| `path` | `nothing` | H5 file path (file-based workflow) |
| `paths` | `nothing` | Vector of H5 paths (one per dataset) |
| `h5_format` | `:auto` | H5 format: `:auto`, `:smart`, `:mic` |

**BoxerConfig** key fields: `boxsize` (ROI size in pixels), `min_photons` (detection threshold), `psf_sigma` (expected PSF width in microns), `backend` (`:auto`/`:gpu`/`:cpu`).

**GaussMLEConfig** key fields: `psf_model` (PSF type, e.g. `GaussianXYNBS()` for variable width), `iterations` (MLE iterations), `backend`, `constraints`, `batch_size`.

**PSF model guidance:** `GaussianXYNBS()` (default) fits per-emitter PSF width -- use for most data. `GaussianXYNB(sigma)` uses a known PSF sigma (faster). `GaussianXYNBSXSY()` fits independent sigma_x and sigma_y (for astigmatic 3D).

### FilterConfig

Quality-based filtering. All criteria use `(min, max)` tuples; `nothing` disables.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `photons` | `nothing` | Photon count range, e.g. `(500.0, Inf)` |
| `precision` | `nothing` | Localization precision range (microns), e.g. `(0.0, 0.007)` |
| `pvalue` | `nothing` | Goodness-of-fit p-value range, e.g. `(1e-3, 1.0)` |
| `psf_sigma` | `nothing` | PSF width filter: `:auto` (mode +/- 10%) or `(min, max)` |

### FrameConnectConfig

Links localizations of the same emitter across consecutive frames. Re-exported from SMLMFrameConnection.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_frame_gap` | `5` | Maximum dark frames allowed in a track |
| `max_sigma_dist` | `5.0` | Spatial matching threshold (in sigma units) |

### CalibrationConfig

Uncertainty calibration from frame connection analysis. Uses connected tracks to estimate motion blur and CRLB scale factors.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `clamp_k_to_one` | `true` | Prevent CRLB scale factor k < 1 |
| `filter_high_chi2` | `false` | Remove tracks with high chi-squared pairs |
| `chi2_filter_threshold` | `6.0` | Chi-squared threshold for track removal |

### DriftConfig

Entropy-based drift correction. Re-exported from SMLMDriftCorrection.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `degree` | `2` | Legendre polynomial degree |
| `dataset_mode` | `:registered` | `:registered` for independent multi-dataset, `:continuous` for one long acquisition |
| `quality` | `:singlepass` | `:singlepass` or `:iterative` |
| `chunk_frames` | `0` | Split continuous data into chunks of N frames (0 = no chunking) |
| `auto_roi` | `true` | Use dense ROI subset for faster estimation |
| `maxn` | `200` | Maximum neighbors for entropy calculation |
| `max_iterations` | `100` | Maximum iterations for optimization |

**Mode guidance:** Use `dataset_mode=:registered` (default) when datasets are independent acquisitions of the same FOV. Use `dataset_mode=:continuous` when data is one long acquisition split across files, with `chunk_frames=4000` for long acquisitions.

### DensityFilterConfig

Removes isolated localizations by neighbor density.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_sigma` | `2.0` | Search radius in localization uncertainty units |
| `min_neighbors` | `:auto` | Minimum neighbor count; `:auto` uses valley detection |

### RenderConfig

Super-resolution image rendering via SMLMRender.

```julia
RenderConfig(zoom=20, colormap=:inferno)                    # Gaussian render (default)
RenderConfig(zoom=20, strategy=HistogramRender())            # Histogram render
RenderConfig(zoom=20, strategy=CircleRender())               # Circle render
RenderConfig(zoom=20, color_by=:absolute_frame, colormap=:turbo)  # Temporal coloring
```

## Multi-Dataset Workflows

Dataset boundaries are encoded in the data structure:

```julia
# Vector of 3D arrays = multiple datasets
image_stacks = [dataset1, dataset2, dataset3]  # each is (height, width, frames)
(result, info) = analyze(image_stacks, config)

# Single 3D array = one dataset
(result, info) = analyze(single_stack, config)

# File-based: MIC format auto-detects blocks as datasets
config = AnalysisConfig(
    camera = cam,
    steps = [DetectFitConfig(path="data.h5", h5_format=:mic), ...],
)
(result, info) = analyze(config)  # No data argument needed
```

## Output Format

`analyze()` returns `(AnalysisResult, AnalysisInfo)`.

| Output | Description |
|--------|-------------|
| `result.smld` | Final drift-corrected BasicSMLD |
| `result.smld_connected` | Connected SMLD with track info (if FrameConnectConfig used) |
| `result.drift_model` | Fitted drift model (if DriftConfig used) |
| `info.elapsed_s` | Total wall time (seconds) |
| `info.steps[:detectfit]` | Per-step upstream info struct |
| `info.step_records` | Vector of StepRecords with timing, config, and summary stats |

When `outdir` is set, each step writes to `outdir/01_detectfit/`, `outdir/02_filter/`, etc., with `config.toml`, `stats.md`, and diagnostic plots.

## I/O

```julia
# Save/load SMLD checkpoints (HDF5 format)
save_smld("results.h5", smld; drift_model=dm, source_file="raw_data.h5")
smld = load_smld("results.h5")
smld_info("results.h5")  # Print file summary without loading

# Load microscope data
images, info = load_smart_h5("smart_data.h5")           # SMART microscope format
images, info = load_mic_h5("mic_data.h5")              # MIC format
block_images = load_mic_h5_block("mic_data.h5", 1)  # Load single block
```

## Multi-Target (Multi-Color)

```julia
mt_config = MultiTargetConfig(
    labels = [:IgG, :C1q],
    colors = [:red, :green],
    render_zoom = 20,
    outdir = "output/cell1/",
)

(result, info) = analyze([
    (images_647, config_647),
    (images_568, config_568),
], mt_config)

result[:IgG].smld    # Per-channel access
result.smlds         # All SMLDs
```

## JuliaSMLM Ecosystem

```
SMLMData (core types)
    +-- SMLMBoxer (ROI detection)
    +-- GaussMLE (GPU-accelerated MLE fitting)
    +-- SMLMFrameConnection (linking across frames)
    +-- SMLMDriftCorrection (entropy-based drift correction)
    +-- SMLMRender (super-resolution rendering)
    +-- SMLMSim (simulation + image generation)
    +-- MicroscopePSFs (PSF models)
    +-- SMLMAnalysis (integrates all)
```

## Related Packages

- **[SMLMData.jl](https://github.com/JuliaSMLM/SMLMData.jl)** - Core types: Emitter, Camera, BasicSMLD
- **[SMLMBoxer.jl](https://github.com/JuliaSMLM/SMLMBoxer.jl)** - ROI detection from raw images
- **[GaussMLE.jl](https://github.com/JuliaSMLM/GaussMLE.jl)** - GPU-accelerated MLE fitting
- **[SMLMFrameConnection.jl](https://github.com/JuliaSMLM/SMLMFrameConnection.jl)** - Linking localizations across frames
- **[SMLMDriftCorrection.jl](https://github.com/JuliaSMLM/SMLMDriftCorrection.jl)** - Entropy-based drift correction
- **[SMLMRender.jl](https://github.com/JuliaSMLM/SMLMRender.jl)** - Super-resolution image rendering
- **[SMLMSim.jl](https://github.com/JuliaSMLM/SMLMSim.jl)** - SMLM data simulation
- **[MicroscopePSFs.jl](https://github.com/JuliaSMLM/MicroscopePSFs.jl)** - PSF models

## License

MIT License - see [LICENSE](LICENSE) file for details.
