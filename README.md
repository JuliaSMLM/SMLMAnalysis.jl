# SMLMAnalysis

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSMLM.github.io/SMLMAnalysis.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/SMLMAnalysis.jl/dev/)
[![Build Status](https://github.com/JuliaSMLM/SMLMAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSMLM/SMLMAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSMLM/SMLMAnalysis.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSMLM/SMLMAnalysis.jl)

SMLM analysis pipeline for the [JuliaSMLM](https://github.com/JuliaSMLM) ecosystem: detection, fitting, filtering, frame connection, drift correction, grouping, clustering, and super-resolution rendering with provenance tracking.

## Analysis steps & primary literature

Every step is a `…Config` you drop into an `AnalysisConfig` (or call standalone via `analyze(state, cfg)`). The method each implements comes from the primary literature cited here; full references with DOIs are in [References](#references).

| Step (config) | What it does | Primary reference(s) |
|---|---|---|
| **Detection** (`DetectFitConfig`→`BoxerConfig`) | Finds candidate emitters / ROIs in the raw frames | Huang 2013 [¹](#references) |
| **Fitting** (`DetectFitConfig`→`GaussMLEConfig`) | Maximum-likelihood localization reaching the Cramér–Rao bound; 2D/3D-astigmatic PSF models | Smith 2010 [²](#references); Huang 2013 [¹](#references) |
| **Filter** (`FilterConfig`) | Quality cuts on photons, background, precision, track length, and the χ²/LLR goodness-of-fit **p-value** | Huang 2011 [³](#references) |
| **Intensity filter** (`IntensityFilterConfig`) | Rejects multi-emitter events via a Poisson upper-tail test against a fitted excitation-field model | *(native to this package)* |
| **Density filter** (`DensityFilterConfig`) | Removes isolated localizations by local k-nearest-neighbor density | *(standard practice)* |
| **Frame connection** (`FrameConnectConfig`) | Links repeated blinks of one fluorophore across frames (spatiotemporal LAP) | Schodt 2021 [⁴](#references) |
| **Drift correction** (`DriftConfig`) | Fiducial-free sample-drift estimation by entropy minimization (DME) | Cnossen 2021 [⁵](#references); Wester 2021 [⁶](#references) |
| **BaGoL** (`BaGoLConfig`) | Bayesian Grouping of Localizations — RJMCMC grouping for sub-nm precision | Fazel 2022 [⁷](#references) |
| **Clustering** (`DBSCANConfig`, `HDBSCANConfig`, `HierarchicalConfig`, `VoronoiConfig`) | Groups localizations into clusters | DBSCAN: Ester 1996 [⁸](#references); HDBSCAN: Campello 2013 [⁹](#references); Voronoi/SR-Tesseler: Levet 2015 [¹⁰](#references) |
| **Spatial statistics** (`HopkinsConfig`, `VoronoiDensityConfig`) | Clustering-tendency and local-density statistics | Hopkins 1954 [¹¹](#references); Levet 2015 [¹⁰](#references) |
| **Cross-correlation** (`CrossCorrConfig`) | Pair-correlation *g(r)* between two channels (co-localization) | Sengupta 2011 [¹²](#references); Veatch 2012 [¹³](#references) |
| **Edge classification** (`edgeclassify`) | Labels localizations as cell interior / membrane / outside | *(native to this package)* |
| **Render / Composite** (`RenderConfig`, `CompositeRenderConfig`) | Super-resolution image; multi-color composite of aligned channels | *(visualization)* |
| **Cross-channel align** (`CrossAlignConfig`) | Registers color channels (entropy / FFT cross-correlation) | Wester 2021 [⁶](#references) |

## Coming soon

Planned steps wrapping in-progress JuliaSMLM packages (config surface not yet exposed here):

| Method | Package | Reference |
|---|---|---|
| **Fourier Ring Correlation** resolution | SMLMResolution | Nieuwenhuizen 2013 [¹⁴](#coming-soon-references) |
| **BaMF** — Bayesian Multi-emitter Fitting (RJMCMC) | SMLMBaMF | Fazel 2019 [¹⁵](#coming-soon-references) |
| **PSF learning** — data-driven PSF modeling | PSFLearning | Liu 2024 [¹⁶](#coming-soon-references) |
| **Deep-learning SR** — U-Net dense localization (DECODE-style) | SMLMDeepFit | Speiser 2021 [¹⁷](#coming-soon-references) |

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

## Loading Data

`image_stacks` above is your **raw camera data** — pixel intensities (ADU), not
localizations. It is one or more 3D arrays of shape `(height, width, frames)`, each
array being one *dataset* (a continuous acquisition). Pass a `Vector` of them for a
multi-dataset experiment, or a single 3D array for one:

```julia
image_stacks = [stack1, stack2]   # 2 datasets, each an (height, width, frames) array
image_stacks = stack              # 1 dataset — a single 3D array is also accepted
```

There are three ways to get it.

**1. Load an H5 file into memory.** The package ships loaders for two microscope
formats (both return `(height, width, frames)`):

```julia
# SMART — one continuous acquisition → one array, wrapped as a 1-element Vector:
stack, _ = smart_h5_to_array("data/experiment.h5")
(result, info) = analyze([stack], config)

# MIC (LidkeLab) — the file holds multiple blocks, each a separate dataset. Load them
# as separate stacks so the block boundaries are preserved (analyzing [stack] from
# load_mic_h5 would merge every block into one dataset):
n = load_mic_h5_info("data/experiment.h5").n_blocks
image_stacks = [load_mic_h5_block("data/experiment.h5", i) for i in 1:n]
(result, info) = analyze(image_stacks, config)
```

**2. Stream from files — no in-memory array.** Point `DetectFitConfig` at the
file(s) and call `analyze(config)` with no data argument. MIC blocks are
auto-detected as separate datasets. This is the memory-efficient path for large
acquisitions:

```julia
# For MIC files you can build the camera from the file's own calibration:
cam = build_camera_from_mic_h5("data/experiment.h5"; pixel_size=0.1)

config = AnalysisConfig(camera = cam, steps = [
    DetectFitConfig(path="data/experiment.h5", h5_format=:mic,
                    boxer=BoxerConfig(boxsize=9, psf_sigma=0.130)),
    FilterConfig(photons=(500.0, Inf)),
    RenderConfig(zoom=20),
])
(result, info) = analyze(config)      # each dataset is loaded from disk in turn
```

Multiple files, one dataset each: `DetectFitConfig(paths=["d1.h5", "d2.h5"], ...)`.

**3. Any other source.** Any `Array{<:Real,3}` of `(height, width, frames)` works —
a TIFF stack read with [TiffImages.jl](https://github.com/tlnagy/TiffImages.jl), or
simulated data from `SMLMSim` (`gen_images`). See
[`examples/loading_data.jl`](examples/loading_data.jl) for a runnable walkthrough of
all three, and the other [examples](examples/) for full simulated pipelines.

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
(smld, _) = analyze(smld, RenderConfig(zoom=20, colormap=:inferno))  # writes an image when an outdir is set; use SMLMRender.render(smld, cfg) for the image in memory

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
| `stepinfo(info, :detectfit).info` | Typed info struct from upstream package |
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

## References

1. Huang, F., Hartwich, T.M.P., Rivera-Molina, F.E. *et al.* "Video-rate nanoscopy using sCMOS camera-specific single-molecule localization algorithms." *Nature Methods* **10**, 653–658 (2013). [doi:10.1038/nmeth.2488](https://doi.org/10.1038/nmeth.2488)
2. Smith, C.S., Joseph, N., Rieger, B., Lidke, K.A. "Fast, single-molecule localization that achieves theoretically minimum uncertainty." *Nature Methods* **7**, 373–375 (2010). [doi:10.1038/nmeth.1449](https://doi.org/10.1038/nmeth.1449)
3. Huang, F., Schwartz, S.L., Byars, J.M., Lidke, K.A. "Simultaneous multiple-emitter fitting for single molecule super-resolution imaging." *Biomedical Optics Express* **2**(5), 1377–1393 (2011). [doi:10.1364/BOE.2.001377](https://doi.org/10.1364/BOE.2.001377)
4. Schodt, D.J., Lidke, K.A. "Spatiotemporal Clustering of Repeated Super-Resolution Localizations via Linear Assignment Problem." *Frontiers in Bioinformatics* **1**, 724325 (2021). [doi:10.3389/fbinf.2021.724325](https://doi.org/10.3389/fbinf.2021.724325)
5. Cnossen, J., Cui, T.J., Joo, C., Smith, C. "Drift correction in localization microscopy using entropy minimization." *Optics Express* **29**(18), 27961–27974 (2021). [doi:10.1364/OE.426620](https://doi.org/10.1364/OE.426620)
6. Wester, M.J., Schodt, D.J., Mazloom-Farsibaf, H., Fazel, M., Pallikkuth, S., Lidke, K.A. "Robust, fiducial-free drift correction for super-resolution imaging." *Scientific Reports* **11**, 23672 (2021). [doi:10.1038/s41598-021-02850-7](https://doi.org/10.1038/s41598-021-02850-7)
7. Fazel, M., Wester, M.J., Schodt, D.J. *et al.* "High-precision estimation of emitter positions using Bayesian grouping of localizations." *Nature Communications* **13**, 7152 (2022). [doi:10.1038/s41467-022-34894-2](https://doi.org/10.1038/s41467-022-34894-2)
8. Ester, M., Kriegel, H.-P., Sander, J., Xu, X. "A density-based algorithm for discovering clusters in large spatial databases with noise." *Proc. 2nd Int. Conf. Knowledge Discovery and Data Mining (KDD-96)*, 226–231 (1996).
9. Campello, R.J.G.B., Moulavi, D., Sander, J. "Density-based clustering based on hierarchical density estimates." *PAKDD 2013*, LNCS **7819**, 160–172. [doi:10.1007/978-3-642-37456-2_14](https://doi.org/10.1007/978-3-642-37456-2_14)
10. Levet, F., Hosy, E., Kechkar, A. *et al.* "SR-Tesseler: a method to segment and quantify localization-based super-resolution microscopy data." *Nature Methods* **12**, 1065–1071 (2015). [doi:10.1038/nmeth.3579](https://doi.org/10.1038/nmeth.3579)
11. Hopkins, B., Skellam, J.G. "A new method for determining the type of distribution of plant individuals." *Annals of Botany* **18**(2), 213–227 (1954). [doi:10.1093/oxfordjournals.aob.a083391](https://doi.org/10.1093/oxfordjournals.aob.a083391)
12. Sengupta, P., Jovanovic-Talisman, T., Skoko, D., Renz, M., Veatch, S.L., Lippincott-Schwartz, J. "Probing protein heterogeneity in the plasma membrane using PALM and pair correlation analysis." *Nature Methods* **8**, 969–975 (2011). [doi:10.1038/nmeth.1704](https://doi.org/10.1038/nmeth.1704)
13. Veatch, S.L., Machta, B.B., Shelby, S.A., Chiang, E.N., Holowka, D.A., Baird, B.A. "Correlation functions quantify super-resolution images and estimate apparent clustering due to over-counting." *PLoS ONE* **7**(2), e31457 (2012). [doi:10.1371/journal.pone.0031457](https://doi.org/10.1371/journal.pone.0031457)

### Coming-soon references

14. Nieuwenhuizen, R.P.J., Lidke, K.A., Bates, M. *et al.* "Measuring image resolution in optical nanoscopy." *Nature Methods* **10**, 557–562 (2013). [doi:10.1038/nmeth.2448](https://doi.org/10.1038/nmeth.2448)
15. Fazel, M., Wester, M.J., Mazloom-Farsibaf, H. *et al.* "Bayesian multiple emitter fitting using reversible jump Markov chain Monte Carlo." *Scientific Reports* **9**, 13791 (2019). [doi:10.1038/s41598-019-50232-x](https://doi.org/10.1038/s41598-019-50232-x)
16. Liu, S., Chen, J., Hellgoth, J. *et al.* "Universal inverse modelling of point spread functions for SMLM localization and microscope characterization." *Nature Methods* **21**, 1082–1093 (2024). [doi:10.1038/s41592-024-02282-x](https://doi.org/10.1038/s41592-024-02282-x)
17. Speiser, A., Müller, L.-R., Hoess, P. *et al.* "Deep learning enables fast and dense single-molecule localization with high accuracy." *Nature Methods* **18**, 1082–1090 (2021). [doi:10.1038/s41592-021-01236-x](https://doi.org/10.1038/s41592-021-01236-x)

## License

MIT License - see [LICENSE](LICENSE) file.
