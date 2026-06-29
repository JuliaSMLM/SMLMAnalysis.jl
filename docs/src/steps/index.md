```@meta
CurrentModule = SMLMAnalysis
```

# Pipeline Steps: Overview

This section documents every analysis step — one page per `analyze()` step. Each
step is selected by the **type** of its config (see [The Pipeline Model](@ref)),
runs in the order you list it in `AnalysisConfig.steps`, and returns an updated
SMLD plus a [`StepInfo`](@ref). This page is the map: the catalog, the typical
order, and where each method comes from in the literature.

## Catalog

| Step | Config | Backed by | What it does |
|------|--------|-----------|--------------|
| [Detection & Fitting](@ref "Detection & Fitting") | `DetectFitConfig` | SMLMBoxer + GaussMLE | Detect ROIs in raw frames and fit a PSF to each → the initial localizations |
| [Quality Filter](@ref) | `FilterConfig` | SMLMAnalysis | Threshold on photons, precision, p-value, PSF width, z |
| [Intensity Filter](@ref) | `IntensityFilterConfig` | SMLMAnalysis | Reject multi-emitter events via a Poisson upper-tail test |
| [Frame Connection](@ref) | `FrameConnectConfig` | SMLMFrameConnection | Link repeated blinks into single, higher-precision localizations |
| [Drift Correction](@ref) | `DriftConfig` | SMLMDriftCorrection | Estimate and remove sample drift (fiducial-free) |
| [Density Filter](@ref) | `DensityFilterConfig` | SMLMAnalysis | Drop isolated localizations by neighbor count |
| [Rendering](@ref "Rendering") | `RenderConfig` | SMLMRender | Render localizations to a super-resolution image |
| [Bayesian Grouping](@ref "Bayesian Grouping (BaGoL)") | `BaGoLConfig` | SMLMBaGoL | Group localizations into true emitters (RJMCMC) |
| [Clustering](@ref clustering-step) | `DBSCANConfig`, `VoronoiConfig`, `HopkinsConfig`, … | SMLMClustering | Label clusters / compute spatial statistics |
| [Composite Render](@ref) | `CompositeRenderConfig` | SMLMAnalysis + SMLMRender | Multi-channel RGB overlay |
| [Cross-Alignment](@ref) | `CrossAlignConfig` | SMLMAnalysis + SMLMDriftCorrection | Register channels to a common frame |
| [Cross-Correlation](@ref) | `CrossCorrConfig` | SMLMAnalysis + NearestNeighbors | Pair cross-correlation g(r) between channels |

The last three operate on **multiple channels** (`Vector{BasicSMLD}`) and are
usually driven by [`MultiTargetConfig`](@ref) — see [Multi-Channel](@ref
"Multi-Channel Analysis").

## A typical pipeline

A standard single-color dSTORM / DNA-PAINT analysis runs roughly in this order:

```julia
steps = [
    DetectFitConfig(boxer = BoxerConfig(boxsize = 9, psf_sigma = 0.130),
                    fitter = GaussMLEConfig(psf_model = GaussianXYNBS())),
    FilterConfig(photons = (500.0, Inf)),     # cut weak / spurious fits
    FrameConnectConfig(max_frame_gap = 5),     # link repeated blinks
    DriftConfig(degree = 2),                   # remove drift
    RenderConfig(zoom = 20, colormap = :inferno),
]
```

Common variations:

- Add an [Intensity Filter](@ref) early to reject multi-emitter events.
- Add a [Density Filter](@ref) after drift correction to clean isolated noise.
- Replace/augment the final render with [Bayesian Grouping](@ref "Bayesian Grouping (BaGoL)")
  for counting and maximum precision, then render the grouped emitters.
- Add a [Clustering](@ref clustering-step) step to label structures, then render colored by label.

Steps are composable (see [The Pipeline Model](@ref)): reorder, repeat, or omit
them as the experiment requires.

## Primary literature

Each method's primary reference, as published by the backing package. Full
citations are on the [References](@ref references-page) page; each step page repeats its own
reference and links to the upstream documentation for algorithm details.

| Step | Method | Primary reference |
|------|--------|-------------------|
| Detection | Gaussian-filter local-max + sCMOS-aware thresholding | Huang et al., *Nat. Methods* **10**, 653 (2013) |
| Fitting | MLE Gaussian/astigmatic PSF; CRLB precision | Smith et al., *Nat. Methods* **7**, 373 (2010); Mortensen et al., *Nat. Methods* **7**, 377 (2010) |
| Frame connection | Spatiotemporal LAP clustering of repeated blinks | Schodt & Lidke, *Front. Bioinform.* (2021) |
| Drift correction | Fiducial-free entropy minimization (Legendre basis) | Cnossen et al., *Opt. Express* **29**, 27961 (2021); Wester et al., *Sci. Rep.* **11**, 23672 (2021) |
| Bayesian grouping | Collapsed reversible-jump MCMC grouping | Fazel et al., *Nat. Commun.* **13**, 7152 (2022) |
| Clustering | DBSCAN; Voronoi tessellation (SR-Tesseler) | Ester et al., *KDD-96* (1996); Levet et al., *Nat. Methods* **12**, 1065 (2015) |

The native SMLMAnalysis steps — [Quality Filter](@ref), [Intensity Filter](@ref),
[Density Filter](@ref), and the multi-channel steps — are described in full on
their own pages, including the method and any reference.
