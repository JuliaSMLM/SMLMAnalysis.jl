```@meta
CurrentModule = SMLMAnalysis
```

# The JuliaSMLM Ecosystem

SMLMAnalysis sits at the top of [JuliaSMLM](https://github.com/JuliaSMLM), a set
of focused packages that each own one part of the SMLM analysis problem. This
page is the map: what each package does, which step it backs, and when you might
call it directly instead of through the pipeline.

All packages share the core data types from **SMLMData** (`BasicSMLD`,
`Emitter2DFit`, cameras), so values pass between them without conversion.

## Dependency map

```
SMLMData ............ core types (Emitter, Camera, BasicSMLD, ROIBatch) — no deps
   │
   ├── SMLMBoxer ............ ROI detection from raw frames
   ├── GaussMLE ............. GPU-accelerated MLE PSF fitting
   ├── SMLMFrameConnection .. link blinks across frames (+ uncertainty calibration)
   ├── SMLMDriftCorrection .. fiducial-free drift correction
   ├── SMLMRender ........... super-resolution image rendering
   ├── SMLMBaGoL ............ Bayesian grouping of localizations (RJMCMC)
   ├── SMLMClustering ....... DBSCAN / Hierarchical / Voronoi / Hopkins
   ├── SMLMSim .............. simulation + image generation
   └── MicroscopePSFs ....... PSF models (Gaussian, Airy, …)
                              │
                       SMLMAnalysis (integrates all of the above)
```

## What each package does

| Package | Backs the step | Role | Primary reference |
|---------|----------------|------|-------------------|
| **SMLMBoxer** | [Detection](@ref "Detection & Fitting") | Finds candidate molecules and cuts ROIs from raw frames | Huang et al. 2013 |
| **GaussMLE** | [Fitting](@ref "Detection & Fitting") | Maximum-likelihood Gaussian/astigmatic PSF fitting with CRLB uncertainties (GPU) | Smith et al. 2010 |
| **SMLMFrameConnection** | [Frame Connection](@ref) | Links repeated blinks of one fluorophore into single, higher-precision localizations; optional uncertainty calibration | Schodt & Lidke 2021 |
| **SMLMDriftCorrection** | [Drift Correction](@ref) | Fiducial-free drift estimation by entropy minimization | Cnossen et al. 2021; Wester et al. 2021 |
| **SMLMRender** | [Rendering](@ref "Rendering") | Turns localizations into super-resolution images (histogram / Gaussian / circle / ellipse) | — |
| **SMLMBaGoL** | [Bayesian Grouping](@ref "Bayesian Grouping (BaGoL)") | Groups localizations into true emitters via RJMCMC, beyond raw precision | Fazel et al. 2022 |
| **SMLMClustering** | [Clustering](@ref clustering-step) | Cluster labeling (DBSCAN/Hierarchical/Voronoi) and spatial statistics (Hopkins, Voronoi density) | Ester et al. 1996; Levet et al. 2015 |
| **SMLMSim** | — | Simulates SMLM data and generates synthetic image stacks (used throughout the docs) | — |
| **MicroscopePSFs** | — | PSF models used by fitting and simulation | — |

Full citations are on the [References](@ref references-page) page; each [step page](@ref
"Pipeline Steps: Overview") repeats its own primary reference and links to the
upstream package's documentation for the algorithm details.

## Steps SMLMAnalysis owns

Not every step comes from an upstream package. SMLMAnalysis implements several of
its own — the parts of a real pipeline that do not belong to any single
lower-level package:

- **[Quality Filter](@ref)** (`FilterConfig`) — threshold filtering on photons,
  precision, p-value, PSF width, and z.
- **[Intensity Filter](@ref)** (`IntensityFilterConfig`) — Poisson upper-tail
  rejection of multi-emitter events against an estimated excitation field.
- **[Density Filter](@ref)** (`DensityFilterConfig`) — removal of isolated
  localizations by neighbor count, with automatic threshold selection.
- **Multi-channel steps** — [Composite Render](@ref), [Cross-Alignment](@ref),
  and [Cross-Correlation](@ref) operate across colors and are built here on top
  of SMLMRender / SMLMDriftCorrection / NearestNeighbors.

## When to use a package directly

The pipeline is the right tool when you want a reproducible, multi-step analysis
with provenance and on-disk diagnostics. Reach past it, straight to an upstream
package, when:

- You need a **capability the pipeline doesn't expose** — an upstream function or
  option SMLMAnalysis doesn't surface. Because all packages share SMLMData types,
  you can pull `result.smld` out of a pipeline and pass it to any upstream
  function.
- You are **developing or debugging one algorithm** in isolation and want its
  package's full API and docs.
- You are doing **post-hoc exploration** on an already-analyzed SMLD (e.g.
  trying several clustering parameters) and don't need it recorded as a pipeline
  step.

Since `using SMLMAnalysis` re-exports the key upstream types and verbs (`cluster`,
`render`, `run_bagol`, `frameconnect`, the config types, the emitter types), you
can usually do both styles from the same session without extra imports.
