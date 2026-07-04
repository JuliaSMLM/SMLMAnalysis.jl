```@meta
CurrentModule = SMLMAnalysis
```

# Detection & Fitting

This is the first step of every pipeline: it turns raw camera frames into
localizations. It detects candidate molecules and cuts a small region of interest
(ROI) around each ([SMLMBoxer](https://github.com/JuliaSMLM/SMLMBoxer.jl)), then
fits a point-spread-function model to every ROI to recover a sub-pixel position,
photon count, background, and per-localization uncertainty
([GaussMLE](https://github.com/JuliaSMLM/GaussMLE.jl)). It is selected by a
`DetectFitConfig`.

```julia
analyze(images, DetectFitConfig(boxer = BoxerConfig(boxsize = 9, psf_sigma = 0.130),
                                fitter = GaussMLEConfig(psf_model = GaussianXYNBS())))
# → (smld, StepInfo)
```

## When to use / prerequisites

- **Always first.** Detection & fitting is the only step that produces a
  `BasicSMLD` from raw image data; every other step consumes and returns an
  existing `smld`. Placing another step before it is a `MethodError` (see
  [The Pipeline Model](@ref)).
- Needs a **camera**. In a full pipeline the camera is injected automatically
  from `AnalysisConfig.camera`; for a standalone `analyze(images, cfg)` call set
  `camera =` on the config (or let it be built from a MIC H5 file, below).

## Inputs, returns & artifacts

- **Input** — the raw image stack(s), supplied one of three ways:
  1. **In-memory**: pass the images to `analyze(images, cfg)` — a single 3D array
     (one dataset) or a `Vector` of 3D arrays (multiple datasets).
  2. **One file, split into datasets**: `DetectFitConfig(path = "data.h5",
     dataset_frames = [...])`.
  3. **Multiple files, one per dataset**: `DetectFitConfig(paths = ["d1.h5", "d2.h5"])`.
- **Returns** — `(smld, StepInfo)`; the `StepInfo.info` is a `DetectFitInfo`
  (`n_rois`, `n_fits`, `n_datasets`, `n_frames_per_dataset`, …).
- **Artifacts** (when `outdir` is set) — a detection/fit overlay on sample frames
  and fit-quality diagnostics; at `DEBUG` verbosity, per-frame detection movies.

## Concept

Two operations run back to back:

1. **Detection (SMLMBoxer).** The frames are filtered to enhance
   diffraction-limited spots and thresholded; local maxima become candidate
   molecules, and a `boxsize`×`boxsize` ROI is cut around each.
2. **Fitting (GaussMLE).** Each ROI is fit by **maximum likelihood** to a PSF
   model, returning the sub-pixel position, photons, background, and the
   Cramér–Rao-bound (CRLB) uncertainty of each — the `σ_x`/`σ_y` that the rest of
   the pipeline relies on (see [Data Model & Provenance](@ref)).

The detailed detection and fitting algorithms live in the SMLMBoxer and GaussMLE
documentation; this step composes them and handles datasets, files, and cameras.

## Configuration

`DetectFitConfig` embeds the two upstream configs and adds the data-source and
camera options:

| field | default | meaning |
|---|---|---|
| `boxer` | `BoxerConfig(boxsize=11, psf_sigma=0.135)` | detection config (ROI size, expected PSF σ in µm) |
| `fitter` | `GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20)` | PSF-fit config (model, iterations) |
| `camera` | `nothing` | the camera; injected by `AnalysisConfig`, or set for standalone use |
| `path` / `paths` | `nothing` | file-based sources (single file / one file per dataset) |
| `dataset_frames` | `nothing` | frame ranges that split a single file into datasets |
| `datasets` | `nothing` | subset of source slots to include |
| `h5_format` | `:auto` | `:auto`, `:smart`, or `:mic` (see [I/O & Resume](@ref)) |
| `pixel_size`, `qe` | `nothing`, `1.0` | build an SCMOS camera from MIC H5 calibration when `camera` is unset |

PSF models (from GaussMLE) include `GaussianXYNB`, `GaussianXYNBS` (free width),
`GaussianXYNBSXSY` (elliptical), and `AstigmaticXYZNB` (3D).

```julia
config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(
            boxer  = BoxerConfig(boxsize = 9, psf_sigma = 0.130),
            fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), iterations = 20)),
        # … subsequent steps …
    ],
)
(result, info) = analyze(image_stacks, config)
info.steps[:detectfit]    # the DetectFitInfo
```

## Output & interpretation

`DetectFitInfo` reports how detection and fitting went:

| field | meaning |
|---|---|
| `n_rois` | candidate ROIs found by detection |
| `n_fits` | localizations surviving the fit |
| `n_datasets`, `n_frames_per_dataset` | dataset/frame structure of the result |

Sanity checks: `n_fits` close to `n_rois` means most detections fit cleanly; a
large gap suggests `boxsize`/`psf_sigma` are mismatched to the data or the
detection threshold is too low. The overlay figure should show boxes centered on
real spots, not noise.

## Notes & caveats

- **Per-dataset processing.** Detection and fitting loop over datasets
  individually, so arbitrarily large acquisitions can be processed
  memory-efficiently. Frame numbers are **per dataset** — see
  [Multi-Dataset](@ref "Multi-Dataset Acquisitions").
- **GPU.** GaussMLE fits on the GPU (CUDA) when available; see its docs for
  GPU/CPU behavior and the per-model parameter sets.
- **Tune `boxsize` and `psf_sigma` to your optics** — they are the parameters
  that most affect detection quality.

## References

- **Detection.** F. Huang, *et al.* "Video-rate nanoscopy using sCMOS
  camera-specific single-molecule localization algorithms." *Nature Methods*
  **10**, 653–658 (2013).
  [doi:10.1038/nmeth.2488](https://doi.org/10.1038/nmeth.2488)
- **MLE fitting.** C. S. Smith, N. Joseph, B. Rieger, K. A. Lidke. "Fast,
  single-molecule localization that achieves theoretically minimum uncertainty."
  *Nature Methods* **7**, 373–375 (2010).
  [doi:10.1038/nmeth.1449](https://doi.org/10.1038/nmeth.1449)

See the [SMLMBoxer](https://github.com/JuliaSMLM/SMLMBoxer.jl) and
[GaussMLE](https://github.com/JuliaSMLM/GaussMLE.jl) documentation for the
detection and fitting algorithms and their full options.
