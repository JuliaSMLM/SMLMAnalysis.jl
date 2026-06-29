```@meta
CurrentModule = SMLMAnalysis
```

# Drift Correction

Sample drift — the slow wander of the stage during a long acquisition — blurs a
super-resolution reconstruction just as motion blur smears a photograph. The
drift-correction step estimates the drift trajectory and subtracts it from every
localization. It is selected by a `DriftConfig` and backed by
[SMLMDriftCorrection](https://github.com/JuliaSMLM/SMLMDriftCorrection.jl), which
corrects drift **fiducial-free** — no beads required.

```julia
analyze(smld, DriftConfig(degree = 2))   # → (corrected_smld, StepInfo)
```

## When to use / prerequisites

- Run on a `BasicSMLD` of localizations (after [detection/fitting](@ref
  "Detection & Fitting"), and usually after a [quality filter](@ref
  "Quality Filter") and [frame connection](@ref "Frame Connection")).
- Drift is estimated from the data's own redundancy — the same fixed structures
  reappear across frames — so it needs enough localizations spread over the
  acquisition to register against. Very sparse data drifts-correct poorly.

## Inputs, returns & artifacts

- **Input:** the current `smld`.
- **Returns:** `(corrected_smld, StepInfo)`. The fitted drift model is on
  `StepInfo.info.model`, and the pipeline also surfaces it on
  [`AnalysisResult`](@ref)`.drift_model`.
- **Artifacts** (when `outdir` is set): `drift_trajectory.png` (X/Y drift vs.
  frame and the XY path), `stats.md` (max drift, inter-dataset shifts, entropy,
  convergence), and at `DETAILED` verbosity a `per_dataset.md` breakdown. A
  corrected SMLD is checkpointed to `smld_corrected.jld2` at `Checkpoint.EXPENSIVE`.

## Concept

The algorithm models the reconstruction as a Gaussian mixture and searches for
the drift trajectory that makes the localizations of fixed structures pile up as
tightly as possible — i.e. the trajectory that **minimizes the entropy** of the
rendered cloud (Cnossen et al. 2021). The trajectory itself is expanded in a
**Legendre-polynomial** basis over normalized time, so `degree` sets how
wiggly a drift path the fit can follow. This refines the fiducial-free,
redundancy-based approach of Wester et al. 2021. For the full derivation see the
SMLMDriftCorrection documentation.

## Two modes

`dataset_mode` chooses how the datasets within an SMLD relate (see
[Multi-Dataset](@ref "Multi-Dataset Acquisitions")):

- **`:continuous`** — one long, uninterrupted movie. For long acquisitions
  (≳4000 frames) split the fit into chunks (`chunk_frames`) or use a single
  higher-degree polynomial for shorter ones.

  ```julia
  DriftConfig(degree = 3, dataset_mode = :continuous, chunk_frames = 4000)
  ```

- **`:registered`** — multiple independent acquisitions of the same field of view
  with stage registration between them (e.g. SeqSRM). Each dataset is corrected
  independently and then aligned to the others; this requires spatial overlap
  between datasets.

  ```julia
  DriftConfig(degree = 2, dataset_mode = :registered)
  ```

## Configuration

`DriftConfig` is re-exported from SMLMDriftCorrection (the "upstream owns the
config" idiom). The fields you set most often:

| field | typical | meaning |
|---|---|---|
| `degree` | `2`–`5` | polynomial degree of the drift trajectory; higher = more flexible |
| `dataset_mode` | `:registered` | `:continuous` (one movie) or `:registered` (independent, overlapping acquisitions) |
| `quality` | `:singlepass` | `:singlepass` (fast) or `:iterative` (refines to convergence) |
| `chunk_frames` | `0` | for continuous mode, frames per chunk (≈4000 is reasonable); `0` = single polynomial |
| `n_chunks` | `0` | alternative to `chunk_frames`: fixed number of chunks |
| `auto_roi` | `false` | estimate drift from a dense ROI subset (faster, stronger signal) |

See the SMLMDriftCorrection documentation for the complete field list and
defaults.

```julia
# Continuous 8000-frame acquisition, chunked
(corrected, info) = analyze(smld,
    DriftConfig(degree = 3, dataset_mode = :continuous,
                chunk_frames = 4000, auto_roi = true))

drift = info.info.model                       # the fitted drift model
SMLMDriftCorrection.drift_trajectory(drift)   # sampled trajectory (µm)
```

## Output & interpretation

The step's `StepInfo.summary` reports the headline numbers:

| field | meaning |
|---|---|
| `max_drift_nm` | largest intra-dataset drift excursion (nm) |
| `max_intershift_nm` | largest inter-dataset shift (nm), multi-dataset only |
| `dataset_mode`, `quality` | the mode/quality actually used |
| `converged`, `iterations` | for `:iterative` quality |
| `entropy`, `backend` | final entropy and compute backend |

Sanity checks: a sensible `max_drift_nm` is tens to a few hundred nm over a long
acquisition; the `drift_trajectory.png` should look like a smooth wander, not
noise. The step **warns** (at `PROGRESS`) when inter-dataset shifts exceed
~500 nm — usually a sign the wrong `dataset_mode` was chosen (e.g.
`:registered` on data acquired as one continuous movie).

## Notes & caveats

- **Mode mismatch is the most common mistake.** Continuous data corrected as
  `:registered` (or vice versa) produces large spurious inter-dataset shifts;
  heed the warning and switch modes.
- **Frame numbering is per-dataset**, which is what lets the Legendre basis
  normalize each dataset's time to ``[-1, 1]`` — see
  [Multi-Dataset](@ref "Multi-Dataset Acquisitions").
- **Threads help.** Per-dataset correction is parallelized; start Julia with
  `-t auto`.

## References

- J. Cnossen, T. J. Cui, C. Joo, C. Smith. "Drift correction in localization
  microscopy using entropy minimization." *Optics Express* **29**, 27961 (2021).
  [doi:10.1364/OE.426620](https://doi.org/10.1364/OE.426620)
- M. J. Wester, *et al.* "Robust, fiducial-free drift correction for
  super-resolution imaging." *Scientific Reports* **11**, 23672 (2021).
  [doi:10.1038/s41598-021-02850-7](https://doi.org/10.1038/s41598-021-02850-7)

See the [SMLMDriftCorrection documentation](https://github.com/JuliaSMLM/SMLMDriftCorrection.jl)
for the algorithm in full and all configuration options.
