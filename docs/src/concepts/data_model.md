```@meta
CurrentModule = SMLMAnalysis
```

# Data Model & Provenance

Everything that flows between steps is built from a handful of types defined in
**SMLMData** and shared across the whole ecosystem. Knowing them makes the
`(result, info)` tuples and the on-disk outputs easy to read.

## Localizations: emitters

A single fitted blinking event is an **emitter**. The common 2D type is
`Emitter2DFit`, which carries a position, a photon count, a background, and the
*fitted uncertainties* of each:

| Field | Meaning |
|-------|---------|
| `x`, `y` | position, in **microns** |
| `photons`, `bg` | integrated photons and per-pixel background |
| `σ_x`, `σ_y` | localization precision (CRLB estimate), in microns |
| `σ_photons`, `σ_bg` | uncertainties of photons / background |
| `frame` | frame index *within its dataset* |
| `dataset` | which dataset (acquisition block) the emitter belongs to |
| `track_id` | track id assigned by [Frame Connection](@ref) (`0` if unconnected) |
| `id` | cluster id assigned by [Clustering](@ref clustering-step) (`0` = noise) |

3D data uses `Emitter3DFit` (adds `z`, `σ_z`); astigmatic and σ-fitting PSF
models add PSF-width fields. **All positions and uncertainties are in microns**
throughout the pipeline.

The fitted `σ_x`/`σ_y` are not incidental — they are the currency several steps
spend: [Frame Connection](@ref) weights links by them, the calibration adjusts
them, [Bayesian Grouping](@ref "Bayesian Grouping (BaGoL)") groups by them, and
[Density Filter](@ref) measures neighbor distances in their units.

## The container: `BasicSMLD`

A whole localization set is a `SMLMData.BasicSMLD`:

```julia
struct BasicSMLD
    emitters       # Vector of Emitter2DFit / Emitter3DFit / …
    camera         # the AbstractCamera (pixel geometry, gain, …)
    n_frames       # frames per dataset
    n_datasets     # number of datasets
    metadata       # Dict{String,Any} of provenance
end
```

This single value is the state threaded through the pipeline: each step takes a
`BasicSMLD` and returns a `BasicSMLD`. Note that `n_frames` is **per dataset**,
not the global total — see [Multi-Dataset](@ref "Multi-Dataset Acquisitions") for
why frame numbering is per-dataset.

## What `analyze()` returns

### `AnalysisResult`

The whole-pipeline result holds the final localizations plus the special state
captured along the way:

```julia
result.smld             # final BasicSMLD after all steps
result.smld_connected   # frame-connected tracks (or nothing)
result.drift_model      # fitted drift trajectory (or nothing)
```

### Provenance: `StepInfo` and `AnalysisInfo`

Provenance is first-class, never hidden in globals. Each step produces a
[`StepInfo`](@ref):

```julia
struct StepInfo
    number       # position in the pipeline
    name         # step name, e.g. "filter", "driftcorrect"
    config       # the exact config used
    timestamp
    elapsed_s
    summary      # Dict of headline stats (counts, rates, …)
    info         # the upstream package's own typed info struct (or nothing)
end
```

The whole run aggregates these into an [`AnalysisInfo`](@ref):

```julia
(result, info) = analyze(data, config)
info.elapsed_s              # total wall-clock time
info.step_infos             # Vector{StepInfo} — full ordered history
info.steps[:driftcorrect]   # the DriftInfo from the drift step
info.steps[:detectfit]      # the DetectFitInfo, etc.
```

`info.steps` is keyed by step name and holds each upstream package's own info
struct, so you reach a step's detailed results (e.g. `FrameConnectInfo`,
`DriftInfo`, `BaGoLInfo`) directly.

## On-disk outputs

When `AnalysisConfig.outdir` is set, each step writes a numbered subdirectory —
`outdir/02_filter/`, `outdir/05_driftcorrect/`, … — containing:

- `config.toml` — the step's configuration (reproducibility);
- `info.toml` — the scalar fields of the upstream info struct;
- `stats.md` — a human-readable summary;
- diagnostic figures (gated by [verbosity](@ref "Running a Pipeline"));
- optionally a checkpointed SMLD (`*.jld2`), gated by the checkpoint level.

The pipeline also writes a top-level `outdir/config.toml` capturing the camera,
ROI, verbosity, and the full ordered step list. For persisting and reloading
analysis state across sessions, see [I/O & Resume](@ref).
