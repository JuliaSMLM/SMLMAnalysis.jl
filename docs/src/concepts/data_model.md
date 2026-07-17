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

## The camera & coordinate system

Every `BasicSMLD` carries a **camera** — an `AbstractCamera` that defines the
pixel grid in physical space. The camera is the bridge between *pixel indices*
(what the detector records) and *microns* (what every downstream step works in):
it stores the pixel edge positions in microns, with the origin `(0, 0)` at the
top-left corner and pixel `(1, 1)` centered at `(pixel_size/2, pixel_size/2)`.
This is why every position and uncertainty in the pipeline is in microns — the
camera fixes the scale once, and detection and fitting express their results in it.

You set the camera **once**, on the [`AnalysisConfig`](@ref) (or per channel in a
[multi-target run](@ref "The Multi-Target Model")). The pipeline injects it into
`DetectFitConfig`, carries it on the resulting SMLD, and reuses it wherever a
pixel↔physical conversion is needed — detection, MLE fitting, and
[rendering](@ref "Rendering") (a render at `zoom = 20` produces an image 20× the
camera's pixel resolution). The camera also rides along through I/O, so a reloaded
SMLD still knows its geometry.

Two camera types are re-exported from SMLMData:

- **`IdealCamera(nx, ny, pixel_size)`** — a uniform pixel grid with Poisson-only
  noise. `pixel_size` is in microns: `IdealCamera(256, 128, 0.1)` is a 256×128
  sensor with 100 nm pixels. The right choice for EMCCD-like data and simulations.
- **`SCMOSCamera(nx, ny, pixel_size, readnoise; offset, gain, qe)`** — adds the
  per-pixel sCMOS noise model (read noise, offset, gain, quantum efficiency), each
  a scalar or a full per-pixel map. The MLE fitter uses these to weight each pixel
  correctly. [`build_camera_from_mic_h5`](@ref) constructs one straight from a MIC
  calibration file.

Pixel size is the master scale: it fixes how pixel coordinates map to microns, and
with them the meaning of every `σ`, every neighbor radius, and every render
dimension. Set it to match your objective and sensor.

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
info.elapsed_s                       # total wall-clock time
info.step_infos                      # Vector{StepInfo} — full ordered history
stepinfo(info, :driftcorrect).info   # the DriftInfo from the drift step
stepinfo(info, :detectfit).info      # the DetectFitInfo, etc.
```

[`stepinfo(info, name)`](@ref) searches `info.step_infos` for the first step with that
name (a `Symbol` or `String`) and returns its `StepInfo`; reach the upstream package's own
info struct (e.g. `FrameConnectInfo`, `DriftInfo`, `BaGoLInfo`) via `.info`. When a step
name repeats (e.g. several render steps), [`stepinfos(info, name)`](@ref) returns them all
in pipeline order.

!!! note "Migration from `info.steps`"
    The `info.steps::Dict{Symbol,Any}` field was removed (from both `AnalysisInfo` and
    `MultiTargetInfo`). Replace `info.steps[:x]` with `stepinfo(info, :x).info`. The old
    Dict returned the *last* occurrence of a duplicated step name and hid repeats behind
    suffixed keys like `:compositerender_3`; the accessors return the *first* match and
    expose every repeat via `stepinfos(info, :compositerender)`. `MultiTargetInfo.channels`
    is unchanged.

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
