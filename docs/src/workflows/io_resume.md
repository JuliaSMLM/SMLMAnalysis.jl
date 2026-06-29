```@meta
CurrentModule = SMLMAnalysis
```

# I/O & Resume

This page covers getting data in and results out: importing raw microscope
acquisitions, saving and reloading localizations, persisting full pipeline state
for cross-session resume, and the per-step checkpoints that let you iterate on
later steps without re-running expensive ones.

## Saving and loading localizations (HDF5)

A finished (or intermediate) `BasicSMLD` is saved to a self-describing HDF5 file
with [`save_smld`](@ref) and read back with [`load_smld`](@ref):

```julia
save_smld("results.h5", result.smld;
          source_file = "/data/experiment.h5",   # recorded in /provenance
          drift_model = result.drift_model)       # optional, embedded in the file
smld = load_smld("results.h5")
```

The file stores emitter columns (positions and uncertainties **in microns**,
photons, background, frame/dataset/track ids), the camera and its calibration,
and optional provenance. `save_smld` also takes `compression` (HDF5 deflate
level 0–9, default `3`). `load_smld` reconstructs the correct emitter type
automatically, including the GaussMLE σ-fitting types (`Emitter2DFitSigma`,
`Emitter2DFitSigmaXY`) and 3D (`Emitter3DFit`).

To inspect a file without loading all emitters, use [`smld_info`](@ref), which
prints the format version, emitter type, counts, PSF model, save timestamp, and
available emitter fields:

```julia
smld_info("results.h5")
```

## Saving and resuming full pipeline state (JLD2)

HDF5 holds one `BasicSMLD`. To snapshot the **whole** pipeline — final SMLD, the
raw SMLD, the frame-connected tracks, the drift model, and step provenance — use
the JLD2-backed [`save_pipeline_state`](@ref) / [`load_pipeline_state`](@ref):

```julia
(result, info) = analyze(image_stacks, config)
save_pipeline_state("output/pipeline.jld2", result;
                    step_infos = info.step_infos,
                    camera     = config.camera)
```

`load_pipeline_state` returns a `NamedTuple` you can read fields from or feed back
into [`analyze`](@ref) to continue from where you left off:

```julia
state = load_pipeline_state("output/pipeline.jld2")
state.smld             # final BasicSMLD
state.smld_connected   # frame-connected tracks (or nothing)
state.drift_model      # fitted drift model (or nothing)
state.step_infos       # Vector{StepInfo}

(smld, _) = analyze(state.smld, FilterConfig(photons = (300.0, Inf)))
```

## Step checkpoints during a run

When `AnalysisConfig.outdir` is set, the pipeline can drop each step's output
`BasicSMLD` as a JLD2 file inside that step's subdirectory, so you can resume from
any step or sweep a later step's parameters without re-running the upstream work.
What gets written is controlled by the `checkpoint` level (see
[Running a Pipeline](@ref) for the level semantics and how to set it):

| Level | Constant | Writes |
|-------|----------|--------|
| 0 | `Checkpoint.NONE` | nothing |
| 1 | `Checkpoint.END` | only the final SMLD-producing step |
| 2 | `Checkpoint.EXPENSIVE` (default) | DetectFit, FrameConnect, Drift, BaGoL |
| 3 | `Checkpoint.ALL` | every SMLD-producing step (filters included) |

Each step uses a stable filename, so you always know where to look:

| Step | Checkpoint file | Written at |
|------|-----------------|------------|
| Detection & Fitting | `smld_raw.jld2` | `EXPENSIVE` |
| Frame Connection | `smld_combined.jld2` | `EXPENSIVE` |
| Drift Correction | `smld_corrected.jld2` (drift model embedded) | `EXPENSIVE` |
| BaGoL | `smld_bagol.jld2` | `EXPENSIVE` |
| Filter | `smld_filtered.jld2` | `ALL` |
| Intensity Filter | `smld_intensity.jld2` | `ALL` |
| Density Filter | `smld_density.jld2` | `ALL` |
| Clustering | `smld_clustered.jld2` | `ALL` |

The SMLD is stored under the `smld` key, so reload it with JLD2 and pick up the
pipeline from there:

```julia
using JLD2
smld = JLD2.load("output/01_detectfit/smld_raw.jld2")["smld"]   # full BasicSMLD
(smld, _) = analyze(smld, FilterConfig(photons = (500.0, Inf)))
```

The default `EXPENSIVE` level guarantees no costly step ever runs producing only
figures — its SMLD is always on disk for downstream iteration.

## Importing microscope data

### SMART microscope H5

The SMART format keeps the image stack under `/Main/data`. Load the full stack
(or a frame range), peek at metadata only, or load a properly oriented array:

```julia
images = load_smart_h5("acquisition.h5")                 # (width, height, frames)
images = load_smart_h5("acquisition.h5"; frame_range = 1:1000)
info   = load_smart_h5_info("acquisition.h5")            # width, height, nframes, dtype, file_size_gb
data, info = smart_h5_to_array("acquisition.h5"; max_frames = 1000)  # transposed to (rows, cols, frames)
```

### MIC H5

The MIC (MATLAB Instrument Control) format stores one or more data **blocks**,
each treated as a dataset, plus optional per-pixel sCMOS calibration. Pixel size
and QE are not stored and must be supplied when building a camera:

```julia
images, dataset_indices = load_mic_h5("experiment.h5")   # images (h×w×frames), block index per frame
info  = load_mic_h5_info("experiment.h5")                # n_frames, n_blocks, frames_per_block, has_calibration
block = load_mic_h5_block("experiment.h5", 1)            # one block, memory-efficient
cam   = build_camera_from_mic_h5("experiment.h5"; pixel_size = 0.1, qe = 0.9)  # SCMOSCamera from calibration
```

### Loading directly in `DetectFitConfig`

For large acquisitions you usually do not load images into memory yourself.
[`DetectFitConfig`](@ref) reads files directly via `path` (one file) or `paths`
(one file per dataset), with `h5_format = :auto` (default), `:smart`, or `:mic`;
MIC blocks are auto-detected as datasets. With `h5_format = :mic`, `camera =
nothing`, and `pixel_size`/`qe` set, the camera is built from the MIC calibration
automatically.

```julia
config = AnalysisConfig(
    camera = cam,
    steps  = [DetectFitConfig(path = "experiment.h5", h5_format = :mic,
                              boxer = BoxerConfig(boxsize = 9)), RenderConfig(zoom = 20)],
)
(result, info) = analyze(config)   # no in-memory data argument needed
```

## Output directory layout

With `outdir` set, every step writes a numbered subdirectory —
`outdir/01_detectfit/`, `outdir/02_filter/`, `outdir/05_driftcorrect/`, … — each
containing `config.toml` (the exact config used), `info.toml` (scalar fields of
the upstream info struct), `stats.md` (a human-readable summary), any
verbosity-gated figures, and the checkpoint SMLD described above. A top-level
`outdir/config.toml` records the camera, ROI, verbosity, checkpoint level, and the
full ordered step list. See [Data Model & Provenance](@ref) for what each holds.
