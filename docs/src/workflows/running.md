```@meta
CurrentModule = SMLMAnalysis
```

# Running a Pipeline

This page is a practical guide to running an analysis: the two ways to drive the
pipeline, the `analyze()` entry forms, how to read what comes back, and the two
dials — verbosity and checkpointing — that control how much is written to disk.
For the *why* behind the design, see [The Pipeline Model](@ref); for the catalog
of available steps, see [Pipeline Steps](@ref "Pipeline Steps: Overview").

## Two ways to drive the pipeline

`analyze()` is the only verb. You run a whole pipeline at once (config-driven) or
one step at a time (step-by-step); both call the same dispatch, so a step behaves
identically either way.

### Config-driven (run it all at once)

Build an [`AnalysisConfig`](@ref) holding the camera, an ordered `steps` vector,
and output settings, then hand it to `analyze()`:

```julia
config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(boxer = BoxerConfig(boxsize = 9, psf_sigma = 0.130),
                        fitter = GaussMLEConfig(psf_model = GaussianXYNBS())),
        FilterConfig(photons = (500.0, Inf)),
        FrameConnectConfig(max_frame_gap = 5),
        DriftConfig(degree = 2),
        RenderConfig(zoom = 20, colormap = :inferno),
    ],
    outdir = "output/",
)
(result, info) = analyze(image_stacks, config)
```

This is the form to reach for when the recipe is settled — it captures the entire
pipeline in one provenance-tracked object.

### Step-by-step (inspect between steps)

Call `analyze()` repeatedly, threading the working `smld` from one call into the
next. Each call returns `(smld, step_info::StepInfo)`, so you can inspect counts,
timing, and intermediate figures before deciding the next step's parameters. The
first call takes the image data and a [`DetectFitConfig`](@ref) (which carries its
own `camera`, since there is no `AnalysisConfig` to inject one); every later call
takes the `smld`:

```julia
(smld, info) = analyze(image_stacks,
    DetectFitConfig(camera = cam, boxer = BoxerConfig(boxsize = 9)))
@show length(smld.emitters)             # decide the photon cut from this

(smld, info) = analyze(smld, FilterConfig(photons = (500.0, Inf)))
(smld, info) = analyze(smld, DriftConfig(degree = 2))
(smld, info) = analyze(smld, RenderConfig(zoom = 20))
```

Because detection/fitting is the expensive step, a common pattern is to run it
once, save the `smld` (see [I/O & Resume](@ref)), then iterate cheaply on filter
and render parameters.

## `analyze()` entry forms

| Call | Use when |
|------|----------|
| `analyze(data, config::AnalysisConfig)` | The primary form; `data` is one image stack or a `Vector` of stacks (one per dataset). |
| `analyze(data, step1, step2, …; camera, outdir, …)` | Varargs convenience — builds the `AnalysisConfig` for you from positional step configs plus keywords. |
| `analyze(config::AnalysisConfig)` | File-based: no `data` argument. The `DetectFitConfig` carries `path=`/`paths=` and images are loaded from disk. |

The varargs form is handy for short, throwaway runs:

```julia
(result, info) = analyze(image_stacks,
    DetectFitConfig(boxer = BoxerConfig(boxsize = 9)),
    FilterConfig(photons = (500.0, Inf)),
    DriftConfig(degree = 2);
    camera = cam, outdir = "output/")
```

The file-based form keeps large acquisitions off the heap — the data argument is
omitted entirely:

```julia
config = AnalysisConfig(
    camera = cam,
    steps = [DetectFitConfig(path = "data.h5", h5_format = :mic,
                             boxer = BoxerConfig(boxsize = 9)),
             FilterConfig(photons = (500.0, Inf)),
             RenderConfig(zoom = 20)],
    outdir = "output/",
)
(result, info) = analyze(config)
```

## Reading the result

A full pipeline run returns an [`AnalysisResult`](@ref) and an [`AnalysisInfo`](@ref):

```julia
result.smld             # final BasicSMLD after all steps
result.smld_connected   # tracks from FrameConnectConfig (nothing if not run)
result.drift_model      # fitted drift polynomial from DriftConfig (nothing if not run)

info.elapsed_s          # total wall time, seconds
info.steps[:driftcorrect]  # upstream info struct keyed by step name
info.step_infos         # Vector{StepInfo}: per-step timing, config, summary
```

The emitters in `result.smld` carry positions and uncertainties **in microns**;
see [Data Model & Provenance](@ref) for the full type layout and what each step
contributes to the `(result, info)` tuple. In the step-by-step style the second
tuple element of each call is a single `StepInfo` rather than the aggregated
`AnalysisInfo`.

## Output directory

Setting `outdir` turns on disk output. Each step writes into its own numbered
subdirectory, `outdir/NN_stepname/`, alongside a top-level `config.toml` (the full
`AnalysisConfig`) and `summary.md` (the run timing/result table):

```
output/
├── config.toml          # AnalysisConfig: camera, ROI, verbosity, step manifest
├── summary.md           # per-step timing and result table
├── log.txt              # teed @info output
├── 01_detectfit/        # config.toml, info.toml, stats.md, figures
├── 02_filter/
├── 04_driftcorrect/
└── 08_render/
```

Each step subdir holds `config.toml` (that step's config), `info.toml` (scalar
fields of the upstream info struct), `stats.md`, and any figures gated by the
verbosity level. With `outdir = nothing` (the default) nothing is written and the
pipeline runs in-memory only. For saving and reloading the `smld` itself across
sessions, see [I/O & Resume](@ref).

## Verbosity levels

The `verbose` field (default `Verbosity.STANDARD`) controls how much diagnostic
output each step produces. Levels are cumulative — each adds to the ones below it:

| Level | Constant | Adds |
|-------|----------|------|
| 0 | `Verbosity.SILENT` | Errors only |
| 1 | `Verbosity.PROGRESS` | Step names, counts, timing |
| 2 | `Verbosity.STANDARD` | + `stats.md` and basic figures (fit quality, overlays, drift plots) |
| 3 | `Verbosity.DETAILED` | + diagnostic plots, per-filter breakdowns, localizations-per-frame |
| 4 | `Verbosity.DEBUG` | + MP4 animations, frame-by-frame, heavy visualization |

```julia
config = AnalysisConfig(camera = cam, steps = [...], verbose = Verbosity.DETAILED)

# or per call, step-by-step
(smld, info) = analyze(smld, RenderConfig(zoom = 20); verbose = Verbosity.DEBUG)
```

## Checkpoint levels

The `checkpoint` field (default `Checkpoint.EXPENSIVE`) controls which steps
persist their output `smld` to disk as a JLD2 checkpoint, so you can resume
downstream work without re-running upstream steps:

| Level | Constant | Writes |
|-------|----------|--------|
| 0 | `Checkpoint.NONE` | No SMLD checkpoints |
| 1 | `Checkpoint.END` | Only the final SMLD-producing step's output |
| 2 | `Checkpoint.EXPENSIVE` | Expensive steps (DetectFit, FrameConnect, Drift, BaGoL); cheap filters skipped *(default)* |
| 3 | `Checkpoint.ALL` | Every SMLD-producing step |

The `EXPENSIVE` default guarantees that no costly step ever runs leaving only
image/diagnostic output — its `smld` is always on disk for the next iteration.
Checkpoints require `outdir` to be set. See [I/O & Resume](@ref) for loading them
back.

```julia
config = AnalysisConfig(camera = cam, steps = [...], outdir = "output/",
                        checkpoint = Checkpoint.ALL)
```
