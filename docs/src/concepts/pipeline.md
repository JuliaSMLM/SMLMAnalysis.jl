```@meta
CurrentModule = SMLMAnalysis
```

# The Pipeline Model

SMLMAnalysis has one organizing idea: **an analysis is an ordered list of typed
step configs, and Julia's method dispatch routes each one to its
implementation.** Understanding that idea explains the whole package — why steps
compose freely, why `analyze()` is the only verb you call, and why adding a step
requires no change to the orchestrator.

![The analyze() pipeline: a single call threads the data through an ordered list of typed step configs, each returning (smld, info); optional steps drop in by adding their config.](../assets/pipeline_schematic.svg)

## A pipeline is a vector of step configs

A localization pipeline is a sequence of transformations on a localization set.
In SMLMAnalysis each transformation is a **step**, and each step is described by
a typed configuration struct that subtypes `SMLMData.AbstractSMLMConfig`:

```julia
config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(boxer = BoxerConfig(boxsize = 9, psf_sigma = 0.130),
                        fitter = GaussMLEConfig(psf_model = GaussianXYNBS())),
        FilterConfig(photons = (500.0, Inf)),
        FrameConnectConfig(max_frame_gap = 5),
        DriftConfig(degree = 2),
        RenderConfig(zoom = 20),
    ],
    outdir = "output/",
)
```

`AnalysisConfig.steps` is just a `Vector{AbstractSMLMConfig}`. The config carries
*what* you want done and the parameters for each step; it carries no logic.

## Dispatch is the router

The single verb is [`analyze`](@ref). Calling `analyze(state, cfg)` dispatches on
the types of **both** arguments — the current data `state` and the step `cfg` —
to the right method. There is no step registry, no string lookup, and no
`if/elseif` chain in the pipeline loop. The orchestrator is a plain fold:

```julia
state = image_stacks
for (i, cfg) in enumerate(steps)
    cfg = _prepare_step(cfg, camera)         # e.g. inject the camera into DetectFitConfig
    (state, step_info) = analyze(state, cfg; outdir, step_number = i, verbose)
    push!(step_infos, step_info)
end
```

Three consequences follow directly:

- **Composability.** The loop does not care what the steps are — only that
  `analyze()` has a method for `(typeof(state), typeof(cfg))`. So steps can be
  reordered, repeated, or omitted freely.
- **Extensibility.** Adding a step needs no change to the orchestrator. Define a
  config type and an `analyze()` method and the step works — even from *another
  package or your own script*. See [Extending the Pipeline](@ref).
- **Error clarity.** A wrong order is a `MethodError`, not a silent wrong answer.
  `FilterConfig` before any `DetectFitConfig` gives
  `no method matching analyze(::Vector{…}, ::FilterConfig)` — because filtering
  expects a `BasicSMLD`, and only detection produces one from raw images.

## The tuple pattern: `(result, info)`

Every `analyze()` call returns a **2-tuple** following the JuliaSMLM convention:

```julia
(result, info) = analyze(data, config)   # whole pipeline
(smld,   info) = analyze(smld, cfg)       # one step
```

- The **result** becomes the input state for the next step. For most steps it is
  a `BasicSMLD`; for the whole pipeline it is an [`AnalysisResult`](@ref).
- The **info** is the provenance: per-step it is a [`StepInfo`](@ref) wrapping the
  upstream package's own info struct; for the whole pipeline it is an
  [`AnalysisInfo`](@ref) aggregating them all. Nothing is hidden in globals — if
  a step computed it, it is on the returned info.

This is why you can run a pipeline in one call *or* step by step and get exactly
the same per-step information either way.

## State threading

The orchestrator threads one working `smld` through the steps: each step
receives the current localizations and returns the updated set. Two pieces of
state are captured from specific steps and surfaced on the final result:

- `result.smld_connected` — the frame-connected tracks, from `FrameConnectConfig`.
- `result.drift_model` — the fitted drift trajectory, from `DriftConfig`.

The first step is special only in that it manufactures the initial `smld`:

```
DetectFitConfig (required first)
        │ produces the initial smld from raw images
        ▼
  ┌────────────────────────────────────┐
  │ FilterConfig          (0+ times)   │
  │ IntensityFilterConfig (0+ times)   │   any order,
  │ FrameConnectConfig    (0–1)        │   any combination —
  │ DriftConfig           (0–1)        │   dispatch doesn't
  │ DensityFilterConfig   (0+ times)   │   care, as long as
  │ BaGoLConfig           (0–1)        │   each step's input
  │ cluster configs       (0+ times)   │   type is satisfied
  │ RenderConfig          (0+ times)   │
  └────────────────────────────────────┘
```

The only hard constraint is that a `DetectFitConfig` comes first, because it is
the only step that turns raw image stacks into a `BasicSMLD`; every other step
consumes and returns an existing `smld`.

## Repeatable and optional steps

Because dispatch — not a registry — drives the loop, multiplicity is free:

- **`FilterConfig`**, **`RenderConfig`**, **`DensityFilterConfig`**,
  **`IntensityFilterConfig`**, and the **cluster configs** can appear any number
  of times (a coarse filter early and a tighter one after connection; several
  renders at different zooms; cluster then render the labels).
- **`FrameConnectConfig`** and **`DriftConfig`** are typically used once, but this
  is convention, not an enforced rule.

Uncertainty calibration is *not* a separate step — it is configured inside
`FrameConnectConfig` via its `calibration =` field, runs as part of frame
connection, and its results appear on `FrameConnectInfo.calibration`.

## Two layers per step

Each step is implemented in two layers, which is worth knowing when reading the
source or [adding your own](@ref "Extending the Pipeline"):

1. An **internal work function** (not exported) that does the computation and
   returns `(result, <StepInfo>)` — either implemented in SMLMAnalysis (e.g.
   `filter_step`) or delegated to an upstream package (e.g.
   `SMLMDriftCorrection.driftcorrect`).
2. A thin **`analyze()` method** that times the call, wraps the result in a
   [`StepInfo`](@ref), and returns the tuple. Its signature ends in `kwargs...`
   so the orchestrator can pass context keywords a step may ignore.

## Where configs come from

Some step configs are defined in SMLMAnalysis; others are re-exported from the
upstream package that owns them (the "upstream owns the config" idiom — a `const`
alias, dispatched on directly).

| Config | Defined in | Notes |
|--------|-----------|-------|
| `DetectFitConfig` | SMLMAnalysis | wraps SMLMBoxer + GaussMLE |
| `FilterConfig` | SMLMAnalysis | native quality filtering |
| `IntensityFilterConfig` | SMLMAnalysis | native Poisson multi-emitter rejection |
| `DensityFilterConfig` | SMLMAnalysis | native neighbor filtering |
| `FrameConnectConfig` / `CalibrationConfig` | SMLMFrameConnection | re-exported (`CalibrationConfig` is a sub-config, not a step) |
| `DriftConfig` | SMLMDriftCorrection | re-exported |
| `RenderConfig` | SMLMRender | re-exported |
| `BaGoLConfig` | SMLMBaGoL | re-exported; step wraps `run_bagol` |
| `DBSCANConfig`, `VoronoiConfig`, `HopkinsConfig`, … | SMLMClustering | re-exported; steps wrap `cluster` / `cluster_statistics` |

See [The JuliaSMLM Ecosystem](@ref) for what each upstream package does, and
[Data Model & Provenance](@ref) for the types that flow between steps.
