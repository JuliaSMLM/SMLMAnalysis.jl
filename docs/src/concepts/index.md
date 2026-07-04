```@meta
CurrentModule = SMLMAnalysis
```

# Concepts: Overview

SMLMAnalysis is an **integration layer**: it does little numerical work itself,
and instead composes the JuliaSMLM ecosystem's specialized packages into one
reproducible pipeline. These pages explain the few ideas you need to use it well
and to reason about what a pipeline is doing. Read them top to bottom, or jump to
a topic:

- [The Pipeline Model](@ref) — how `analyze()` dispatch turns a list of step
  configs into an analysis, the `(result, info)` tuple pattern, and state
  threading. **Start here.**
- [The JuliaSMLM Ecosystem](@ref) — what each upstream package does, which steps
  SMLMAnalysis owns itself, and when to reach for a package directly.
- [Data Model & Provenance](@ref) — the `BasicSMLD` localization container and
  emitter types that flow between steps, the camera and micron coordinate system,
  and the `StepInfo` / `AnalysisInfo` records that capture what happened.
- [The Multi-Target Model](@ref) — how the single-channel pipeline generalizes to
  multi-color experiments: a `Vector{BasicSMLD}` state and cross-channel steps.

For day-to-day tasks (installing, running, multi-dataset and multi-channel
acquisitions, I/O, adding a step) see the [Workflows](@ref "Installation & Setup")
section; for the per-step reference see [Pipeline Steps](@ref "Pipeline Steps: Overview").

## The mental model in one paragraph

You hand `analyze()` some data and a description of the analysis. The description
is an [`AnalysisConfig`](@ref): a camera, an output directory, and an ordered
`steps` vector of typed configs. The pipeline applies each step in turn, passing
the localization set from one to the next, and returns the final localizations
([`AnalysisResult`](@ref)) together with a complete provenance record
([`AnalysisInfo`](@ref)). Each step is backed by one ecosystem package (or by
SMLMAnalysis itself), selected purely by the config's type.

## Glossary

- **Localization** — a single fitted blinking event: a position with a fitted
  uncertainty (and photons, background, …). Stored as an `Emitter2DFit` /
  `Emitter3DFit`. Positions are in **microns** throughout.
- **SMLD** — Single-Molecule Localization Data: the `SMLMData.BasicSMLD`
  container holding a vector of emitters plus the camera, frame count, and
  dataset count. This is the value that flows between steps.
- **Camera** — an `AbstractCamera` (`IdealCamera` or `SCMOSCamera`) that defines
  the pixel grid in microns and so the coordinate system every step works in. Set
  once on the [`AnalysisConfig`](@ref); see
  [The camera & coordinate system](@ref).
- **Emitter** — after [Bayesian grouping](@ref "Bayesian Grouping (BaGoL)") the
  word also means a *true fluorophore* reconstructed from many localizations; in
  raw data "emitter" and "localization" coincide.
- **Step** — one transformation in the pipeline, named by its config type
  (`DetectFitConfig`, `FilterConfig`, …). Each step has an `analyze(state, cfg)`
  method.
- **Config** — a typed struct describing a step's parameters
  (`<: SMLMData.AbstractSMLMConfig`). Configs carry parameters, not logic.
- **Info** — the second element of every `analyze()` return tuple: a
  [`StepInfo`](@ref) per step (wrapping the upstream package's own info struct),
  aggregated into an [`AnalysisInfo`](@ref) for the whole run.
- **Dataset** — an independent acquisition block within one SMLD (e.g. a cell, an
  ROI, or a registered segment). Frames are numbered *per dataset*. See
  [Multi-Dataset](@ref "Multi-Dataset Acquisitions").
- **Channel / target** — one color in a multi-color experiment. Several SMLDs are
  analyzed together by [`MultiTargetConfig`](@ref). See
  [The Multi-Target Model](@ref) for the concept and
  [Multi-Channel](@ref "Multi-Channel Analysis") for the workflow.
- **Tuple pattern** — the JuliaSMLM convention that analysis functions return
  `(result, info)` rather than mutating arguments or hiding state in globals.
- **CRLB** — the Cramér–Rao lower bound, the theoretical minimum localization
  uncertainty; the per-localization `σ` reported by the [fitter](@ref
  "Detection & Fitting") is a CRLB estimate, and several steps use it as a weight.
