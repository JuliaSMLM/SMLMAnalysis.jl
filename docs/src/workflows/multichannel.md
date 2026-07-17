```@meta
CurrentModule = SMLMAnalysis
```

# Multi-Channel Analysis

A multi-color experiment images two or more targets in separate channels — for
example an IgG variant and complement C1q, each on its own dye. This page shows
how to analyze every channel with its own pipeline and then run cross-channel
steps (composite renders, alignment, cross-correlation) that operate on all
channels together. The orchestrator is `MultiTargetConfig`; if you are new to the
single-channel flow, read [The Pipeline Model](@ref) and [Running a Pipeline](@ref)
first.

## The `MultiTargetConfig`

`MultiTargetConfig` is the single object that describes a multi-channel run — the
channels, their colors, and the cross-channel steps to run once each channel is
analyzed. It is the orchestrator you construct and hand to `analyze`:

```@docs
MultiTargetConfig
```

Its fields in practical terms:

| Field | Type | Meaning |
|-------|------|---------|
| `labels` | `Vector{Symbol}` | Channel names, e.g. `[:IgG, :C1q]`. Must be unique and match the number of channels. |
| `colors` | `Vector{Symbol}` | One color per channel. Defaults to cyan/magenta for 2 channels, CMY for 3 (up to 6; pass explicit colors beyond that). |
| `steps` | `Vector{AbstractMultiTargetStep}` | Ordered cross-channel steps (composite renders, alignment, cross-correlation). |
| `outdir` | `String` | Root output directory. |
| `verbose` | `Int` | Verbosity level (default `Verbosity.STANDARD`). |

```julia
mt = MultiTargetConfig(
    labels = [:IgG, :C1q],
    colors = [:cyan, :magenta],
    steps = [
        CompositeRenderConfig(zoom = 20.0, strategy = GaussianRender()),
        CrossAlignConfig(method = :entropy),
        CompositeRenderConfig(zoom = 20.0, strategy = GaussianRender()),  # post-alignment
        CrossCorrConfig(r_max = 0.5, dr = 0.005),
    ],
    outdir = "output/cell1/",
)
```

The `steps` run in order, threading a `Vector{BasicSMLD}` through. Most steps pass
the SMLDs through unchanged; [Cross-Alignment](@ref) replaces them with aligned
copies, so a second composite render placed after it shows the corrected overlay.

## The shape of a multi-channel run

A multi-channel analysis has two phases, both described by the one
`MultiTargetConfig` above and run by a single `analyze` call:

1. **Per-channel pipelines.** Each channel is analyzed independently with its own
   [`AnalysisConfig`](@ref) — exactly the pipeline you would run on a single color.
2. **Cross-channel steps.** Once every channel has produced a `BasicSMLD`, the
   ordered multi-target steps run across the resulting `Vector{BasicSMLD}`.

## Running it

Call `analyze` with a vector of `(data, AnalysisConfig)` tuples — one per channel,
in the same order as `labels` — together with the `MultiTargetConfig`:

```julia
(result, info) = analyze([
    (image_stacks_647, config_647),
    (image_stacks_568, config_568),
], mt)
```

Each `data` is whatever the single-channel `analyze` accepts (an image stack, a
vector of stacks, or a file path). Each channel's `AnalysisConfig` carries its own
camera, steps, ROI, and verbosity; the orchestrator reuses all of those but
redirects the channel's output under `outdir/<label>/`, so you do not set a
per-channel `outdir` yourself.

## Working with the result

`analyze` returns a `(MultiTargetResult, MultiTargetInfo)` tuple, following the
JuliaSMLM convention (see [Data Model & Provenance](@ref)).

`MultiTargetResult` fields and indexing:

```julia
result.smlds        # Vector{BasicSMLD}, one per channel (aligned if CrossAlign ran)
result.channels     # Dict{Symbol, AnalysisResult}
result.step_infos   # Vector{StepInfo} for the cross-channel steps
result.outdir       # root output directory
result[:IgG]        # AnalysisResult for one channel (indexing == result.channels[:IgG])
result[:IgG].smld   # that channel's final SMLD
keys(result)        # the channel labels, in order
```

`MultiTargetInfo` carries the metadata:

```julia
info.elapsed_s      # total wall time, seconds
info.channels[:IgG] # AnalysisInfo for one channel's pipeline
info.step_infos     # Vector{StepInfo} for the cross-channel steps
stepinfo(info, :crossalign).info  # look up a cross-channel step by name (searches step_infos, not channels)
```

## Cross-channel steps

Each multi-target step is a `<: AbstractMultiTargetStep` config dispatched on the
channel `Vector{BasicSMLD}`. They are documented on their own pages — add them to
`MultiTargetConfig.steps` in the order you want them to run:

- **[Composite Render](@ref)** — overlay the channels into a single multi-color
  image, each tinted by its `colors` entry (per-step `colors` override the
  `MultiTargetConfig` default). Accepts the same strategies as [Rendering](@ref "Rendering")
  (`GaussianRender`, `HistogramRender`, `CircleRender`). Pass-through.

- **[Cross-Alignment](@ref)** — remove residual channel-to-channel shift
  (chromatic offset, registration error) by entropy or FFT cross-correlation.
  State-modifying: it returns aligned SMLDs, so steps after it see the corrected
  positions.

- **[Cross-Correlation](@ref)** — quantify co-localization between a pair of
  channels via the pair cross-correlation `g(r)`. Pass-through.

## Output layout

The run mirrors the two-phase structure on disk: a subdirectory per channel, a
`composite/` directory for the cross-channel steps, and per-channel SMLD files
plus a serialized config at the root.

```
output/cell1/
├── IgG/                    # per-channel pipeline output (01_detectfit/, 02_filter/, ...)
├── C1q/
├── composite/              # cross-channel step outputs, numbered in step order
│   ├── 01_compositerender/
│   ├── 02_crossalign/
│   ├── 03_compositerender/
│   ├── 04_crosscorr/
│   └── README.md           # color scheme, per-channel counts, step summary
├── smld_IgG.h5             # per-channel saved SMLDs (with drift model)
├── smld_C1q.h5
└── multi_target_config.toml
```
