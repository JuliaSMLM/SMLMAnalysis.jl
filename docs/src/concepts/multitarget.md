```@meta
CurrentModule = SMLMAnalysis
```

# The Multi-Target Model

A multi-color experiment images two or more targets — say an IgG variant and
complement C1q — each in its own channel. SMLMAnalysis handles these by
**generalizing the single-channel pipeline**, not by bolting on a separate system.
If you understand [The Pipeline Model](@ref), you already understand most of this.

## Two axes of generalization

The single-channel pipeline threads one `BasicSMLD` through a list of step configs,
dispatching `analyze(state, config)` on their types. The multi-target pipeline
keeps exactly that idea and changes only two things:

|         | Single-channel | Multi-target |
|---------|----------------|--------------|
| **State** | one `BasicSMLD` | a `Vector{BasicSMLD}` — one entry per channel |
| **Steps** | `<: AbstractSMLMConfig` | `<: AbstractMultiTargetStep` |

Everything else — the `analyze(state, step)` dispatch, the `(result, info)` tuple
per step, the threaded state — is unchanged. A cross-channel step is simply an
`analyze(::Vector{BasicSMLD}, ::AbstractMultiTargetStep)` method.

## Two phases

One [`MultiTargetConfig`](@ref) and one `analyze` call describe a run with two
phases:

1. **Per-channel pipelines.** Each channel is analyzed *independently* by a full
   single-channel pipeline — the orchestrator literally calls
   `analyze(data, AnalysisConfig)` once per channel, reusing that channel's own
   camera, steps, ROI, and verbosity. The resulting per-channel SMLDs become the
   channel vector.
2. **Cross-channel steps.** The ordered `AbstractMultiTargetStep`s then run over
   that `Vector{BasicSMLD}`, threading it from one step to the next.

Phase 1 shares no state between channels — each is exactly the analysis you would
run on a single color. All channel *interaction* happens in phase 2.

## Pass-through vs. state-modifying steps

As with single-channel steps, a cross-channel step either transforms the state or
just reads it for its outputs:

- **Pass-through** steps return the channel vector unchanged —
  [Composite Render](@ref) (the multi-color overlay image) and
  [Cross-Correlation](@ref) (the co-localization `g(r)`).
- **State-modifying** steps replace the vector — [Cross-Alignment](@ref) returns
  *aligned* SMLDs. Because the state threads in order, a composite render placed
  **after** an alignment step shows the corrected overlay, while one placed before
  shows the raw registration. Step order is meaningful.

## No provenance is lost

`analyze` returns a `(MultiTargetResult, MultiTargetInfo)` tuple that exposes both
the merged channel vector and each channel's complete single-channel records:

```julia
result.smlds         # Vector{BasicSMLD} (aligned if Cross-Alignment ran)
result[:IgG]         # the channel's full AnalysisResult
info.channels[:IgG]  # the channel's full AnalysisInfo
```

Because every channel keeps its own [`AnalysisResult`](@ref) and
[`AnalysisInfo`](@ref), you can reach any per-channel step's info exactly as in a
single-channel run — see [Data Model & Provenance](@ref).

For the concrete configuration, the `analyze` call shape, the on-disk output
layout, and each cross-channel step in detail, see [Multi-Channel Analysis](@ref).
