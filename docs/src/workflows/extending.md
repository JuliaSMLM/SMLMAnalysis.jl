```@meta
CurrentModule = SMLMAnalysis
```

# Extending the Pipeline

Because the pipeline is routed entirely by [Julia's multiple dispatch](@ref
"The Pipeline Model"), **you can add a new step without modifying SMLMAnalysis at
all.** Define a config type and an `analyze(smld, cfg)` method — in your own
package, in a script, or even at the REPL — and your step composes into
`AnalysisConfig.steps` like any built-in one. There is no step registry to
register with and no orchestrator to edit; the fold simply calls `analyze(state,
cfg)` and dispatch finds your method.

This is the same mechanism the built-in steps use. The only requirements are:

1. Your config subtypes `SMLMData.AbstractSMLMConfig` (so it is a valid element
   of the `steps` vector).
2. You add a method `analyze(smld::BasicSMLD, cfg::YourConfig; kwargs...)` that
   returns a `(result, StepInfo)` tuple.
3. The method accepts (and may ignore) the orchestrator's keyword arguments via a
   trailing `kwargs...`.

## Minimal external step

The smallest useful step is a few lines. Put this in your own module or script —
nothing here lives inside SMLMAnalysis:

```julia
using SMLMAnalysis
using SMLMAnalysis: analyze, StepInfo, Verbosity   # extend `analyze`, build a StepInfo
using SMLMData

# 1. A config type — parameters only, no logic.
struct SpatialFilterConfig <: SMLMData.AbstractSMLMConfig
    x_range::Tuple{Float64,Float64}   # microns
    y_range::Tuple{Float64,Float64}
end

# 2. An analyze() method. Dispatch on (BasicSMLD, SpatialFilterConfig) routes
#    the pipeline here automatically.
function SMLMAnalysis.analyze(smld::BasicSMLD, cfg::SpatialFilterConfig;
                              step_number::Int = 0, verbose::Int = Verbosity.STANDARD,
                              kwargs...)
    keep = [cfg.x_range[1] <= e.x <= cfg.x_range[2] &&
            cfg.y_range[1] <= e.y <= cfg.y_range[2] for e in smld.emitters]
    out = BasicSMLD(smld.emitters[keep], smld.camera,
                    smld.n_frames, smld.n_datasets, smld.metadata)

    summary = Dict{Symbol,Any}(:n_before => length(smld.emitters),
                               :n_after  => length(out.emitters))
    (out, StepInfo(step_number, cfg, 0.0, summary))
end
```

That is enough to use it both standalone and in a pipeline:

```julia
# Standalone
(cropped, info) = analyze(smld, SpatialFilterConfig((1.0, 5.0), (1.0, 5.0)))

# As a pipeline step — it slots in next to the built-ins
config = AnalysisConfig(camera = cam, steps = [
    DetectFitConfig(boxer = BoxerConfig(boxsize = 9, psf_sigma = 0.130)),
    SpatialFilterConfig((1.0, 5.0), (1.0, 5.0)),   # ← your step
    RenderConfig(zoom = 20),
])
(result, info) = analyze(image_stacks, config)
```

If you list `SpatialFilterConfig` before any `DetectFitConfig`, you get a clear
`MethodError` (no `analyze(::Vector{…}, ::SpatialFilterConfig)`), not a silent
wrong result — see [error clarity](@ref "The Pipeline Model").

## The recommended structure

For anything beyond a few lines, follow the same two-layer shape the built-in
steps use: an **internal work function** that does the computation and returns
`(result, <YourInfo>)`, and a thin **`analyze()` wrapper** that times it and
builds the [`StepInfo`](@ref). This keeps the work reusable and testable
independently of the pipeline.

### 1. Config and info types

```julia
@kwdef struct SpatialFilterConfig <: SMLMData.AbstractSMLMConfig
    x_range::Tuple{Float64,Float64}
    y_range::Tuple{Float64,Float64}
end

struct SpatialFilterInfo <: SMLMData.AbstractSMLMInfo
    n_before::Int
    n_after::Int
end
```

Subtyping `AbstractSMLMInfo` lets your info struct be stored on `StepInfo.info`
and surfaced in [`AnalysisInfo`](@ref) alongside the upstream info structs.

### 2. Internal work function

```julia
function spatialfilter_step(smld::BasicSMLD, cfg::SpatialFilterConfig;
                            outdir = nothing, step_number::Int = 0,
                            verbose::Int = Verbosity.STANDARD)
    v = verbose
    dir = SMLMAnalysis.step_outdir(outdir, step_number, cfg)
    v >= Verbosity.PROGRESS && @info "[$step_number] $(SMLMAnalysis.step_name(cfg))"

    n_before = length(smld.emitters)
    keep = [cfg.x_range[1] <= e.x <= cfg.x_range[2] &&
            cfg.y_range[1] <= e.y <= cfg.y_range[2] for e in smld.emitters]
    out = BasicSMLD(smld.emitters[keep], smld.camera,
                    smld.n_frames, smld.n_datasets, smld.metadata)

    if dir !== nothing
        mkpath(dir)
        SMLMAnalysis._save_config!(dir, cfg)     # writes config.toml for provenance
    end
    v >= Verbosity.PROGRESS && @info "  → $(length(out.emitters)) / $n_before"
    (out, SpatialFilterInfo(n_before, length(out.emitters)))
end
```

### 3. Summary dispatch

`_step_summary` turns your info struct into the `Dict` shown in the pipeline's
summary tables and `stats.md`:

```julia
SMLMAnalysis._step_summary(info::SpatialFilterInfo) = Dict{Symbol,Any}(
    :n_before   => info.n_before,
    :n_after    => info.n_after,
    :acceptance => round(info.n_after / max(1, info.n_before), digits = 3),
)
```

### 4. `analyze()` wrapper

```julia
function SMLMAnalysis.analyze(smld::BasicSMLD, cfg::SpatialFilterConfig;
                              outdir = nothing, step_number::Int = 0,
                              verbose::Int = Verbosity.STANDARD, kwargs...)
    t = @elapsed (out, sf_info) = spatialfilter_step(smld, cfg;
        outdir = outdir, step_number = step_number, verbose = verbose)
    (out, StepInfo(step_number, cfg, t, SMLMAnalysis._step_summary(sf_info); info = sf_info))
end
```

The trailing `kwargs...` is required: the orchestrator passes context keywords
(`outdir`, `step_number`, `verbose`, `checkpoint`, and occasionally others) to
every step, and yours must accept the ones it does not use.

## Wrapping an upstream package

A step that wraps an existing `(result, info)`-returning function is even
thinner — the work function *is* the upstream call. This is exactly how the
[Clustering](@ref clustering-step) step wraps `SMLMClustering.cluster`, [Drift Correction](@ref)
wraps `SMLMDriftCorrection.driftcorrect`, and [Bayesian Grouping](@ref "Bayesian
Grouping (BaGoL)") wraps `SMLMBaGoL.run_bagol`. Follow the "upstream owns the
config" idiom: alias the upstream config (`const YourConfig = ThatPkg.Config`)
and dispatch your `analyze()` method on it directly rather than re-declaring its
fields.

## Contributing a step upstream

If your step is generally useful and you want it in SMLMAnalysis itself, the
in-repo recipe is the same five pieces, plus: put the work function in
`src/steps/yourstep.jl`, `include` and `export` it from `src/SMLMAnalysis.jl`,
add a docstring and an [API Reference](@ref) entry, and add a step page under
[Pipeline Steps](@ref "Pipeline Steps: Overview"). See `CONTRIBUTING` for the
full checklist.
