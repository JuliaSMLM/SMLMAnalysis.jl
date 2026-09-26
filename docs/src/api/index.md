# [API Reference](@id API-Reference)

```@meta
CurrentModule = SMLMAnalysis
```

SMLMAnalysis re-exports a large API from across the JuliaSMLM ecosystem. To keep
each page readable, the reference is split into four parts:

- **Overview** (this page) — the module, the `analyze` verb, and the core pipeline
  types.
- [Step Configs & Info](@ref) — per-step configuration and info structs.
- [Multi-Target & I/O](@ref) — multi-channel types, file import/export, and utilities.
- [Internals](@ref) — non-exported helpers.

```@docs
SMLMAnalysis
```

## Core Functions

```@docs
analyze
stepinfo
stepinfos
```

## Types

```@docs
AnalysisConfig
SMLMAnalysis.AnalysisResult
SMLMAnalysis.AnalysisInfo
SMLMAnalysis.StepInfo
Checkpoint
```

`AnalysisResult`, `AnalysisInfo`, and `StepInfo` are not exported (a user only
*receives* them, never types them to build a pipeline) — reach them as
`SMLMAnalysis.Name` when a type annotation is needed.

### Verbosity Levels

`Verbosity` is a module with integer constants controlling output detail:

| Level | Constant | Output |
|-------|----------|--------|
| 0 | `Verbosity.SILENT` | Errors only |
| 1 | `Verbosity.PROGRESS` | Step names, counts, timing |
| 2 | `Verbosity.STANDARD` | + stats.md, basic figures |
| 3 | `Verbosity.DETAILED` | + diagnostic plots, per-filter breakdowns |
| 4 | `Verbosity.DEBUG` | + MP4 animations, frame-by-frame analysis |
