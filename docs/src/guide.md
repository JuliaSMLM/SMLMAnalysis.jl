# [User Guide](@id Guide)

This page is an orientation map for the rest of the documentation. Each topic now
has its own dedicated page under **Concepts**, **Workflows**, and **Pipeline
Steps** — use the task-oriented links below to jump straight to what you need.

New to the package? Start with the [Tutorial](@ref Tutorial), which runs a complete
pipeline on simulated data and shows the output of every step.

## Find what you need

| I want to… | Go to |
|------------|-------|
| Understand how `analyze()` routes steps by dispatch | [The Pipeline Model](@ref) |
| See where SMLMAnalysis sits in the JuliaSMLM ecosystem | [The JuliaSMLM Ecosystem](@ref) |
| Know what an SMLD stores and how provenance is tracked | [Data Model & Provenance](@ref) |
| Install the package and its unregistered dependencies | [Installation & Setup](@ref) |
| Run a full pipeline from an `AnalysisConfig` | [Running a Pipeline](@ref) |
| Set the output detail level (verbosity) | [Running a Pipeline](@ref) |
| Process multiple files, or chunk one long movie | [Multi-Dataset Acquisitions](@ref) |
| Analyze multiple color channels and overlay them | [Multi-Channel Analysis](@ref) |
| Save intermediate state and resume a session later | [I/O & Resume](@ref) |
| Import SMART / MIC microscope `.h5` data | [I/O & Resume](@ref) |
| Add a custom step to the pipeline | [Extending the Pipeline](@ref) |
| Diagnose a failure or an unexpected result | [Troubleshooting](@ref) |

## Configure a step

Every pipeline step has a reference page documenting its config fields, the
algorithm it dispatches to, and the diagnostic outputs it writes:

[Detection & Fitting](@ref) · [Quality Filter](@ref) · [Intensity Filter](@ref) ·
[Frame Connection](@ref) · [Drift Correction](@ref) · [Density Filter](@ref) ·
[Rendering](@ref) · [Bayesian Grouping (BaGoL)](@ref) · [Clustering](@ref clustering-step)

Multi-channel steps: [Composite Render](@ref) · [Cross-Alignment](@ref) ·
[Cross-Correlation](@ref)

!!! note "Uncertainty calibration is not a separate step"
    CRLB-vs-observed uncertainty calibration is configured *inside* frame
    connection via `FrameConnectConfig(calibration = CalibrationConfig(...))`, so
    it is documented on the [Frame Connection](@ref) page rather than as its own
    step.

See [Pipeline Steps: Overview](@ref) for the full list and the ordering rules
(`DetectFitConfig` must come first; most other steps are freely orderable,
repeatable, or optional), and the [API Reference](@ref API-Reference) for every
exported config, info, and function.
