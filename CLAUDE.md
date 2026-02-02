# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

SMLMAnalysis.jl is the high-level integration package for the JuliaSMLM ecosystem. It orchestrates all SMLM (Single Molecule Localization Microscopy) analysis packages into unified workflows with checkpointing and provenance tracking.

## JuliaSMLM Ecosystem

SMLMAnalysis depends on and orchestrates these JuliaSMLM packages. Each package has its own agent that can be contacted for coordination.

### Package Hierarchy

```
SMLMData (core types - no deps)
    ↓
┌───┴───┬───────────┬──────────────┬─────────────────┬─────────┐
│       │           │              │                 │         │
SMLMBoxer  GaussMLE  SMLMRender  SMLMFrameConnection  SMLMSim  SMLMBaGoL
│           │                                         │
└─────┬─────┘                                    MicroscopePSFs
      │
SMLMDriftCorrection (also depends on SMLMSim)
      │
SMLMAnalysis (integrates all)
```

### Package Details

| Package | Path | Agent | Expertise |
|---------|------|-------|-----------|
| SMLMData | `../SMLMData` | @data | Core types: Emitter, Camera, BasicSMLD, ROIBatch |
| SMLMBoxer | `../SMLMBoxer` | @boxer | ROI detection from images |
| GaussMLE | `../GaussMLE` | @gaussmle | GPU-accelerated MLE fitting |
| SMLMFrameConnection | `../SMLMFrameConnection` | @frameconnect | Linking localizations across frames |
| SMLMDriftCorrection | `../SMLMDriftCorrection` | @drift | Sample drift correction |
| SMLMRender | `../SMLMRender` | @render | Super-resolution image rendering |
| SMLMSim | `../SMLMSim` | @sim | SMLM data simulation, fluorophore kinetics |
| MicroscopePSFs | `../MicroscopePSFs` | @psf | PSF models (Gaussian, Airy, etc.) |
| SMLMBaGoL | `../SMLMBaGoL` | @bagol | Bayesian grouping of localizations |

### Agent Communication

To coordinate changes across packages:

```bash
# Send message to specific agent
agent send @data "Question about Emitter2DFit fields"

# Broadcast to all JuliaSMLM agents
agent send --broadcast "Need SMLMData 0.6 compat update"

# Check who's online
agent list
```

Common coordination scenarios:
- **Type changes**: @data announces, all packages update compat
- **API design**: Convene multi-party discussion via @data
- **Breaking changes**: Use `[breaking]` tag, each package assesses impact

## Development Commands

```bash
# Run tests
julia --project=. -e 'using Pkg; Pkg.test()'

# Build documentation
julia --project=docs docs/make.jl

# Develop with local JuliaSMLM packages
julia --project=. -e 'using Pkg; Pkg.develop(path="../SMLMData")'

# Run examples (have their own Project.toml)
cd examples && julia --project=. stepwise_example.jl
```

## Architecture

### Step-Based Pipeline

The package uses a step-based architecture where each analysis step is configured via a typed `StepConfig` struct that maps directly to upstream package kwargs:

```julia
# Interactive step-by-step (multi-dataset support)
a = Analysis(images, camera; n_datasets=4, outdir="output/")
run_step!(a, DetectFitConfig(boxsize=9, psf_model=:variable))
run_step!(a, FilterConfig(photons=(500.0, Inf)))
run_step!(a, FrameConnectConfig(maxframegap=5))
run_step!(a, DriftCorrectConfig(degree=2))
run_step!(a, RenderConfig(zoom=20))

# One-liner with defaults
result = analyze(images, camera; outdir="output/", n_datasets=4)

# Reset and try different params (checkpoints auto-created after detectfit)
reset!(a, 1)  # Go back to after detectfit
run_step!(a, FilterConfig(photons=(300.0, Inf)))  # Try looser filter

# Resume from disk checkpoint (cross-session)
a = resume_analysis("output/"; images=images)
```

### Module Structure

```
src/
├── SMLMAnalysis.jl      # Main module, re-exports ecosystem types
├── types.jl             # Analysis, StepConfig, StepRecord, Verbosity, DataSource
├── analysis.jl          # run_step!, reset!, checkpoint!, debug!, analyze()
├── steps/               # One file per step type
│   ├── common.jl        # Shared helpers (_save_box_overlay, _calculate_mode)
│   ├── detectfit.jl     # DetectFitConfig → combined SMLMBoxer.getboxes + GaussMLE.fit
│   ├── filter.jl        # FilterConfig → photon/precision/pvalue/psf filtering
│   ├── frameconnect.jl  # FrameConnectConfig → SMLMFrameConnection.frameconnect
│   ├── driftcorrect.jl  # DriftCorrectConfig → SMLMDriftCorrection.driftcorrect
│   ├── isolated.jl      # IsolatedConfig → isolated emitter filtering
│   └── render.jl        # RenderConfig → SMLMRender.render
├── calibration.jl       # Uncertainty calibration from frame connection
└── io/
    ├── smld_io.jl       # HDF5 serialization (save_smld, load_smld)
    ├── smart_h5.jl      # SMART microscope HDF5 import
    ├── lidkelab_h5.jl   # LidkeLab/MIC H5 format import (block-based loading)
    └── checkpoint_io.jl # JLD2 checkpoint persistence for cross-session resume
```

### Key Types

- **`Analysis`**: Mutable state container holding DataSource, camera, multi-dataset info (`n_datasets`, `n_frames_per_dataset`), intermediate products (roi_batch, roi_datasets, smld_raw, smld, smld_connected, drift_model), checkpoints, and step history
- **`DataSource`**: Lazy loading wrapper - can hold images directly or a file path for deferred loading
- **`StepConfig`**: Abstract type; each step has a concrete config with kwargs mirroring upstream packages
- **`StepRecord`**: Logged after each step with timing, config, summary statistics, and upstream info struct
- **`AnalysisInfo`**: Aggregated metadata from all steps, containing per-step info structs from upstream packages
- **`Verbosity`**: Output detail levels (SILENT=0, PROGRESS=1, STANDARD=2, DETAILED=3, DEBUG=4)

### Tuple-Pattern API

SMLMAnalysis follows the JuliaSMLM tuple-pattern where functions return `(result, Info)` tuples:

```julia
# analyze() returns (Analysis, AnalysisInfo)
(result, info) = analyze(images, camera; outdir="output/")
result.smld               # Final SMLD
info.elapsed_ns           # Total time in nanoseconds
info.steps[:detectfit]    # Per-step info from upstream packages
info.steps[:driftcorrect] # DriftInfo from SMLMDriftCorrection

# Interactive usage - extract info after running steps
a = Analysis(images, camera)
run_step!(a, DetectFitConfig())
run_step!(a, FilterConfig(photons=(500.0, Inf)))
info = get_analysis_info(a)  # Builds AnalysisInfo from step records
```

Each step internally handles tuple returns from upstream packages:
- `getboxes()` → `(ROIBatch, BoxesInfo)`
- `fit()` → `(BasicSMLD, FitInfo)`
- `frameconnect()` → `(combined, ConnectInfo)` with `.connected` in info
- `driftcorrect()` → `(corrected_smld, DriftInfo)` with `.model` in info
- `render()` → `(image, RenderInfo)`

Step records store upstream info in `step.info` field for later access.

### Data Flow

- All packages use SMLMData.jl types (BasicSMLD, Emitter2DFit, etc.)
- Coordinates are in microns
- Each step's `run_step!` mutates the Analysis in place and optionally writes outputs to `outdir/{step_number}_{step_name}/`
- Checkpoints auto-created after expensive steps (detectfit) for interactive `reset!`
- Checkpoints can be persisted to disk via `checkpoint=true` constructor arg for cross-session resume

### Multi-Dataset Architecture

For large acquisitions split into multiple datasets (e.g., 4 datasets × 2000 frames = 8000 total frames):

- **Per-dataset processing**: Detection and fitting loop over datasets individually, enabling memory-efficient analysis of arbitrarily large acquisitions
- **Frame numbering**: Frames are per-dataset (1 to `n_frames_per_dataset`), NOT global. This is required for LegendrePoly drift correction which normalizes frames to [-1, 1]
- **Dataset tracking**: `emitter.dataset` field tracks which dataset each emitter belongs to
- **Global frames for plots only**: Drift correction plots convert to global frame indices for visualization, but internal data stays per-dataset
- **SMLD structure**: `smld.n_frames` = frames per dataset, `smld.n_datasets` = number of datasets

### DetectFitConfig Data Sources

The combined detection+fitting step supports three data source modes:

1. **In-memory images**: Pass images to `Analysis(images, camera; n_datasets=N)` constructor
2. **Single file, multiple datasets**: `DetectFitConfig(path="data.h5", n_datasets=4)` splits frames evenly
3. **Multiple files**: `DetectFitConfig(paths=["d1.h5", "d2.h5"])` loads one file per dataset

H5 formats auto-detected: `:smart` (SMART microscope), `:mic` (LidkeLab MIC format with block-based loading)

### Adding a New Step

1. Create `src/steps/yourstep.jl` with:
   - `@kwdef struct YourStepConfig <: StepConfig` with `name::String` and `verbose::Int` fields
   - `run_step!(a::Analysis, cfg::YourStepConfig)` that calls upstream package and updates `a`
   - Optional `_save_step_outputs!` for verbosity-gated figures/stats
2. Include it in `SMLMAnalysis.jl` and export the config
3. Each step must:
   - Increment `a.step_counter` first
   - Call `_record!(a, cfg, timing, summary)` after execution
   - Call `_checkpoint!(a)` after expensive operations

### Re-exported Types

From ecosystem packages (available after `using SMLMAnalysis`):
- **SMLMData**: AbstractCamera, IdealCamera, SCMOSCamera, Emitter2D/3D, BasicSMLD, ROIBatch
- **GaussMLE**: GaussMLEFitter, GaussianXYNB/S/SXSY, LocalizationResult, fit
- **SMLMBoxer**: getboxes
- **SMLMFrameConnection**: frameconnect
- **SMLMDriftCorrection**: driftcorrect
- **SMLMRender**: render, save_image, HistogramRender, GaussianRender, CircleRender
- **SMLMSim**: StaticSMLMParams, DiffusionSMLMParams, simulate, gen_images, Nmer2D, Nmer3D, GenericFluor

## Calibration Module

The calibration.jl module provides uncertainty calibration from frame connection analysis:
- `analyze_frameconnect_drift(smld_connected)`: Analyzes frame-to-frame drift, returns χ² statistics and calibration model
- `apply_uncertainty_calibration(smld, σ_motion, k_scale)`: Applies calibration: σ²_corrected = σ²_motion + k² × σ²_CRLB
- `recombine_tracks(smld_connected)`: Weighted average of track localizations with calibrated uncertainties
