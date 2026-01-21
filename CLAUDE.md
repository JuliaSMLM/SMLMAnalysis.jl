# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

SMLMAnalysis.jl is the high-level integration package for the JuliaSMLM ecosystem. It orchestrates all SMLM (Single Molecule Localization Microscopy) analysis packages into unified workflows with checkpointing and provenance tracking.

## Development Commands

```bash
# Run tests
julia --project=. -e 'using Pkg; Pkg.test()'

# Build documentation
julia --project=docs docs/make.jl

# Develop with local JuliaSMLM packages
julia --project=. -e 'using Pkg; Pkg.develop(path="../SMLMData")'

# Run examples (have their own Project.toml)
cd examples && julia --project=. basic_workflow.jl
```

## Architecture

### Step-Based Pipeline

The package uses a step-based architecture where each analysis step is configured via a typed `StepConfig` struct that maps directly to upstream package kwargs:

```julia
# Interactive step-by-step (multi-dataset support)
a = Analysis(images, camera; n_datasets=4, outdir="output/")
run_step!(a, DetectConfig(boxsize=11, min_photons=500))
run_step!(a, FitConfig(psf_model=:variable))
run_step!(a, FilterConfig(min_photons=500))
run_step!(a, DriftCorrectConfig(degree=2))
run_step!(a, RenderConfig(zoom=20))

# Or batch via recipe
recipe = [DetectConfig(...), FitConfig(...), FilterConfig(...)]
result = run_recipe(recipe, images, camera; outdir="output/")

# Reset and try different params (checkpoints auto-created after detect/fit)
reset!(a, 2)  # Go back to after fit
run_step!(a, FilterConfig(min_photons=300))  # Try looser filter

# Resume from disk checkpoint (cross-session)
a = resume_analysis("output/"; step=2)
```

### Module Structure

```
src/
├── SMLMAnalysis.jl      # Main module, re-exports ecosystem types
├── types.jl             # Analysis, StepConfig, StepRecord, Verbosity, DataSource
├── analysis.jl          # run_step!, reset!, checkpoint!, run_recipe, analyze()
├── steps/               # One file per step type
│   ├── common.jl        # Shared helpers (_save_box_overlay, _grid_figure_size)
│   ├── detect.jl        # DetectConfig → SMLMBoxer.getboxes
│   ├── fit.jl           # FitConfig → GaussMLE.fit
│   ├── filter.jl        # FilterConfig → photon/precision/pvalue filtering
│   ├── frameconnect.jl  # FrameConnectConfig → SMLMFrameConnection.frameconnect
│   ├── driftcorrect.jl  # DriftCorrectConfig → SMLMDriftCorrection.driftcorrect
│   ├── isolated.jl      # IsolatedConfig → isolated emitter filtering
│   └── render.jl        # RenderConfig → SMLMRender.render
├── calibration.jl       # Uncertainty calibration from frame connection
└── io/
    ├── smld_io.jl       # HDF5 serialization (save_smld, load_smld)
    ├── smart_h5.jl      # SMART microscope HDF5 import
    └── checkpoint_io.jl # JLD2 checkpoint persistence for cross-session resume
```

### Key Types

- **`Analysis`**: Mutable state container holding DataSource, camera, multi-dataset info (`n_datasets`, `n_frames_per_dataset`), intermediate products (roi_batch, roi_datasets, smld_raw, smld, smld_connected, drift_model), checkpoints, and step history
- **`DataSource`**: Lazy loading wrapper - can hold images directly or a file path for deferred loading
- **`StepConfig`**: Abstract type; each step has a concrete config (DetectConfig, FitConfig, etc.) with kwargs mirroring upstream packages
- **`StepRecord`**: Logged after each step with timing, config, and summary statistics
- **`Verbosity`**: Output detail levels (SILENT=0, PROGRESS=1, STANDARD=2, DETAILED=3, DEBUG=4)

### Data Flow

- All packages use SMLMData.jl types (BasicSMLD, Emitter2DFit, etc.)
- Coordinates are in microns
- Each step's `run_step!` mutates the Analysis in place and optionally writes outputs to `outdir/{step_number}_{step_name}/`
- Checkpoints auto-created after expensive steps (detect, fit) for interactive `reset!`
- Checkpoints can be persisted to disk via `checkpoint=true` constructor arg for cross-session resume

### Multi-Dataset Architecture

For large acquisitions split into multiple datasets (e.g., 4 datasets × 2000 frames = 8000 total frames):

- **Per-dataset processing**: Detection and fitting loop over datasets individually, enabling memory-efficient analysis of arbitrarily large acquisitions
- **Frame numbering**: Frames are per-dataset (1 to `n_frames_per_dataset`), NOT global. This is required for LegendrePoly drift correction which normalizes frames to [-1, 1]
- **Dataset tracking**: `roi_datasets::Vector{Int}` tracks which dataset each ROI belongs to; fit step uses this for correct `emitter.dataset` assignment
- **Global frames for plots only**: Drift correction plots convert to global frame indices for visualization, but internal data stays per-dataset
- **SMLD structure**: `smld.n_frames` = frames per dataset, `smld.n_datasets` = number of datasets

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
- **GaussMLE**: GaussMLEFitter, GaussianXYNB/S/SXSY, LocalizationResult
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
