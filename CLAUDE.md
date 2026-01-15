# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

SMLMAnalysis.jl is the high-level integration package for the JuliaSMLM ecosystem. It orchestrates all SMLM (Single Molecule Localization Microscopy) analysis packages into unified workflows with provenance tracking.

## Development Commands

```bash
# Run tests
julia --project=. -e 'using Pkg; Pkg.test()'

# Build documentation
julia --project=docs docs/make.jl

# Develop with local JuliaSMLM packages
julia --project=. -e 'using Pkg; Pkg.develop(path="../SMLMData")'
```

## Architecture

### Module Structure (`src/`)

```
src/
├── SMLMAnalysis.jl      # Main module, re-exports types from ecosystem packages
├── analyze.jl           # Main analyze() pipeline function
├── config.jl            # AnalysisConfig and AnalysisResult types
├── provenance.jl        # SMLMWorkflow and ProcessingStep for tracking
├── filtering.jl         # filter_smld(), filter_isolated()
├── calibration.jl       # Uncertainty calibration from frame connection
├── helpers.jl           # Data format converters between packages
├── io/
│   ├── smld_io.jl       # HDF5 serialization (save_smld, load_smld)
│   └── smart_h5.jl      # SMART microscope HDF5 import
├── figures/
│   └── figures.jl       # All visualization functions
└── stats/
    └── writers.jl       # Markdown stats file writers
```

### Main Workflow

The `analyze()` function runs this pipeline:

```
Images + Camera → Detection (SMLMBoxer) → Fitting (GaussMLE) → Filtering →
Frame Connection → Uncertainty Calibration → Drift Correction → Rendering
```

Each step is configurable via `AnalysisConfig`:
```julia
config = AnalysisConfig(
    fit_model = :variable,      # :fixed, :variable, :anisotropic
    drift = true,
    render = true,
    outdir = "output/"
)
result = analyze(images, camera, config)
```

### Data Flow

- All packages use SMLMData.jl types (BasicSMLD, Emitter2DFit, etc.)
- Coordinates are in microns
- `SMLMWorkflow` tracks each processing step with parameters and timestamps
- HDF5 format (v1.1) stores emitters, camera, drift model, provenance

### Key Files

- **analyze.jl**: The main `analyze()` function (~300 lines) - start here for pipeline logic
- **config.jl**: `AnalysisConfig` struct with all pipeline parameters
- **filtering.jl**: Photon, precision, PSF sigma, and pvalue filtering
- **calibration.jl**: Uncertainty calibration from frame-to-frame analysis
- **figures/figures.jl**: Detection, fitting, drift, and frame connection visualizations
- **stats/writers.jl**: Markdown file generation for diagnostics

### Re-exported Types

From ecosystem packages (available directly after `using SMLMAnalysis`):
- **SMLMData**: `AbstractCamera`, `IdealCamera`, `SCMOSCamera`, `Emitter2D/3D`, `BasicSMLD`
- **SMLMSim**: `StaticSMLMParams`, `simulate`, `gen_images`
- **SMLMBoxer**: `getboxes`
- **GaussMLE**: `GaussMLEFitter`, `fit`, `GaussianXYNB/S/SXSY`, `ROIBatch`
- **SMLMRender**: `render`, `save_image`, `HistogramRender`, `GaussianRender`
- **SMLMFrameConnection**: `frameconnect`
- **SMLMDriftCorrection**: `driftcorrect`
