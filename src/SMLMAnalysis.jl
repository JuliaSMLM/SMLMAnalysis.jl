"""
    SMLMAnalysis

High-level integration package for the JuliaSMLM ecosystem.

Provides a step-based pipeline for SMLM analysis with:
- Typed step configs that mirror upstream package kwargs
- Checkpointing and reset for interactive exploration
- Verbosity levels for controlling output detail
- Recipe-based batch execution for reproducibility

# Quick Start
```julia
using SMLMAnalysis

# Interactive step-by-step
a = Analysis(images, camera; outdir="output/")
run_step!(a, DetectConfig(boxsize=11, min_photons=500))
run_step!(a, FitConfig(psf_model=:variable))
run_step!(a, FilterConfig(min_photons=500))
run_step!(a, DriftCorrectConfig(degree=2))
run_step!(a, RenderConfig(zoom=20))

# Or batch via recipe
recipe = [
    DetectConfig(boxsize=11),
    FitConfig(psf_model=:variable),
    FilterConfig(min_photons=500),
    DriftCorrectConfig(degree=2),
    RenderConfig(zoom=20),
]
result = run_recipe(recipe, images, camera; outdir="output/")

# Reset and try different params
reset!(a, 2)  # Go back to after fit
run_step!(a, FilterConfig(min_photons=300))  # Try looser filter
```

# Re-exported Types
Key types from ecosystem packages are re-exported for convenience:
- SMLMData: AbstractCamera, IdealCamera, SCMOSCamera, BasicSMLD, Emitter types
- GaussMLE: GaussMLEFitter, PSF models, ROIBatch
- SMLMRender: render strategies
"""
module SMLMAnalysis

using Dates
using Statistics
using TOML

# Core dependencies
using SMLMData
using SMLMSim
using SMLMBoxer
using GaussMLE
using SMLMFrameConnection
using SMLMRender
using SMLMDriftCorrection
using MicroscopePSFs
using HDF5
using JLD2
using CairoMakie
using NearestNeighbors

# Re-export from SMLMData
export AbstractCamera, IdealCamera, SCMOSCamera
export AbstractEmitter, Emitter2D, Emitter3D, Emitter2DFit, Emitter3DFit
export BasicSMLD, ROIBatch

# Re-export from SMLMSim
export StaticSMLMParams, DiffusionSMLMParams
export simulate, gen_images, gen_image
export Nmer2D, Nmer3D, GenericFluor

# Re-export from SMLMBoxer
export getboxes

# Re-export from GaussMLE
export GaussMLEFitter
export GaussianXYNB, GaussianXYNBS, GaussianXYNBSXSY, AstigmaticXYZNB
export LocalizationResult

# Re-export from SMLMFrameConnection
export frameconnect

# Re-export from SMLMDriftCorrection
export driftcorrect

# Re-export from SMLMRender
export render, save_image
export HistogramRender, GaussianRender, CircleRender

# ============================================================
# Core types
# ============================================================
include("types.jl")
export Verbosity
export DataSource, get_images
export StepConfig, StepRecord
export AnalysisCheckpoint, Analysis

# ============================================================
# Step configs and run_step! implementations
# ============================================================
include("steps/common.jl")  # Shared helpers for steps

include("steps/detect.jl")
export DetectConfig

include("steps/fit.jl")
export FitConfig

include("steps/filter.jl")
export FilterConfig

include("steps/frameconnect.jl")
export FrameConnectConfig

include("steps/driftcorrect.jl")
export DriftCorrectConfig

include("steps/isolated.jl")
export IsolatedConfig

include("steps/render.jl")
export RenderConfig, RenderSpec, DEFAULT_RENDERS

# ============================================================
# I/O (before analysis.jl - checkpoint_io is used by analysis.jl)
# ============================================================
include("io/smld_io.jl")
export save_smld, load_smld, smld_info

include("io/smart_h5.jl")
export load_smart_h5, load_smart_h5_info, load_smart_h5_frame, smart_h5_to_array

include("io/checkpoint_io.jl")
export resume_analysis

# ============================================================
# Analysis functions
# ============================================================
include("analysis.jl")
export run_step!, reset!, checkpoint!, debug!
export run_recipe, analyze

# ============================================================
# Calibration (used by frameconnect step)
# ============================================================
include("calibration.jl")
export analyze_frameconnect_drift, apply_uncertainty_calibration, recombine_tracks

end # module
