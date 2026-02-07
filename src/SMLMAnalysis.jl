"""
    SMLMAnalysis

High-level integration package for the JuliaSMLM ecosystem.

Provides a step-based pipeline for SMLM analysis with:
- Typed step configs that mirror upstream package kwargs
- Checkpointing and reset for interactive exploration
- Verbosity levels for controlling output detail

# Quick Start
```julia
using SMLMAnalysis

# One-liner with defaults
result = analyze(images, camera; outdir="output/", n_datasets=4)

# Or interactive step-by-step
a = Analysis(images, camera; outdir="output/", n_datasets=4)
run_step!(a, DetectFitConfig(boxsize=9, psf_model=:variable))
run_step!(a, FilterConfig(photons=(500.0, Inf)))
run_step!(a, FrameConnectConfig(maxframegap=5))
run_step!(a, DriftCorrectConfig(degree=2))
run_step!(a, RenderConfig(zoom=20, colormap=:inferno))

# Reset and try different params
reset!(a, 1)  # Go back to after detectfit
run_step!(a, FilterConfig(photons=(300.0, Inf)))  # Try looser filter
```

# Re-exported Types
Key types from ecosystem packages are re-exported for convenience:
- SMLMData: AbstractCamera, IdealCamera, SCMOSCamera, BasicSMLD, Emitter types
- GaussMLE: GaussMLEConfig, PSF models, ROIBatch
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
using SMLMBaGoL
using HDF5
using JLD2
using CairoMakie
using NearestNeighbors
using Optim

# Re-export from SMLMData
export AbstractCamera, IdealCamera, SCMOSCamera
export AbstractEmitter, Emitter2D, Emitter3D, Emitter2DFit, Emitter3DFit
export BasicSMLD, ROIBatch
export AbstractSMLMConfig, AbstractSMLMInfo

# Re-export from SMLMSim
export StaticSMLMConfig, DiffusionSMLMConfig
export simulate, gen_images, gen_image
export Nmer2D, Nmer3D, GenericFluor

# Re-export from SMLMBoxer
export getboxes

# Re-export from GaussMLE
export GaussMLEConfig
export GaussianXYNB, GaussianXYNBS, GaussianXYNBSXSY, AstigmaticXYZNB
export GaussMLEFitInfo
# Re-export fit - use GaussMLE's fit for fitters
using GaussMLE: fit
export fit

# Re-export from SMLMFrameConnection
export frameconnect

# Re-export from SMLMDriftCorrection
export driftcorrect

# Re-export from SMLMRender
export render, save_image
export HistogramRender, GaussianRender, CircleRender, EllipseRender
# Re-export RenderConfig from SMLMRender (used directly as step config)
const RenderConfig = SMLMRender.RenderConfig
export RenderConfig

# Re-export from SMLMBaGoL
export MAPNResult, run_bagol, estimate_mapn

# ============================================================
# Core types
# ============================================================
include("types.jl")
export Verbosity
export DataSource, get_images
export AnalysisConfig, StepRecord, AnalysisInfo
export AnalysisCheckpoint, Analysis
export crop_camera, crop_images
export step_name

# ============================================================
# Step configs and run_step! implementations
# ============================================================
include("steps/common.jl")  # Shared helpers for steps

include("steps/detectfit.jl")
export DetectFitConfig

include("steps/filter.jl")
export FilterConfig

include("steps/frameconnect.jl")
export FrameConnectConfig

include("steps/driftcorrect.jl")
export DriftCorrectConfig

include("steps/densityfilter.jl")
export DensityFilterConfig

include("steps/render.jl")

include("steps/bagol.jl")
export BaGoLConfig

# ============================================================
# I/O (before analysis.jl - checkpoint_io is used by analysis.jl)
# ============================================================
include("io/smld_io.jl")
export save_smld, load_smld, smld_info

include("io/smart_h5.jl")
export load_smart_h5, load_smart_h5_info, load_smart_h5_frame, smart_h5_to_array

include("io/lidkelab_h5.jl")
export load_lidkelab_h5, load_lidkelab_h5_info, load_lidkelab_h5_block
export load_lidkelab_h5_calibration, load_lidkelab_h5_calibration_for_scmos

include("io/checkpoint_io.jl")
export resume_analysis

# ============================================================
# Analysis functions
# ============================================================
include("analysis.jl")
export run_step!, reset!, checkpoint!, debug!
export analyze, get_analysis_info, get_config

# ============================================================
# Calibration (used by frameconnect step)
# ============================================================
include("calibration.jl")
export analyze_frameconnect_drift, apply_uncertainty_calibration, recombine_tracks

end # module
