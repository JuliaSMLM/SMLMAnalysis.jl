"""
    SMLMAnalysis

High-level integration package for the JuliaSMLM ecosystem.

Provides a functional pipeline for SMLM analysis with:
- Pure step functions: detectfit, filter_step, frameconnect_step, etc.
- AnalysisConfig for declarative pipeline description
- analyze() for one-shot execution

# Quick Start
```julia
using SMLMAnalysis

# One-liner with AnalysisConfig
config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(boxsize=9, psf_model=:variable),
        FilterConfig(photons=(500.0, Inf)),
        DriftCorrectConfig(degree=2),
        RenderConfig(zoom=20, colormap=:inferno),
    ],
    outdir = "output/",
)
(result, info) = analyze(image_stacks, config)

# Or use pure step functions directly
(smld, detectfit_info) = detectfit(image_stacks, camera, DetectFitConfig(boxsize=9))
(smld, filter_info) = filter_step(smld, FilterConfig(photons=(500.0, Inf)))
(smld, drift_info) = driftcorrect_step(smld, DriftCorrectConfig(degree=2))
```

# Re-exported Types
Key types from ecosystem packages are re-exported for convenience:
- SMLMData: AbstractCamera, IdealCamera, SCMOSCamera, BasicSMLD, Emitter types
- GaussMLE: GaussMLEConfig, PSF models, ROIBatch
- SMLMRender: render strategies
"""
module SMLMAnalysis

using Dates
using Logging
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
using Optim

# Re-export from SMLMData
export AbstractCamera, IdealCamera, SCMOSCamera
export AbstractEmitter, Emitter2D, Emitter3D, Emitter2DFit, Emitter3DFit
export BasicSMLD, ROIBatch
export AbstractSMLMConfig, AbstractSMLMInfo

# Re-export from SMLMSim
export StaticSMLMConfig, DiffusionSMLMConfig
export simulate, gen_images, gen_image
export Nmer2D, Nmer3D, Line2D, GenericFluor

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

# ============================================================
# Core types
# ============================================================
include("types.jl")
export Verbosity
export DataSource, get_images, n_datasets, n_frames_per_dataset
export AnalysisConfig, AnalysisResult, AnalysisInfo, StepRecord
export MultiTargetConfig, MultiTargetResult, MultiTargetInfo
export crop_camera, crop_images
export step_name

# ============================================================
# Step configs and pure step functions
# ============================================================
include("steps/common.jl")  # Shared helpers for steps
export step_outdir

include("steps/detectfit.jl")
export DetectFitConfig, detectfit

include("steps/filter.jl")
export FilterConfig, filter_step

include("steps/frameconnect.jl")
export FrameConnectConfig, frameconnect_step

include("steps/driftcorrect.jl")
export DriftCorrectConfig, driftcorrect_step

include("steps/densityfilter.jl")
export DensityFilterConfig, densityfilter_step

include("steps/render.jl")
export render_step

# ============================================================
# I/O
# ============================================================
include("io/smld_io.jl")
export save_smld, load_smld, smld_info

include("io/smart_h5.jl")
export load_smart_h5, load_smart_h5_info, load_smart_h5_frame, smart_h5_to_array

include("io/lidkelab_h5.jl")
export load_lidkelab_h5, load_lidkelab_h5_info, load_lidkelab_h5_block
export load_lidkelab_h5_calibration, load_lidkelab_h5_calibration_for_scmos

include("io/checkpoint_io.jl")

# ============================================================
# Analysis orchestrator
# ============================================================
include("analysis.jl")
export analyze

# ============================================================
# Multi-target orchestration
# ============================================================
include("multitarget.jl")

# ============================================================
# Calibration (used by frameconnect step)
# ============================================================
include("calibration.jl")
export analyze_frameconnect_drift, apply_uncertainty_calibration, recombine_tracks

end # module
