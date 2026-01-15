"""
    SMLMAnalysis

High-level integration package for the JuliaSMLM ecosystem.

Provides unified workflows that chain together simulation, detection, fitting,
and analysis tools from multiple packages. Tracks analysis provenance for
reproducibility.

# Main entry point
- `analyze(data, camera, config)`: Complete SMLM analysis pipeline

# Core types
- `AnalysisConfig`: Configuration for analysis pipeline
- `AnalysisResult`: Results container with provenance tracking
- `SMLMWorkflow`: Tracks processing steps for reproducibility

# Re-exported types
Key types from SMLMData, SMLMSim, SMLMBoxer, GaussMLE, and SMLMRender are
re-exported for convenience.
"""
module SMLMAnalysis

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
using CairoMakie
using NearestNeighbors
using Statistics
using StatsBase: countmap
using TOML
using LinearAlgebra: det, inv

# Re-export critical types from SMLMData
export AbstractCamera, IdealCamera, SCMOSCamera
export AbstractEmitter, Emitter2D, Emitter3D, Emitter2DFit, Emitter3DFit
export BasicSMLD, SmiteSMLD

# Re-export simulation types and functions
export StaticSMLMParams, DiffusionSMLMParams
export simulate, gen_images, gen_image
export Nmer2D, Nmer3D, GenericFluor

# Re-export detection functions
export getboxes

# Re-export fitting types and functions
export GaussMLEFitter, fit
export GaussianXYNB, GaussianXYNBS, GaussianXYNBSXSY, AstigmaticXYZNB
export ROIBatch, LocalizationResult

# Re-export frame connection
export frameconnect

# Re-export drift correction
export driftcorrect

# Re-export SMLMRender functions and types
export render, save_image
export HistogramRender, GaussianRender, CircleRender

# Include provenance tracking system
include("provenance.jl")
export SMLMWorkflow, ProcessingStep, add_step!, summarize_output

# Include configuration types
include("config.jl")
export AnalysisConfig, AnalysisResult
export save_config, load_config

# Include filtering functions
include("filtering.jl")
export filter_smld, filter_isolated, adaptive_clip_percentile

# Include calibration functions
include("calibration.jl")
export analyze_frameconnect_drift, apply_uncertainty_calibration, recombine_tracks

# Include I/O modules
include("io/smld_io.jl")
export save_smld, load_smld, smld_info

include("io/smart_h5.jl")
export load_smart_h5, load_smart_h5_info, load_smart_h5_frame, smart_h5_to_array

# Include helper functions
include("helpers.jl")
export boxer_to_roi_batch, localization_result_to_smld
export summarize_boxer_result, summarize_fit_result, summarize_smld

# Include output modules as submodules
module Figures
    using CairoMakie
    using Statistics
    using StatsBase: countmap
    include("figures/figures.jl")
end

module Stats
    using Statistics
    include("stats/writers.jl")
end

# Include main analysis pipeline
include("analyze.jl")
export analyze

end
