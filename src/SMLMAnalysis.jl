"""
    SMLMAnalysis

High-level integration package for the JuliaSMLM ecosystem.

Provides unified workflows that chain together simulation, detection, fitting,
and analysis tools from multiple packages. Tracks analysis provenance for
reproducibility.

# Main workflows
- `simulate_detect_fit_workflow`: Complete pipeline from simulation to fitted localizations
- `standard_localization_workflow`: Detection and fitting for experimental data

# Core types
- `SMLMWorkflow`: Tracks processing steps for reproducibility
- Helper functions for data format conversions

# Re-exported types
Key types from SMLMData, SMLMSim, SMLMBoxer, and GaussMLE are re-exported
for convenience.
"""
module SMLMAnalysis

# Core dependencies
using SMLMData
using SMLMSim
using SMLMBoxer
using GaussMLE
# using SMLMFrameConnection  # temporarily disabled - StatsBase conflict
using SMLMRender
using HDF5

# Re-export critical types from SMLMData
export AbstractCamera, IdealCamera, SCMOSCamera
export AbstractEmitter, Emitter2D, Emitter3D, Emitter2DFit, Emitter3DFit
export BasicSMLD, SmiteSMLD

# Re-export simulation types and functions
export StaticSMLMParams, DiffusionSMLMParams
export simulate, gen_images, gen_image
export Nmer2D, Nmer3D, GenericFluor  # Pattern and molecule types

# Re-export detection functions
export getboxes

# Re-export fitting types and functions
export GaussMLEFitter, fit
export GaussianXYNB, GaussianXYNBS, GaussianXYNBSXSY, AstigmaticXYZNB
export ROIBatch, LocalizationResult

# Re-export frame connection function
# export frameconnect  # temporarily disabled

# Include workflow tracking system
include("workflow.jl")
export SMLMWorkflow, ProcessingStep, add_step!, summarize_output

# Include helper functions
include("helpers.jl")
export boxer_to_roi_batch, localization_result_to_smld
export summarize_boxer_result, summarize_fit_result, summarize_smld

# Include data import functions
include("import_smart_h5.jl")
export load_smart_h5, load_smart_h5_info, load_smart_h5_frame, smart_h5_to_array

# Include SMLD I/O functions
include("io.jl")
export save_smld, load_smld, smld_info

# Include high-level workflows
include("workflows.jl")
export simulate_detect_fit_workflow, standard_localization_workflow
export AnalysisConfig, AnalysisResult, analyze
export save_config, load_config

# Re-export SMLMRender functions and types
export render, save_image
export HistogramRender, GaussianRender, CircleRender  # Rendering strategies

end
