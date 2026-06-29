"""
    SMLMAnalysis

High-level integration package for the JuliaSMLM ecosystem.

Provides a unified `analyze()` API for SMLM analysis. The config type
determines the operation via multiple dispatch:

# Quick Start
```julia
using SMLMAnalysis

# Full pipeline with AnalysisConfig
config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(
            boxer=BoxerConfig(boxsize=9, psf_sigma=0.130),
            fitter=GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20)),
        FilterConfig(photons=(500.0, Inf)),
        FrameConnectConfig(max_frame_gap=5),
        DriftConfig(degree=2, dataset_mode=:registered),
        RenderConfig(zoom=20, colormap=:inferno),
    ],
    outdir = "output/",
)
(result, info) = analyze(image_stacks, config)

# Individual steps via analyze() dispatch
(smld, info) = analyze(image_stacks, DetectFitConfig(
    camera=cam, boxer=BoxerConfig(boxsize=9, psf_sigma=0.130)))
(smld, info) = analyze(smld, FilterConfig(photons=(500.0, Inf)))
(smld, info) = analyze(smld, FrameConnectConfig(max_frame_gap=5))
(smld, info) = analyze(smld, DriftConfig(degree=2))
(img, info)  = analyze(smld, RenderConfig(zoom=20, colormap=:inferno))
```

# Re-exported Types
Key types from ecosystem packages are re-exported for convenience:
- SMLMData: AbstractCamera, IdealCamera, SCMOSCamera, BasicSMLD, Emitter types
- GaussMLE: GaussMLEConfig, PSF models, ROIBatch
- SMLMFrameConnection: FrameConnectConfig
- SMLMDriftCorrection: DriftConfig
- SMLMRender: render strategies
"""
module SMLMAnalysis

using Dates
using Logging
using Random
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
using SMLMBaGoL
using SMLMClustering
using MicroscopePSFs
using HDF5
using JLD2
using CairoMakie
using NearestNeighbors
using Optim
using Distributions: Poisson, ccdf

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
export getboxes, BoxerConfig

# Re-export from GaussMLE
export GaussMLEConfig
export GaussianXYNB, GaussianXYNBS, GaussianXYNBSXSY, AstigmaticXYZNB
export GaussMLEFitInfo
# Re-export fit - use GaussMLE's fit for fitters
using GaussMLE: fit
export fit

# Re-export from SMLMFrameConnection
export frameconnect
# Re-export FrameConnectConfig (used directly as step config)
const FrameConnectConfig = SMLMFrameConnection.FrameConnectConfig
export FrameConnectConfig
# Re-export CalibrationConfig and CalibrationResult (used via FrameConnectConfig.calibration)
const CalibrationConfig = SMLMFrameConnection.CalibrationConfig
const CalibrationResult = SMLMFrameConnection.CalibrationResult
export CalibrationConfig, CalibrationResult

# Re-export from SMLMDriftCorrection
export driftcorrect
# Re-export alignment API (used by CrossAlignConfig step)
const AlignConfig = SMLMDriftCorrection.AlignConfig
const AlignInfo = SMLMDriftCorrection.AlignInfo
export align_smld, AlignConfig, AlignInfo

# Re-export from SMLMBaGoL
export run_bagol, BaGoLDiagnostics

# Re-export from SMLMRender
export render, save_image
export HistogramRender, GaussianRender, CircleRender, EllipseRender
# Re-export RenderConfig from SMLMRender (used directly as step config)
const RenderConfig = SMLMRender.RenderConfig
export RenderConfig

# Re-export from SMLMClustering
export cluster, cluster_statistics
const AbstractClusterConfig = SMLMClustering.AbstractClusterConfig
const AbstractStatisticsConfig = SMLMClustering.AbstractStatisticsConfig
const ClusterInfo = SMLMClustering.ClusterInfo
const ClusterStatisticsInfo = SMLMClustering.ClusterStatisticsInfo
const DBSCANConfig = SMLMClustering.DBSCANConfig
const HierarchicalConfig = SMLMClustering.HierarchicalConfig
const VoronoiConfig = SMLMClustering.VoronoiConfig
const HopkinsConfig = SMLMClustering.HopkinsConfig
const VoronoiDensityConfig = SMLMClustering.VoronoiDensityConfig
export AbstractClusterConfig, AbstractStatisticsConfig
export ClusterInfo, ClusterStatisticsInfo
export DBSCANConfig, HierarchicalConfig, VoronoiConfig
export HopkinsConfig, VoronoiDensityConfig

# ============================================================
# Core types
# ============================================================
include("types.jl")
export Verbosity, Checkpoint
export DataSource, get_images, n_datasets, n_frames_per_dataset
export AnalysisConfig, AnalysisResult, AnalysisInfo, StepInfo
export DetectFitInfo, FilterInfo, DensityFilterInfo, IntensityFilterInfo, BaGoLInfo
export CompositeRenderInfo, CrossAlignInfo, CrossCorrInfo
export AbstractMultiTargetStep
export MultiTargetConfig, MultiTargetResult, MultiTargetInfo
export crop_camera, crop_images
export step_name

# ============================================================
# Step configs and pure step functions
# ============================================================

# Forward-declare analyze so step files can add dispatch methods
function analyze end

include("steps/common.jl")  # Shared helpers for steps
export step_outdir

include("steps/detectfit.jl")
export DetectFitConfig

include("steps/filter.jl")
export FilterConfig

include("steps/frameconnect.jl")
# FrameConnectConfig is re-exported above (from SMLMFrameConnection)

# CalibrationConfig is re-exported above (from SMLMFrameConnection)

include("steps/driftcorrect.jl")
# DriftConfig is defined as const alias in driftcorrect.jl and exported below
export DriftConfig

include("steps/densityfilter.jl")
export DensityFilterConfig

include("steps/intensityfilter.jl")
export IntensityFilterConfig

include("steps/render.jl")

include("steps/composite_render.jl")
export CompositeRenderConfig

include("steps/cross_align.jl")
export CrossAlignConfig

include("steps/crosscorr.jl")
export CrossCorrConfig

include("steps/bagol.jl")
export BaGoLConfig

include("steps/clustering.jl")
# Clustering config types (DBSCANConfig/HopkinsConfig/…) are re-exported above
# from SMLMClustering; this file only adds analyze() dispatch — no new exports.

# ============================================================
# I/O
# ============================================================
include("io/smld_io.jl")
export save_smld, load_smld, smld_info

include("io/smart_h5.jl")
export load_smart_h5, load_smart_h5_info, load_smart_h5_frame, smart_h5_to_array

include("io/mic_h5.jl")
export load_mic_h5, load_mic_h5_info, load_mic_h5_block
export load_mic_h5_calibration, load_mic_h5_calibration_for_scmos
export build_camera_from_mic_h5

include("io/checkpoint_io.jl")
export save_pipeline_state, load_pipeline_state

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
# Precompilation workload (PrecompileTools)
# ============================================================
# Runs a tiny end-to-end pipeline on synthetic CPU-only data at build time,
# caching the orchestration glue + fit/render specializations into the
# pkgimage. This is the high-leverage spot: `using SMLMAnalysis` is the lab's
# entry point, and the first `analyze()` in a fresh session otherwise pays
# ~1 min of JIT.
#
# Invariants that keep the workload safe to run during precompilation:
#   - backend = :cpu      → no GPU kernels (uncacheable, and no device on CI)
#   - outdir  = nothing   → no disk writes
#   - GaussianXYNBS       → Emitter2DFitSigma, the path the examples exercise
#                           (GaussianXYNB/Emitter2DFitGaussMLE lacks a
#                           `_with_dataset` method — see steps/common.jl)
#   - verbose = SILENT    → no build-time log spam
#   - seeded, dense data  → deterministic and never empty; sparse localization
#                           sets crash downstream reductions over emitter arrays
#
# Disable during active development to skip the workload on every rebuild:
#   using Preferences; set_preferences!(SMLMAnalysis, "precompile_workload" => false; force=true)
using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    # Setup (NOT cached): synthesize a small single-dataset image stack.
    Random.seed!(1)
    cam = IdealCamera(32, 32, 0.1)
    sim = StaticSMLMConfig(density = 5.0, σ_psf = 0.13, nframes = 50, ndatasets = 1)
    (_, si) = simulate(sim;
        pattern  = Nmer2D(n = 8, d = 0.05),
        molecule = GenericFluor(photons = 5.0e4, k_off = 20.0, k_on = 0.04),
        camera   = cam)
    (imgs, _) = gen_images(si.smld_model, MicroscopePSFs.GaussianPSF(0.13);
        dataset = 1, bg = 20.0, poisson_noise = true)

    @compile_workload begin
        # Cached: the detect/fit → filter → frame-connect → render pipeline.
        cfg = AnalysisConfig(
            DetectFitConfig(boxer  = BoxerConfig(boxsize = 7, psf_sigma = 0.13),
                            fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu)),
            FilterConfig(photons = (100.0, Inf)),
            FrameConnectConfig(max_frame_gap = 2),
            RenderConfig(zoom = 10);
            camera  = cam,
            verbose = Verbosity.SILENT,
        )
        analyze([imgs], cfg)
    end
end


end # module
