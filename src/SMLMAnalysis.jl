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
# RenderConfig is a pass-through step: it returns (smld, StepInfo) and writes the
# image to outdir. For the image in-memory, call SMLMRender.render(smld, cfg) directly.
(smld, info) = analyze(smld, RenderConfig(zoom=20, colormap=:inferno))
```

# Re-exported Types
Only the names a user must type to build/run a pipeline are exported; see
`api_overview.md` for the exact list. Highlights:
- SMLMData: IdealCamera, SCMOSCamera (cameras)
- GaussMLE: GaussMLEConfig, PSF models
- SMLMFrameConnection: FrameConnectConfig, CalibrationConfig
- SMLMDriftCorrection: DriftConfig, AlignConfig
- SMLMBaGoL: BaGoLConfig
- SMLMRender: RenderConfig, render strategies
- SMLMClustering: clustering/spatial-statistics/edge-classify config types

Result/info structs (e.g. `AnalysisResult`, `DetectFitInfo`) and upstream verbs
(`simulate`, `render`, `cluster`, `fit`, …) are not exported — reach them via
`SMLMAnalysis.Name` or their owning package.
"""
module SMLMAnalysis

using Dates: Dates, DateTime, now
using Logging: Logging
using Random: Random
using Statistics: Statistics, mean, median, quantile, std, var
using TOML: TOML

# Core dependencies. Every name is imported explicitly (ExplicitImports checks this in
# test/qa/qa.jl); qualified access (SMLMDriftCorrection.driftcorrect) needs only the module.
using SMLMData: SMLMData, AbstractEmitter, BasicSMLD, Emitter2DFit, Emitter3DFit,
    IdealCamera, ROIBatch, SCMOSCamera
using SMLMSim: SMLMSim, GenericFluor, Nmer2D, StaticSMLMConfig, gen_images, simulate
using SMLMBoxer: SMLMBoxer, BoxerConfig
using GaussMLE: GaussMLE, AstigmaticXYZNB, GaussMLEConfig, GaussianXYNB, GaussianXYNBS,
    GaussianXYNBSXSY
using SMLMFrameConnection: SMLMFrameConnection
using SMLMRender: SMLMRender, CircleRender, EllipseRender, GaussianRender, HistogramRender
using SMLMDriftCorrection: SMLMDriftCorrection
using SMLMBaGoL: SMLMBaGoL
using SMLMClustering: SMLMClustering
using MicroscopePSFs: MicroscopePSFs
using HDF5: HDF5, create_group, h5open
using JLD2: JLD2
# Makie's plotting names (and FileIO's save) come through CairoMakie; Makie itself is not a
# direct dependency.
using CairoMakie: CairoMakie, Axis, Colorbar, DataAspect, Figure, GridLayout, Label, Legend,
    Point2f, axislegend, band!, barplot!, heatmap!, hidedecorations!, hidespines!,
    hideydecorations!, hist!, hlines!, lines!, poly!, save, scatter!, text, text!, vlines!,
    vspan!, xlims!, ylims!
using NearestNeighbors: NearestNeighbors, KDTree, inrange
using Optim: Optim, NelderMead, optimize
using Distributions: Poisson, ccdf, Gamma, pdf

# Upstream names reachable as SMLMAnalysis.Name without being exported (api_overview.md,
# "Non-exported but public"; test/exports.jl checks each one resolves). Nothing in src uses
# them, so test/qa/qa.jl exempts them from ExplicitImports' stale-import check.
# The three SMLMData abstract types are extension hooks.
using SMLMData: AbstractCamera, AbstractSMLMConfig, AbstractSMLMInfo
using SMLMSim: Line2D
using GaussMLE: fit
using SMLMFrameConnection: CalibrationResult, frameconnect
using SMLMDriftCorrection: AlignInfo, align_smld, driftcorrect
using SMLMBaGoL: BaGoLDiagnostics, run_bagol
using SMLMRender: render
using SMLMClustering: AbstractClusterConfig, AbstractEdgeClassifyConfig,
    AbstractStatisticsConfig, CellPolygon, ClusterInfo, ClusterStatisticsInfo,
    EdgeClassifyInfo, MultiCellMask, cluster, cluster_statistics, in_cell, interior_fraction,
    interior_mask

# Re-export from SMLMData (cameras only — a user must type these to build a
# pipeline; Emitter2DFit/Emitter3DFit/BasicSMLD/ROIBatch are receive-only and
# reached via SMLMData.Name; AbstractCamera/AbstractSMLMConfig/AbstractSMLMInfo
# are extension hooks, not exported — see workflows/extending.md)
export IdealCamera, SCMOSCamera

# Re-export from SMLMBoxer
export BoxerConfig

# Re-export from GaussMLE
export GaussMLEConfig
export GaussianXYNB, GaussianXYNBS, GaussianXYNBSXSY, AstigmaticXYZNB

# Re-export from SMLMFrameConnection
# Re-export FrameConnectConfig (used directly as step config)
const FrameConnectConfig = SMLMFrameConnection.FrameConnectConfig
export FrameConnectConfig
# Re-export CalibrationConfig (used via FrameConnectConfig.calibration)
const CalibrationConfig = SMLMFrameConnection.CalibrationConfig
export CalibrationConfig

# Re-export from SMLMDriftCorrection
# Re-export DriftConfig (used directly as a pipeline step, like RenderConfig)
const DriftConfig = SMLMDriftCorrection.DriftConfig
export DriftConfig
# Re-export AlignConfig (used by CrossAlignConfig step)
const AlignConfig = SMLMDriftCorrection.AlignConfig
export AlignConfig

# Re-export from SMLMBaGoL
# Re-export BaGoLConfig (upstream owns the config, used directly as a pipeline step)
const BaGoLConfig = SMLMBaGoL.BaGoLConfig
export BaGoLConfig

# Re-export from SMLMRender
export HistogramRender, GaussianRender, CircleRender, EllipseRender
# Re-export RenderConfig from SMLMRender (used directly as step config)
const RenderConfig = SMLMRender.RenderConfig
export RenderConfig

# Re-export from SMLMClustering
const DBSCANConfig = SMLMClustering.DBSCANConfig
const HDBSCANConfig = SMLMClustering.HDBSCANConfig
const HierarchicalConfig = SMLMClustering.HierarchicalConfig
const VoronoiConfig = SMLMClustering.VoronoiConfig
const HopkinsConfig = SMLMClustering.HopkinsConfig
const VoronoiDensityConfig = SMLMClustering.VoronoiDensityConfig
export DBSCANConfig, HDBSCANConfig, HierarchicalConfig, VoronoiConfig
export HopkinsConfig, VoronoiDensityConfig

# Re-export from SMLMClustering — edge classification
const OuterPolygonConfig = SMLMClustering.OuterPolygonConfig
const KdeValleyConfig = SMLMClustering.KdeValleyConfig
export OuterPolygonConfig, KdeValleyConfig

# ============================================================
# Core types
# ============================================================
include("types.jl")
export Verbosity, Checkpoint
export AnalysisConfig
export stepinfo, stepinfos
export MultiTargetConfig
# AnalysisResult/AnalysisInfo/StepInfo/*Info structs, AbstractMultiTargetStep,
# MultiTargetResult/MultiTargetInfo, and the step_name/step_outdir extension
# hooks are not exported — reach them via SMLMAnalysis.Name (see
# api_overview.md's "Non-exported but public" section).

# ============================================================
# Step configs and pure step functions
# ============================================================

# Forward-declare analyze so step files can add dispatch methods
function analyze end

include("steps/common.jl")  # Shared helpers for steps (step_outdir is an
# extension hook, not exported — see workflows/extending.md)

include("steps/detectfit.jl")
export DetectFitConfig

include("steps/filter.jl")
export FilterConfig

include("steps/frameconnect.jl")
# FrameConnectConfig / CalibrationConfig are re-exported above (from SMLMFrameConnection)

include("steps/driftcorrect.jl")
# DriftConfig is re-exported above (from SMLMDriftCorrection)

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
# BaGoLConfig is re-exported above (from SMLMBaGoL)

include("steps/clustering.jl")
# Clustering config types (DBSCANConfig/HopkinsConfig/…) are re-exported above
# from SMLMClustering; this file only adds analyze() dispatch — no new exports.

include("steps/edgeclassify.jl")
# Edge-classify config types (OuterPolygonConfig/KdeValleyConfig) are re-exported
# above from SMLMClustering; this file only adds analyze() dispatch + a step_name
# override — no new exports here.

# ============================================================
# I/O
# ============================================================
include("io/atomic.jl")
include("io/smld_io.jl")
export save_smld, load_smld, smld_info

include("io/smart_h5.jl")
export load_smart_h5, load_smart_h5_info, smart_h5_to_array

include("io/mic_h5.jl")
export load_mic_h5, load_mic_h5_info, load_mic_h5_block,
    load_mic_h5_calibration, load_mic_h5_calibration_for_scmos
export build_camera_from_mic_h5

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
# AI coding-assistant guide installer
# ============================================================
include("agent_guide.jl")
export install_agent_guide, uninstall_agent_guide, agent_guide_status

# ============================================================
# Compact display for SMLMAnalysis-owned config types
# ============================================================
# Scoped to owned types only (see _show_config in types.jl). One method covers the
# owned AbstractMultiTargetStep tree (CompositeRender/CrossAlign/CrossCorr); the rest
# are the single-target step configs. Upstream const-aliases (FrameConnectConfig,
# DriftConfig, RenderConfig, BaGoLConfig, clustering configs) are deliberately absent
# — extending Base.show for them would be type piracy.
Base.show(io::IO, cfg::AbstractMultiTargetStep) = _show_config(io, cfg)
for T in (
        AnalysisConfig, MultiTargetConfig, DetectFitConfig,
        FilterConfig, DensityFilterConfig, IntensityFilterConfig,
    )
    @eval Base.show(io::IO, cfg::$T) = _show_config(io, cfg)
end

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
#   - backend = :cpu      → on BOTH the boxer and the fitter: no GPU kernels
#                           (uncacheable, no device on CI, and the boxer's GPU
#                           path polls NVML, which some GPUs do not support)
#   - outdir  = nothing   → no disk writes
#   - GaussianXYNBS       → Emitter2DFitSigma, the path the examples exercise
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
    (_, si) = simulate(
        sim;
        pattern = Nmer2D(n = 8, d = 0.05),
        molecule = GenericFluor(photons = 5.0e4, k_off = 20.0, k_on = 0.04),
        camera = cam
    )
    (imgs, _) = gen_images(
        si.smld_model, MicroscopePSFs.GaussianPSF(0.13);
        dataset = 1, bg = 20.0, poisson_noise = true
    )

    @compile_workload begin
        # Cached: the detect/fit → filter → frame-connect → render pipeline.
        cfg = AnalysisConfig(
            DetectFitConfig(
                boxer = BoxerConfig(boxsize = 7, psf_sigma = 0.13, backend = :cpu),
                fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu)
            ),
            FilterConfig(photons = (100.0, Inf)),
            FrameConnectConfig(max_frame_gap = 2),
            RenderConfig(zoom = 10);
            camera = cam,
            verbose = Verbosity.SILENT,
        )
        analyze([imgs], cfg)
    end
end


end # module
