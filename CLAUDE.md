# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

SMLMAnalysis.jl is the high-level integration package for the JuliaSMLM ecosystem. It orchestrates all SMLM (Single Molecule Localization Microscopy) analysis packages into unified workflows with provenance tracking.

## JuliaSMLM Ecosystem

SMLMAnalysis depends on and orchestrates these JuliaSMLM packages. Each package has its own agent that can be contacted for coordination.

### Package Hierarchy

```
SMLMData (core types - no deps)
    |
+---+---+-----------+--------------+-----------------+---------+
|       |           |              |                 |         |
SMLMBoxer  GaussMLE  SMLMRender  SMLMFrameConnection  SMLMSim
|           |                                         |
+-----+-----+                                    MicroscopePSFs
      |
SMLMDriftCorrection (also depends on SMLMSim)
      |
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

**Always start Julia with `-t auto`** for multithreading (drift correction, frame connection, etc.):

```bash
# Run tests
julia -t auto --project=. -e 'using Pkg; Pkg.test()'

# Build documentation
julia -t auto --project=docs docs/make.jl

# Develop with local JuliaSMLM packages
julia --project=. -e 'using Pkg; Pkg.develop(path="../SMLMData")'

# Run examples (have their own Project.toml)
cd examples && julia -t auto --project=. stepwise_example.jl
```

**Test coverage note**: Tests are minimal (type construction and `analyze()` dispatch method existence checks). No integration tests with actual image data or GPU. The examples serve as de facto integration tests.

## Architecture

### Functional Pipeline

The package uses a unified `analyze()` dispatch architecture where the config type determines the operation. Steps are configured via typed config structs (`<: SMLMData.AbstractSMLMConfig`):

```julia
# Full pipeline with AnalysisConfig
config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(
            boxer=BoxerConfig(boxsize=9, psf_sigma=0.130),
            fitter=GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20)),
        FilterConfig(photons=(500.0, Inf)),
        FrameConnectConfig(max_frame_gap=5,
            calibration=CalibrationConfig(clamp_k_to_one=true)),
        DriftConfig(degree=2, dataset_mode=:registered),
        RenderConfig(zoom=20, colormap=:inferno),
    ],
    outdir = "output/",
)
(result, info) = analyze(image_stacks, config)

# Step-by-step with analyze() dispatch
(smld, df_info) = analyze(image_stacks, DetectFitConfig(
    camera=cam, boxer=BoxerConfig(boxsize=9, psf_sigma=0.130)))
(smld, f_info) = analyze(smld, FilterConfig(photons=(500.0, Inf)))
(smld, fc_info) = analyze(smld, FrameConnectConfig(max_frame_gap=5,
    calibration=CalibrationConfig(clamp_k_to_one=true)))
(smld, dc_info) = analyze(smld, DriftConfig(degree=2))
(img, r_info) = analyze(smld, RenderConfig(zoom=20))

# Save intermediate state for later resume
save_smld("output/after_detectfit.h5", smld)
smld = load_smld("output/after_detectfit.h5")
```

### Module Structure

```
src/
├── SMLMAnalysis.jl      # Main module, re-exports ecosystem types
├── types.jl             # AnalysisConfig, AnalysisResult, AnalysisInfo, DataSource, Verbosity
├── analysis.jl          # analyze() pipeline orchestrator
├── multitarget.jl       # MultiTargetConfig, multi-channel composite analysis
├── steps/               # One file per step type (analyze() dispatch + internal functions)
│   ├── common.jl        # Shared helpers (step_outdir, _save_config!, _save_info!)
│   ├── detectfit.jl     # DetectFitConfig, analyze(data, cfg) → (smld, info)
│   ├── filter.jl        # FilterConfig, analyze(smld, cfg) → (smld, info)
│   ├── frameconnect.jl  # Uses SMLMFrameConnection.FrameConnectConfig directly (calibration integrated)
│   ├── driftcorrect.jl  # Uses SMLMDriftCorrection.DriftConfig directly
│   ├── densityfilter.jl # DensityFilterConfig, analyze(smld, cfg) → (smld, info)
│   ├── render.jl        # analyze(smld, RenderConfig) → (image, info)
│   ├── bagol.jl         # BaGoLConfig (NOT included in module - legacy, not yet ported)
└── io/
    ├── smld_io.jl       # HDF5 serialization (save_smld, load_smld)
    ├── smart_h5.jl      # SMART microscope HDF5 import
    ├── mic_h5.jl        # MIC (MATLAB Instrument Control) H5 format import (block-based loading)
    └── checkpoint_io.jl # save_pipeline_state/load_pipeline_state (JLD2)
```

### Key Types

- **`AnalysisConfig`**: Complete pipeline description with camera, ordered steps vector, optional ROI, outdir, verbosity. Primary input to `analyze(data, config)`
- **`AnalysisResult`**: Immutable result from `analyze()` holding final `smld`, optional `smld_connected`, and `drift_model`
- **`AnalysisInfo`**: Aggregated metadata from all steps with `elapsed_s`, `steps` dict (step name → upstream info), and `step_infos` vector
- **`StepInfo <: AbstractSMLMInfo`**: Logged after each step with `elapsed_s`, `config`, `summary` dict, and typed `info::Union{AbstractSMLMInfo, Nothing}` field
- **`DataSource`**: Lazy loading wrapper - holds `images` (single array), `images_vec` (Vector of arrays for multi-dataset), or `path` (file for deferred loading)
- **`DetectFitInfo <: AbstractSMLMInfo`**: Info from detectfit step (boxes_info, fit_info, n_datasets, n_rois, n_fits)
- **`FilterInfo <: AbstractSMLMInfo`**: Info from filter step (n_before, n_after, elapsed_s)
- **`DensityFilterInfo <: AbstractSMLMInfo`**: Info from density filter step (n_before, n_after, threshold)
- **`Verbosity`**: Output detail levels (SILENT=0, PROGRESS=1, STANDARD=2, DETAILED=3, DEBUG=4)
- **`MultiTargetConfig`**: Multi-channel analysis config with per-channel labels, colors, composite rendering
- **`MultiTargetResult`**: Multi-channel result with per-channel `AnalysisResult` access via `result[:label]`

### Tuple-Pattern API

SMLMAnalysis follows the JuliaSMLM tuple-pattern where functions return `(result, info)` tuples:

```julia
# Pipeline analyze() returns (AnalysisResult, AnalysisInfo)
(result, info) = analyze(image_stacks, config)
result.smld               # Final SMLD
info.elapsed_s            # Total elapsed time in seconds
info.steps[:detectfit]    # DetectFitInfo from detectfit step
info.steps[:driftcorrect] # DriftInfo from SMLMDriftCorrection
info.step_infos           # Vector{StepInfo} with per-step timing/config/summary

# Step dispatch returns (result, StepInfo) with typed info
(smld, step_info) = analyze(image_stacks, DetectFitConfig(
    camera=cam, boxer=BoxerConfig(boxsize=9)))
step_info.info            # DetectFitInfo
step_info.elapsed_s       # Step timing

(smld, step_info) = analyze(smld, FrameConnectConfig())
step_info.info.connected  # Connected SMLD (FrameConnectInfo.connected)

(smld, step_info) = analyze(smld, DriftConfig(degree=2))
step_info.info.model      # Drift model (DriftInfo.model)
```

Each step uses config dispatch to upstream packages:
- `SMLMBoxer.getboxes(images, camera, cfg.boxer)` → `(ROIBatch, BoxesInfo)`
- `GaussMLE.fit(roi_batch, cfg.fitter)` → `(BasicSMLD, FitInfo)`
- `SMLMFrameConnection.frameconnect(smld, cfg)` → `(combined, FrameConnectInfo)` with `.connected` in info
- `SMLMDriftCorrection.driftcorrect(smld, cfg)` → `(corrected_smld, DriftInfo)` with `.model` in info
- `SMLMRender.render(smld, cfg)` → `(image, RenderInfo)`

### Data Flow

- All packages use SMLMData.jl types (BasicSMLD, Emitter2DFit, etc.)
- Coordinates are in microns
- Dataset boundaries encoded in data structure: `Vector{Array}` = N datasets, single `Array` = 1 dataset
- Each step function optionally writes outputs to `outdir/{step_number}_{step_name}/`
- Save/resume via `save_smld`/`load_smld` (HDF5) or `save_pipeline_state`/`load_pipeline_state` (JLD2)
- `smld_info(path)` prints file summary without loading data
- FilterConfig `precision` filter uses `max(e.σ_x, e.σ_y)`, not average or RMS

### Multi-Dataset Architecture

Dataset boundaries are encoded in the data structure, not user-specified integers:

- **`Vector{Array}`**: Pass `[dataset1, dataset2, ...]` where each element is a 3D image stack. The length determines `n_datasets`.
- **Single `Array`**: Treated as 1 dataset
- **MIC H5 format**: Blocks auto-detected as separate datasets
- **Multiple file paths**: `DetectFitConfig(paths=["d1.h5", "d2.h5"])` loads one file per dataset

```julia
# In-memory: data structure defines dataset boundaries
image_stacks = [images[:,:,1:2000], images[:,:,2001:4000]]  # 2 datasets
(result, info) = analyze(image_stacks, config)

# File-based: MIC format auto-detects blocks
config = AnalysisConfig(
    camera = cam,
    steps = [DetectFitConfig(path="data.h5", h5_format=:mic,
        boxer=BoxerConfig(boxsize=9, psf_sigma=0.130)), ...],
)
(result, info) = analyze(config)  # No data arg needed
```

For multi-dataset data:
- **Per-dataset processing**: Detection and fitting loop over datasets individually
- **Frame numbering**: Frames are per-dataset (1 to `n_frames_per_dataset`), NOT global. Required for LegendrePoly drift correction normalization to [-1, 1]
- **Dataset tracking**: `emitter.dataset` field tracks which dataset each emitter belongs to
- **SMLD structure**: `smld.n_frames` = frames per dataset, `smld.n_datasets` = number of datasets

### Drift Correction Modes

DriftConfig (from SMLMDriftCorrection) supports two primary use cases via the `dataset_mode` field.

**Continuous single acquisition** (one long movie):
```julia
# Chunked for long acquisitions (>4000 frames)
DriftConfig(
    degree = 3,
    dataset_mode = :continuous,
    chunk_frames = 4000,   # frames per chunk (alternative: n_chunks=N)
    auto_roi = true
)
```
- Consider chunking when >4000 frames; use `chunk_frames=4000` as reasonable max
- Chunking: `chunk_frames` (frames per chunk) OR `n_chunks` (number of chunks), not both
- `auto_roi=true` selects dense regions for better entropy signal

**Registered multi-dataset** (multiple files with stage registration):
```julia
DriftConfig(
    degree = 2,
    dataset_mode = :registered,
    quality = :singlepass
)
```
- Each dataset treated independently, then aligned via entropy optimization
- Requires spatial overlap between datasets for inter-dataset alignment
- Warns at PROGRESS verbosity if inter-dataset shifts exceed 500nm

### DetectFitConfig Data Sources

The combined detection+fitting step supports three data source modes:

1. **In-memory images**: Pass `Vector{Array}` to `analyze(image_stacks, config)` or `analyze(image_stacks, DetectFitConfig(camera=cam, boxer=BoxerConfig(boxsize=9)))`
2. **Single file with blocks**: `DetectFitConfig(path="data.h5", h5_format=:mic)` auto-detects MIC blocks as datasets
3. **Multiple files**: `DetectFitConfig(paths=["d1.h5", "d2.h5"])` loads one file per dataset

H5 formats auto-detected: `:smart` (SMART microscope), `:mic` (LidkeLab MIC format with block-based loading)

### Two-Layer Step Architecture

Each step has two layers:

**Internal function** (NOT exported) — does the work, returns `(result, info::AbstractSMLMInfo)`:
- `detectfit(data, camera, cfg)` → `(smld, DetectFitInfo)`
- `filter_step(smld, cfg)` → `(filtered, FilterInfo)`
- `frameconnect_step(smld, cfg)` → `(combined, FrameConnectInfo)` (upstream, calibration integrated via `cfg.calibration`)
- `driftcorrect_step(smld, cfg)` → `(corrected, DriftInfo)` (upstream, config dispatch)
- `densityfilter_step(smld, cfg)` → `(filtered, DensityFilterInfo)`

**analyze() dispatch** — thin wrapper that times the call and creates StepInfo:
```julia
function analyze(smld::BasicSMLD, cfg::DriftConfig; outdir=nothing, step_number=0, verbose=Verbosity.STANDARD)
    t = @elapsed (corrected, drift_info) = driftcorrect_step(smld, cfg; outdir, step_number, verbose)
    (corrected, StepInfo(step_number, cfg, t, _step_summary(drift_info); info=drift_info))
end
```

The orchestrator in `analysis.jl` calls `analyze()` for each step, collecting `StepInfo`s.

### Adding a New Step

1. Create `src/steps/yourstep.jl` with:
   - `@kwdef struct YourStepConfig <: SMLMData.AbstractSMLMConfig` with relevant fields
   - `struct YourStepInfo <: SMLMData.AbstractSMLMInfo` with step-specific fields
   - Internal function: `yourstep(smld, cfg; outdir=nothing, step_number=0, verbose=Verbosity.STANDARD)` returning `(result, YourStepInfo(...))`
   - `_step_summary(info::YourStepInfo)` returning `Dict{Symbol,Any}` for summary display
   - Dispatch method: `analyze(smld::BasicSMLD, cfg::YourStepConfig; kwargs...)` that times the call and returns `(result, StepInfo(...))`. The `kwargs...` is required so the pipeline can pass context like `smld_raw`.
2. Include it in `SMLMAnalysis.jl` and export the config and info types

No changes to `analysis.jl` needed — the pipeline loop is pure dispatch on `(state_type, config_type)`. Wrong step ordering gives a MethodError.

### Re-exported Types

From ecosystem packages (available after `using SMLMAnalysis`):
- **SMLMData**: AbstractCamera, IdealCamera, SCMOSCamera, Emitter2D/3D/2DFit/3DFit, BasicSMLD, ROIBatch, AbstractSMLMConfig, AbstractSMLMInfo
- **GaussMLE**: GaussMLEConfig, GaussianXYNB/S/SXSY, AstigmaticXYZNB, GaussMLEFitInfo, fit
- **SMLMBoxer**: getboxes, BoxerConfig (embedded in DetectFitConfig)
- **SMLMFrameConnection**: frameconnect, FrameConnectConfig, CalibrationConfig, CalibrationResult (all aliased via `const`)
- **SMLMDriftCorrection**: driftcorrect, DriftConfig (aliased as `const DriftConfig = SMLMDriftCorrection.DriftConfig`)
- **SMLMRender**: render, save_image, HistogramRender, GaussianRender, CircleRender, EllipseRender, RenderConfig (aliased as `const RenderConfig = SMLMRender.RenderConfig`)
- **SMLMSim**: StaticSMLMConfig, DiffusionSMLMConfig, simulate, gen_images, gen_image, Nmer2D, Nmer3D, Line2D, GenericFluor

## Uncertainty Calibration

Calibration is integrated into the frame connection step via `FrameConnectConfig.calibration`:
- Set `calibration=CalibrationConfig(clamp_k_to_one=true)` in `FrameConnectConfig` to enable
- Calibration runs inside `SMLMFrameConnection.frameconnect()`: link -> calibrate -> combine
- Results available via `FrameConnectInfo.calibration::CalibrationResult` (k_scale, sigma_motion_nm, mean_chi2, r_squared)
- Model: `sigma_corrected^2 = sigma_motion^2 + k^2 * sigma_CRLB^2`
- `CalibrationConfig` and `CalibrationResult` are re-exported from SMLMFrameConnection (not standalone steps)
