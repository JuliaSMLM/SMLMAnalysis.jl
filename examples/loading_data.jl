"""
    loading_data.jl

How to get your raw camera data into `image_stacks` — the argument the Quick Start
and every other example passes to `analyze()`.

`image_stacks` is **raw camera data**: pixel intensities (ADU), *not* localizations.
It is one or more 3D arrays of shape `(height, width, frames)`. Each 3D array is one
*dataset* (a continuous acquisition):

    image_stacks = [stack1, stack2]   # 2 datasets, each (height, width, frames)
    image_stacks = stack              # 1 dataset — a single 3D array is also accepted

This script walks the three ways to obtain it:

  1. Stream from an H5 file (no in-memory array)  — the memory-efficient path
  2. Load an H5 file into memory                  — SMART / MIC formats
  3. Any other 3D array                           — TIFF, or simulated data

Only section 3 runs out of the box (it simulates its own data). Sections 1 and 2 are
guarded by `isfile` and print what to do once you point them at a real file.
"""

import Pkg
Pkg.activate(@__DIR__)

using SMLMAnalysis
using SMLMSim
import MicroscopePSFs

# A camera for the file-based examples. For real data, match your detector's pixel
# count and physical pixel size (microns). MIC files can build this for you (below).
cam = IdealCamera(256, 256, 0.1)   # 256x256, 100 nm pixels

# ============================================================================
# 1. Stream from files — no in-memory array (recommended for large acquisitions)
# ============================================================================
# Point DetectFitConfig at the file and call analyze(config) with NO data argument.
# The pipeline loads one dataset at a time, so a multi-GB acquisition never has to
# fit in RAM all at once. MIC files store multiple "blocks" that are auto-detected as
# separate datasets. AnalysisConfig always needs a camera; for MIC files you can build
# one from the file's own calibration with build_camera_from_mic_h5(path; pixel_size=).

const SMART_FILE = joinpath(@__DIR__, "data", "experiment_smart.h5")
const MIC_FILE   = joinpath(@__DIR__, "data", "experiment_mic.h5")

if isfile(MIC_FILE)
    println("[1] Streaming MIC file: $MIC_FILE")
    config = AnalysisConfig(
        camera = build_camera_from_mic_h5(MIC_FILE; pixel_size=0.1),  # or your own IdealCamera/SCMOSCamera
        steps = [
            DetectFitConfig(path=MIC_FILE, h5_format=:mic,
                            boxer=BoxerConfig(boxsize=9, psf_sigma=0.130),
                            fitter=GaussMLEConfig(psf_model=GaussianXYNBS())),
            FilterConfig(photons=(500.0, Inf)),
            RenderConfig(zoom=20, colormap=:inferno),
        ],
        outdir = joinpath(@__DIR__, "output", "loading_data_mic"),
    )
    (result, _) = analyze(config)           # no image_stacks argument
    println("    → $(length(result.smld.emitters)) localizations")
else
    println("[1] (skipped) No file at $MIC_FILE.")
    println("    For real data, stream it straight from disk:")
    println("      DetectFitConfig(path=\"your.h5\", h5_format=:mic, boxer=...)")
    println("      DetectFitConfig(paths=[\"d1.h5\", \"d2.h5\"], ...)   # one dataset per file")
    println("      analyze(config)   # no image_stacks argument")
end

# ============================================================================
# 2. Load an H5 file into memory
# ============================================================================
# Two microscope formats ship with loaders. Both return an (height, width, frames)
# array; wrap it in a Vector to form `image_stacks`. Use this when you want the raw
# frames in hand (e.g. to crop or inspect) before analysis.

if isfile(SMART_FILE)
    println("[2] Loading SMART file into memory: $SMART_FILE")
    stack, info = smart_h5_to_array(SMART_FILE)     # → (height, width, frames), info
    println("    stack size = $(size(stack))  ($(info.nframes) frames)")
    image_stacks = [stack]                          # one dataset
    # (then: analyze(image_stacks, config))
elseif isfile(MIC_FILE)
    println("[2] Loading MIC file into memory: $MIC_FILE")
    # A MIC file holds multiple blocks, each a separate dataset. load_mic_h5 returns
    # every block concatenated (+ a per-frame block id); wrapping that as [stack] would
    # merge all blocks into one dataset. Load blocks individually to keep them distinct:
    n_blocks = load_mic_h5_info(MIC_FILE).n_blocks
    image_stacks = [load_mic_h5_block(MIC_FILE, i) for i in 1:n_blocks]
    println("    loaded $(length(image_stacks)) block(s) as separate datasets")
else
    println("[2] (skipped) In-memory H5 loaders:")
    println("      # SMART — single acquisition:")
    println("      stack, _ = smart_h5_to_array(\"your.h5\");  image_stacks = [stack]")
    println("      # MIC — one stack per block (preserves dataset boundaries):")
    println("      n = load_mic_h5_info(\"your.h5\").n_blocks")
    println("      image_stacks = [load_mic_h5_block(\"your.h5\", i) for i in 1:n]")
end

# ============================================================================
# 3. Any 3D array — here, simulated data (runs out of the box)
# ============================================================================
# `analyze()` accepts any Array{<:Real,3} of (height, width, frames). That includes a
# TIFF stack read with TiffImages.jl, or — as here — frames synthesized by SMLMSim.
# This section actually runs so you can see the shape of `image_stacks` and confirm
# the end-to-end call.

println("[3] Simulating a small stack with SMLMSim ...")
sim_cam = IdealCamera(64, 64, 0.1)
sim = StaticSMLMConfig(density=2.0, σ_psf=0.13, nframes=200, ndatasets=1)
(_, si) = simulate(sim;
    pattern  = Nmer2D(n=8, d=0.05),
    molecule = GenericFluor(photons=5.0e4, k_off=20.0, k_on=0.05),
    camera   = sim_cam)
(frames, _) = gen_images(si.smld_model, MicroscopePSFs.GaussianPSF(0.13);
    dataset=1, bg=20.0, poisson_noise=true)

# `frames` is a plain 3D array — exactly what a TIFF or H5 loader would hand you.
println("    frames::$(typeof(frames))  size = $(size(frames))  → (height, width, frames)")

image_stacks = [frames]           # wrap as a 1-dataset Vector (a bare `frames` also works)
config = AnalysisConfig(
    camera = sim_cam,
    steps  = [
        DetectFitConfig(boxer=BoxerConfig(boxsize=7, psf_sigma=0.13),
                        fitter=GaussMLEConfig(psf_model=GaussianXYNBS())),
        FilterConfig(photons=(500.0, Inf)),
        RenderConfig(zoom=10, colormap=:inferno),
    ],
    outdir  = joinpath(@__DIR__, "output", "loading_data_sim"),
    verbose = Verbosity.PROGRESS,
)
(result, analysis_info) = analyze(image_stacks, config)

println()
println("Done. image_stacks = Vector of $(length(image_stacks)) dataset(s); " *
        "each element is a $(ndims(frames))D $(eltype(frames)) array of (height, width, frames).")
println("→ $(length(result.smld.emitters)) localizations in " *
        "$(round(analysis_info.elapsed_s, digits=2))s.")
