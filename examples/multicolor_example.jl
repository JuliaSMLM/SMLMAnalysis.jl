"""
    multicolor_example.jl

Multi-color (multi-target) SMLM analysis using the MultiTargetConfig API.

This example demonstrates:
1. Loading two simulated channels: Nmer2D clusters + Line2D filaments
2. Defining per-channel AnalysisConfig pipelines
3. Running analyze(channels, MultiTargetConfig) for composite rendering
4. Accessing per-channel results and composite outputs

Pipeline per channel:
  DetectFit -> Filter -> FrameConnect -> DriftCorrect -> Render

For single-channel workflows, see analysis_example.jl
"""

import Pkg
Pkg.activate(@__DIR__)

using SMLMAnalysis

include("generate_data.jl")

# ============================================================================
# Load or generate two-channel data
# ============================================================================

println("="^60)
println("Multi-color SMLM Analysis")
println("="^60)
println()

# Channel 1: Nmer2D octamer clusters (same data as single-target examples)
data_clusters = load_or_generate("single_target")
image_stacks_clusters = data_clusters["image_stacks"]  # Vector{SubArray}, dataset boundaries from data structure
camera = IdealCamera(data_clusters["camera_nx"], data_clusters["camera_ny"], data_clusters["camera_pixelsize"])
psf_sigma = data_clusters["psf_sigma"]

# Channel 2: Line2D filaments
data_lines = load_or_generate("lines")
image_stacks_lines = data_lines["image_stacks"]

println("Channel 1 (clusters): $(length(image_stacks_clusters)) datasets x $(size(image_stacks_clusters[1]))")
println("Channel 2 (lines):    $(length(image_stacks_lines)) datasets x $(size(image_stacks_lines[1]))")
println()

# ============================================================================
# Per-channel pipeline configs
# ============================================================================

config_clusters = AnalysisConfig(
    camera = camera,
    steps = [
        DetectFitConfig(
            boxer = BoxerConfig(boxsize=7, backend=:cpu),
            fitter = GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20, backend=:cpu),
        ),
        FilterConfig(
            photons = (500.0, Inf),
            precision = (0.0, 0.007),
            pvalue = (1e-3, 1.0),
        ),
        FrameConnectConfig(max_frame_gap = 5),
        CalibrationConfig(),
        DriftConfig(degree = 2),
        RenderConfig(zoom = 20, colormap = :inferno),
    ],
    verbose = Verbosity.STANDARD,
)

config_lines = AnalysisConfig(
    camera = camera,
    steps = [
        DetectFitConfig(
            boxer = BoxerConfig(boxsize=7, backend=:cpu),
            fitter = GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20, backend=:cpu),
        ),
        FilterConfig(
            photons = (500.0, Inf),
            precision = (0.0, 0.007),
            pvalue = (1e-3, 1.0),
        ),
        FrameConnectConfig(max_frame_gap = 5),
        CalibrationConfig(),
        DriftConfig(degree = 2),
        RenderConfig(zoom = 20, colormap = :inferno),
    ],
    verbose = Verbosity.STANDARD,
)

# ============================================================================
# Multi-target analysis
# ============================================================================

println("="^60)
println("Running multi-target analysis")
println("="^60)

mt = MultiTargetConfig(
    labels = [:clusters, :lines],
    colors = [:red, :green],
    render_zoom = 20,
    render_strategies = [GaussianRender(), HistogramRender()],
    outdir = joinpath(@__DIR__, "output", "multicolor_example"),
)

(result, info) = analyze([
    (image_stacks_clusters, config_clusters),
    (image_stacks_lines, config_lines),
], mt)

# ============================================================================
# Summary
# ============================================================================

println()
println("="^60)
println("Multi-Target Analysis Complete")
println("="^60)
println()
println("Results:")
for label in result.labels
    ch = result[label]
    println("  $label: $(length(ch.smld.emitters)) localizations")
end
total = sum(length(s.emitters) for s in result.smlds)
println("  Total: $total localizations")
println("  Time: $(round(info.elapsed_s, digits=2))s")
println()
println("Per-channel info:")
for (label, ch_info) in info.channels
    println("  $label: $(round(ch_info.elapsed_s, digits=2))s")
end
println()
println("Output: $(mt.outdir)")
println("  Per-channel: clusters/, lines/")
println("  Composite:   composite/")
println("  SMLDs:       smld_clusters.h5, smld_lines.h5")
