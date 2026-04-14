# Cell2 HEK DNA-PAINT 2-Channel Analysis
#
# SMART-format H5 data: 860x256x100k (2 channels side-by-side, 430px each)
# Camera: ORCA-Fusion (C14440-20UP) from gatta ruler calibration
# Channels: split at x=430, each gets a 256x256 ROI crop
#
# Usage: julia -t auto --project=dev dev/cell2_dnapaint_2ch.jl
#
# Runs each channel independently with ROI cropping.

import Pkg
Pkg.activate(@__DIR__)

using SMLMAnalysis
using Statistics

println("="^60)
println("CELL2 HEK DNA-PAINT 2-CHANNEL")
println("="^60)

# =============================================================================
# Data
# =============================================================================
include(joinpath(@__DIR__, "paths.jl"))
h5file = joinpath(dataroot("adapt"), "projects", "cells-labeling", "Data",
    "DNA PAINT smart", "2026-01-22",
    "Cell2_HEK_75pMimager_EGF_2ch--2026-01-22_18-09-34.h5")
isfile(h5file) || error("File not found: $h5file")

info = load_smart_h5_info(h5file)
println("Data: $(info.nframes) frames, $(info.width)x$(info.height)")
println("File size: $(round(info.file_size_gb, digits=2)) GB")

# =============================================================================
# Camera - ORCA-Fusion from gatta ruler characterization
# =============================================================================
# Full-frame camera covering the 860-wide 2ch image
camera = SCMOSCamera(info.width, info.height, 0.078f0, 0.7f0;
    offset=100.0f0, gain=0.24f0, qe=1.0f0)

# =============================================================================
# Load images (single continuous acquisition)
# =============================================================================
println("Loading images...")
t_load = @elapsed begin
    images, _ = smart_h5_to_array(h5file)
end
println("  Loaded $(size(images)) in $(round(t_load, digits=1))s")

image_stacks = [images]

# =============================================================================
# Channel geometry
# =============================================================================
# 860px wide = 2 channels x 430px, side-by-side in x (columns)
# Each channel gets a 256x256 ROI crop (matching papers-vortex-sr FOV config)
ch_width = info.width ÷ 2  # 430

# Channel ROIs: x ranges within the 860-wide image, y = full height (256)
# starty=100 from analysis_config.json → rows 101:356... but height is only 256
# So the y crop from the existing pipeline must be relative to the 430-wide channel
# For now: use full 256 height, 256 center crop in x within each 430-wide half
x_margin = (ch_width - 256) ÷ 2  # (430-256)/2 = 87

roi_ch1 = (x = (1 + x_margin):(256 + x_margin),
           y = 1:info.height)
roi_ch2 = (x = (ch_width + 1 + x_margin):(ch_width + 256 + x_margin),
           y = 1:info.height)

println("Channel 1 ROI: x=$(roi_ch1.x), y=$(roi_ch1.y)")
println("Channel 2 ROI: x=$(roi_ch2.x), y=$(roi_ch2.y)")

# =============================================================================
# Shared pipeline steps
# =============================================================================
# Split into pre-connection (for dual-view merge) and full pipeline
STEPS_PRE_CONNECTION = [
    DetectFitConfig(
        boxer=BoxerConfig(boxsize=11, min_photons=1000.0, psf_sigma=0.135),
        fitter=GaussMLEConfig(psf_model=GaussianXYNBSXSY(), iterations=20),
    ),
    FilterConfig(
        photons=(500.0, Inf),
        precision=(0.0, 0.010),
        pvalue=(1e-6, 1.0),
        psf_sigma=:auto,
    ),
]

STEPS_POST_CONNECTION = [
    FrameConnectConfig(
        max_frame_gap=5,
        max_sigma_dist=5.0,
        calibration=CalibrationConfig(clamp_k_to_one=true),
    ),
    DriftConfig(
        degree=2,
        dataset_mode=:continuous,
        n_chunks=20,
        quality=:iterative,
    ),
    DensityFilterConfig(
        n_sigma=2.0,
        min_neighbors=:auto,
    ),
    RenderConfig(zoom=20, colormap=:inferno, clip_percentile=0.999),
    RenderConfig(strategy=HistogramRender(), zoom=10, colormap=:turbo,
                 color_by=:absolute_frame, clip_percentile=nothing),
    RenderConfig(strategy=CircleRender(), zoom=50, colormap=:turbo,
                 color_by=:absolute_frame),
]

STEPS = vcat(STEPS_PRE_CONNECTION, STEPS_POST_CONNECTION)

# =============================================================================
# Run each channel
# =============================================================================
base_outdir = joinpath(@__DIR__, "output", "cell2_dnapaint_2ch")

channel_results = Dict{String, AnalysisResult}()

for (label, roi) in [("ch1", roi_ch1), ("ch2", roi_ch2)]
    println("\n" * "="^60)
    println("CHANNEL: $label  ROI: x=$(roi.x)")
    println("="^60)

    outdir = joinpath(base_outdir, label)
    rm(outdir; force=true, recursive=true)

    config = AnalysisConfig(
        camera=camera,
        steps=STEPS,
        roi=roi,
        outdir=outdir,
        verbose=Verbosity.DETAILED,
    )

    # Phase 1: detect + filter (saves pre-connection SMLD for dual-view merge)
    config_pre = AnalysisConfig(
        camera=camera, steps=STEPS_PRE_CONNECTION, roi=roi,
        outdir=outdir, verbose=Verbosity.DETAILED)
    (result_pre, _) = analyze(image_stacks, config_pre)
    save_smld(joinpath(outdir, "smld_preconnect.h5"), result_pre.smld)
    println("  Pre-connection: $(length(result_pre.smld.emitters)) localizations")

    # Phase 2: frameconnect + drift + densityfilter + render (step-by-step from SMLD)
    smld = result_pre.smld
    step_offset = length(STEPS_PRE_CONNECTION)
    local drift_model_ch = nothing
    for (i, cfg) in enumerate(STEPS_POST_CONNECTION)
        (smld_or_img, step_info) = analyze(smld, cfg;
            outdir=outdir, step_number=step_offset + i, verbose=Verbosity.DETAILED)
        if smld_or_img isa BasicSMLD
            smld = smld_or_img
        end
        if cfg isa DriftConfig && step_info.info !== nothing
            drift_model_ch = step_info.info.model
        end
    end
    result = AnalysisResult(smld, nothing, drift_model_ch)
    channel_results[label] = result

    save_smld(joinpath(outdir, "smld.h5"), result.smld; drift_model=result.drift_model)

    n = length(result.smld.emitters)
    println("\n  $label: $n localizations")
    if n > 0
        emitters = result.smld.emitters
        println("  Photons: median=$(round(median([e.photons for e in emitters]), digits=0))")
        println("  Precision: σ_x=$(round(median([e.σ_x for e in emitters])*1000, digits=1))nm")
    end
end

# =============================================================================
# Dual-View Channel Alignment and Merge
# =============================================================================
println("\n" * "="^60)
println("DUAL-VIEW MERGE")
println("="^60)

smld_ch1 = channel_results["ch1"].smld
smld_ch2 = channel_results["ch2"].smld

# --- Step 1: Flip ch2 x-coordinates (beam splitter mirror) ---
# Coordinates are in the full-frame camera system. The beam splitter mirror axis
# is at the center of the 860-pixel image (pixel 430). Flipping maps ch2 coords
# to ch1's coordinate space: x_flipped = 2 * mirror_center - x
pixel_size = 0.078  # μm/pixel
mirror_center_um = (860 ÷ 2) * pixel_size  # 430 * 0.078 = 33.54 μm
println("Flipping ch2 x-coordinates (mirror at $(round(mirror_center_um, digits=2)) μm)")
for e in smld_ch2.emitters
    e.x = 2 * mirror_center_um - e.x
end

# Debug: check coordinate ranges after flip
ch1_xs = [e.x for e in smld_ch1.emitters]
ch2_xs = [e.x for e in smld_ch2.emitters]
ch1_ys = [e.y for e in smld_ch1.emitters]
ch2_ys = [e.y for e in smld_ch2.emitters]
println("  ch1 x range: $(round(minimum(ch1_xs), digits=3)) - $(round(maximum(ch1_xs), digits=3)) μm")
println("  ch2 x range: $(round(minimum(ch2_xs), digits=3)) - $(round(maximum(ch2_xs), digits=3)) μm")
println("  ch1 y range: $(round(minimum(ch1_ys), digits=3)) - $(round(maximum(ch1_ys), digits=3)) μm")
println("  ch2 y range: $(round(minimum(ch2_ys), digits=3)) - $(round(maximum(ch2_ys), digits=3)) μm")

# Debug: check frame ranges
ch1_frames = [e.frame for e in smld_ch1.emitters]
ch2_frames = [e.frame for e in smld_ch2.emitters]
println("  ch1 frames: $(minimum(ch1_frames)) - $(maximum(ch1_frames))")
println("  ch2 frames: $(minimum(ch2_frames)) - $(maximum(ch2_frames))")

# --- Step 2: Align ch2 to ch1 via entropy cross-correlation ---
println("Aligning ch2 → ch1...")
(aligned_smlds, align_info) = align_smld([smld_ch1, smld_ch2], AlignConfig(method=:entropy, verbose=1))
smld_ch2_aligned = aligned_smlds[2]
shift = align_info.shifts[2]
println("  Shift: Δx=$(round(shift[1]*1000, digits=1))nm, Δy=$(round(shift[2]*1000, digits=1))nm")
println("  Total shift: $(round(sqrt(sum(shift.^2))*1000, digits=1))nm")

# --- Step 3: Match localizations across channels (frame + spatial proximity) ---
max_dist = 0.050  # 50nm matching radius in μm
max_frame_gap = 5  # Allow frame mismatch (connected tracks may start on different frames)
println("Matching localizations (max_dist=$(round(max_dist*1000, digits=0))nm, max_frame_gap=$max_frame_gap)...")

# Build per-frame index for ch2
ch2_by_frame = Dict{Int, Vector{Int}}()
for (i, e) in enumerate(smld_ch2_aligned.emitters)
    if !haskey(ch2_by_frame, e.frame)
        ch2_by_frame[e.frame] = Int[]
    end
    push!(ch2_by_frame[e.frame], i)
end

# Greedy nearest-neighbor matching with frame tolerance window
ch2_used = falses(length(smld_ch2_aligned.emitters))
matches = Tuple{Int,Int}[]  # (ch1_idx, ch2_idx)
unmatched_ch1 = Int[]

for (i, e1) in enumerate(smld_ch1.emitters)
    best_j = 0
    best_dist = max_dist
    # Search across frame window
    for df in -max_frame_gap:max_frame_gap
        for j in get(ch2_by_frame, e1.frame + df, Int[])
            ch2_used[j] && continue
            e2 = smld_ch2_aligned.emitters[j]
            d = sqrt((e1.x - e2.x)^2 + (e1.y - e2.y)^2)
            if d < best_dist
                best_dist = d
                best_j = j
            end
        end
    end
    if best_j > 0
        push!(matches, (i, best_j))
        ch2_used[best_j] = true
    else
        push!(unmatched_ch1, i)
    end
end

unmatched_ch2 = findall(.!ch2_used)

n_matched = length(matches)
n_ch1 = length(smld_ch1.emitters)
n_ch2 = length(smld_ch2_aligned.emitters)
println("  Matched: $n_matched pairs")
println("  Unmatched ch1: $(length(unmatched_ch1)), ch2: $(length(unmatched_ch2))")
println("  Match rate: $(round(100 * n_matched / max(n_ch1, 1), digits=1))% of ch1")

# --- Step 4: Precision-weighted merge of matched pairs ---
# Positions: inverse-variance weighted average
# Photons: sum (both channels see photons independently)
# Background: average
# Uncertainties: combined via quadrature
merged_emitters = Emitter2DFit{Float64}[]
sizehint!(merged_emitters, n_matched + length(unmatched_ch1) + length(unmatched_ch2))

# Helper to convert any emitter to Emitter2DFit
function _to_emitter2dfit(e)
    Emitter2DFit{Float64}(
        Float64(e.x), Float64(e.y), Float64(e.photons), Float64(e.bg),
        Float64(e.σ_x), Float64(e.σ_y), Float64(e.σ_photons), Float64(e.σ_bg);
        σ_xy=Float64(e.σ_xy), frame=e.frame, dataset=e.dataset,
        track_id=e.track_id, id=e.id)
end

for (i, j) in matches
    e1 = smld_ch1.emitters[i]
    e2 = smld_ch2_aligned.emitters[j]

    # Inverse-variance weights
    w1x = 1.0 / e1.σ_x^2
    w2x = 1.0 / e2.σ_x^2
    w1y = 1.0 / e1.σ_y^2
    w2y = 1.0 / e2.σ_y^2

    x_merged = (e1.x * w1x + e2.x * w2x) / (w1x + w2x)
    y_merged = (e1.y * w1y + e2.y * w2y) / (w1y + w2y)
    σ_x_merged = 1.0 / sqrt(w1x + w2x)
    σ_y_merged = 1.0 / sqrt(w1y + w2y)

    photons_merged = Float64(e1.photons + e2.photons)
    bg_merged = Float64((e1.bg + e2.bg) / 2)
    σ_photons_merged = sqrt(Float64(e1.σ_photons)^2 + Float64(e2.σ_photons)^2)
    σ_bg_merged = sqrt(Float64(e1.σ_bg)^2 + Float64(e2.σ_bg)^2) / 2

    push!(merged_emitters, Emitter2DFit{Float64}(
        x_merged, y_merged, photons_merged, bg_merged,
        σ_x_merged, σ_y_merged, σ_photons_merged, σ_bg_merged;
        σ_xy=0.0, frame=e1.frame, dataset=e1.dataset,
        track_id=e1.track_id, id=e1.id))
end

# Keep unmatched emitters from both channels
for i in unmatched_ch1
    push!(merged_emitters, _to_emitter2dfit(smld_ch1.emitters[i]))
end
for j in unmatched_ch2
    push!(merged_emitters, _to_emitter2dfit(smld_ch2_aligned.emitters[j]))
end

# Build merged SMLD using ch1's camera (both channels have identical cropped cameras)
smld_merged = BasicSMLD(merged_emitters, smld_ch1.camera,
    smld_ch1.n_frames, smld_ch1.n_datasets,
    Dict{String,Any}("source" => "dual-view merge",
                      "match_radius_um" => max_dist,
                      "align_shift" => shift))

# --- Step 5: Print precision improvement statistics ---
if n_matched > 0
    σ_x_ch1 = median([smld_ch1.emitters[i].σ_x for (i, _) in matches])
    σ_x_merged_med = median([merged_emitters[k].σ_x for k in 1:n_matched])
    improvement = σ_x_ch1 / σ_x_merged_med

    println("\nPrecision improvement (matched pairs):")
    println("  ch1 median σ_x:    $(round(σ_x_ch1*1000, digits=2)) nm")
    println("  merged median σ_x: $(round(σ_x_merged_med*1000, digits=2)) nm")
    println("  improvement:       $(round(improvement, digits=2))× (theoretical √2 ≈ 1.41)")
end

println("\nMerged SMLD: $(length(merged_emitters)) localizations")

# --- Step 6: Save merged SMLD ---
merge_outdir = joinpath(base_outdir, "merged")
mkpath(merge_outdir)
save_smld(joinpath(merge_outdir, "merged.h5"), smld_merged)
println("  Saved: $(joinpath(merge_outdir, "merged.h5"))")

# --- Step 7: Render merged result ---
println("\nRendering merged result...")
for (i, render_cfg) in enumerate([
    RenderConfig(zoom=20, colormap=:inferno, clip_percentile=0.999),
    RenderConfig(strategy=HistogramRender(), zoom=10, colormap=:turbo,
                 color_by=:absolute_frame, clip_percentile=nothing),
])
    (_, render_info) = analyze(smld_merged, render_cfg;
        outdir=merge_outdir, step_number=i, verbose=Verbosity.STANDARD)
end

println("\n" * "="^60)
println("COMPLETE - Output: $base_outdir")
println("="^60)
