# Cell2 Dual-View Merge (standalone)
#
# Loads saved per-channel SMLDs (post-pipeline, connected tracks) and merges.
# For DNA-PAINT dual-view, both channels image the same docking sites.
# Post-frame-connection tracks represent unique binding events; different
# temporal visits to the same site are at nearly identical positions.
#
# Matching strategy: spatial-only (no frame constraint) because connected
# tracks from different binding events at the same docking site have very
# different frame numbers (median |Δframe| ~28k). Spatial matching pairs
# each ch1 binding event with the nearest ch2 binding event at the same site.
#
# For true √2 precision improvement (same-photon merge), use pre-connection
# raw localizations where frame = actual camera frame. See cell2_dnapaint_2ch.jl
# which saves smld_preconnect.h5 for this purpose.
#
# Run cell2_dnapaint_2ch.jl first to generate the per-channel results.
#
# Usage: julia --project=dev dev/cell2_merge.jl

import Pkg
Pkg.activate(@__DIR__)

using SMLMAnalysis
using Statistics

base_outdir = joinpath(@__DIR__, "output", "cell2_dnapaint_2ch")

# =============================================================================
# Load per-channel results
# =============================================================================
println("Loading channel SMLDs...")
smld_ch1 = load_smld(joinpath(base_outdir, "ch1", "smld.h5"))
smld_ch2 = load_smld(joinpath(base_outdir, "ch2", "smld.h5"))
println("  ch1: $(length(smld_ch1.emitters)) localizations")
println("  ch2: $(length(smld_ch2.emitters)) localizations")

# =============================================================================
# Dual-View Channel Alignment and Merge
# =============================================================================
println("\n" * "="^60)
println("DUAL-VIEW MERGE")
println("="^60)

# --- Step 1: Flip ch2 x-coordinates (beam splitter mirror) ---
pixel_size = 0.078  # μm/pixel
mirror_center_um = (860 ÷ 2) * pixel_size  # 430 * 0.078 = 33.54 μm
println("Flipping ch2 x-coordinates (mirror at $(round(mirror_center_um, digits=2)) μm)")
for e in smld_ch2.emitters
    e.x = 2 * mirror_center_um - e.x
end

# Debug: coordinate ranges
println("  ch1 x: $(round(minimum(e.x for e in smld_ch1.emitters), digits=3)) - $(round(maximum(e.x for e in smld_ch1.emitters), digits=3)) μm")
println("  ch2 x: $(round(minimum(e.x for e in smld_ch2.emitters), digits=3)) - $(round(maximum(e.x for e in smld_ch2.emitters), digits=3)) μm")
println("  ch1 y: $(round(minimum(e.y for e in smld_ch1.emitters), digits=3)) - $(round(maximum(e.y for e in smld_ch1.emitters), digits=3)) μm")
println("  ch2 y: $(round(minimum(e.y for e in smld_ch2.emitters), digits=3)) - $(round(maximum(e.y for e in smld_ch2.emitters), digits=3)) μm")

# --- Step 2: Align ch2 to ch1 ---
# Give both SMLDs a shared camera so findshift's CC histogram works
# (the cropped ROI cameras have different pixel_edges, causing CC to fail)
println("Aligning ch2 → ch1...")
all_xs = vcat([e.x for e in smld_ch1.emitters], [e.x for e in smld_ch2.emitters])
all_ys = vcat([e.y for e in smld_ch1.emitters], [e.y for e in smld_ch2.emitters])
x_lo = floor(minimum(all_xs) / pixel_size) * pixel_size
x_hi = ceil(maximum(all_xs) / pixel_size) * pixel_size
y_lo = floor(minimum(all_ys) / pixel_size) * pixel_size
y_hi = ceil(maximum(all_ys) / pixel_size) * pixel_size
shared_cam = IdealCamera(
    Float32.(collect(range(x_lo, x_hi, step=Float64(pixel_size)))),
    Float32.(collect(range(y_lo, y_hi, step=Float64(pixel_size)))))

smld_ch1_align = BasicSMLD(smld_ch1.emitters, shared_cam, smld_ch1.n_frames, smld_ch1.n_datasets)
smld_ch2_align = BasicSMLD(smld_ch2.emitters, shared_cam, smld_ch2.n_frames, smld_ch2.n_datasets)

(aligned_smlds, align_info) = align_smld([smld_ch1_align, smld_ch2_align],
    AlignConfig(method=:fft, verbose=1))
smld_ch2_aligned = aligned_smlds[2]
shift = align_info.shifts[2]
println("  Best shift: Δx=$(round(shift[1]*1000, digits=1))nm, Δy=$(round(shift[2]*1000, digits=1))nm")
println("  Magnitude: $(round(sqrt(sum(shift.^2))*1000, digits=1))nm")

# Save AlignInfo
merge_outdir = joinpath(base_outdir, "merged")
rm(merge_outdir; force=true, recursive=true)
mkpath(merge_outdir)
open(joinpath(merge_outdir, "align_info.txt"), "w") do io
    println(io, "method: $(align_info.method)")
    println(io, "shift_x_nm: $(round(shift[1]*1000, digits=2))")
    println(io, "shift_y_nm: $(round(shift[2]*1000, digits=2))")
    println(io, "magnitude_nm: $(round(sqrt(sum(shift.^2))*1000, digits=2))")
    println(io, "elapsed_s: $(round(align_info.elapsed_s, digits=2))")
end

# --- Step 3: Frame + spatial matching ---
# Dual-view: both channels see the same binding event at the same time.
# Match by spatial proximity within a frame tolerance window.
max_dist = 0.050  # 50nm spatial matching radius
max_frame_gap = 5  # frame tolerance for connected tracks

println("\nMatching (max_dist=$(round(max_dist*1000))nm, max_frame_gap=$max_frame_gap)...")

function match_channels(emitters1, emitters2, max_dist, max_frame_gap)
    # Per-frame index for ch2
    ch2_by_frame = Dict{Int, Vector{Int}}()
    for (j, e) in enumerate(emitters2)
        flist = get!(Vector{Int}, ch2_by_frame, e.frame)
        push!(flist, j)
    end

    ch2_used = falses(length(emitters2))
    matched = Tuple{Int,Int}[]
    unmatched1 = Int[]

    for (i, e1) in enumerate(emitters1)
        best_j = 0
        best_d = max_dist
        for df in -max_frame_gap:max_frame_gap
            for j in get(ch2_by_frame, e1.frame + df, Int[])
                ch2_used[j] && continue
                e2 = emitters2[j]
                d = sqrt((e1.x - e2.x)^2 + (e1.y - e2.y)^2)
                if d < best_d
                    best_d = d
                    best_j = j
                end
            end
        end
        if best_j > 0
            push!(matched, (i, best_j))
            ch2_used[best_j] = true
        else
            push!(unmatched1, i)
        end
    end

    unmatched2 = findall(.!ch2_used)
    (matched, unmatched1, unmatched2)
end

(matches, unmatched_ch1, unmatched_ch2) = match_channels(
    smld_ch1.emitters, smld_ch2_aligned.emitters, max_dist, max_frame_gap)

n_matched = length(matches)
println("  Matched: $n_matched pairs")
println("  Unmatched ch1: $(length(unmatched_ch1)), ch2: $(length(unmatched_ch2))")
println("  Match rate: $(round(100 * n_matched / max(length(smld_ch1.emitters), 1), digits=1))% of ch1")

# --- Step 4: Precision-weighted merge ---
function merge_emitters(emitters1, emitters2, matches, unmatched1, unmatched2)
    merged = Emitter2DFit{Float64}[]
    sizehint!(merged, length(matches) + length(unmatched1) + length(unmatched2))

    for (i, j) in matches
        e1, e2 = emitters1[i], emitters2[j]
        w1x, w2x = 1.0 / e1.σ_x^2, 1.0 / e2.σ_x^2
        w1y, w2y = 1.0 / e1.σ_y^2, 1.0 / e2.σ_y^2

        push!(merged, Emitter2DFit{Float64}(
            (e1.x * w1x + e2.x * w2x) / (w1x + w2x),
            (e1.y * w1y + e2.y * w2y) / (w1y + w2y),
            Float64(e1.photons + e2.photons),
            Float64((e1.bg + e2.bg) / 2),
            1.0 / sqrt(w1x + w2x),
            1.0 / sqrt(w1y + w2y),
            sqrt(Float64(e1.σ_photons)^2 + Float64(e2.σ_photons)^2),
            sqrt(Float64(e1.σ_bg)^2 + Float64(e2.σ_bg)^2) / 2;
            σ_xy=0.0, frame=e1.frame, dataset=e1.dataset,
            track_id=e1.track_id, id=e1.id))
    end

    to_e2dfit(e) = Emitter2DFit{Float64}(
        Float64(e.x), Float64(e.y), Float64(e.photons), Float64(e.bg),
        Float64(e.σ_x), Float64(e.σ_y), Float64(e.σ_photons), Float64(e.σ_bg);
        σ_xy=Float64(e.σ_xy), frame=e.frame, dataset=e.dataset,
        track_id=e.track_id, id=e.id)

    # Unmatched emitters excluded — only matched pairs in final result
    merged
end

merged_emitters = merge_emitters(
    smld_ch1.emitters, smld_ch2_aligned.emitters,
    matches, unmatched_ch1, unmatched_ch2)

smld_merged = BasicSMLD(merged_emitters, smld_ch1.camera,
    smld_ch1.n_frames, smld_ch1.n_datasets,
    Dict{String,Any}("source" => "dual-view merge",
                      "match_radius_um" => max_dist,
                      "align_shift_nm" => shift .* 1000))

# --- Step 5: Statistics ---
if n_matched > 0
    σ_x_ch1 = median([smld_ch1.emitters[i].σ_x for (i, _) in matches])
    σ_x_merged = median([merged_emitters[k].σ_x for k in 1:n_matched])
    println("\nPrecision improvement (matched pairs):")
    println("  ch1 median σ_x:    $(round(σ_x_ch1*1000, digits=2)) nm")
    println("  merged median σ_x: $(round(σ_x_merged*1000, digits=2)) nm")
    println("  improvement:       $(round(σ_x_ch1/σ_x_merged, digits=2))× (√2 ≈ 1.41)")

    # Match distance distribution
    match_dists = [sqrt((smld_ch1.emitters[i].x - smld_ch2_aligned.emitters[j].x)^2 +
                        (smld_ch1.emitters[i].y - smld_ch2_aligned.emitters[j].y)^2) for (i,j) in matches]
    println("  Match distances: median=$(round(median(match_dists)*1000, digits=1))nm, " *
            "mean=$(round(mean(match_dists)*1000, digits=1))nm")
end

println("\nMerged SMLD: $(length(merged_emitters)) localizations")

# --- Step 6: Drift correction on merged data ---
println("\nDrift correcting merged SMLD...")
drift_cfg = DriftConfig(
    degree=2,
    dataset_mode=:continuous,
    n_chunks=20,
    quality=:iterative,
)
(smld_merged, dc_info) = analyze(smld_merged, drift_cfg;
    outdir=merge_outdir, step_number=1, verbose=Verbosity.STANDARD)

# --- Step 7: Save and render ---
save_smld(joinpath(merge_outdir, "merged.h5"), smld_merged)
println("  Saved: $(joinpath(merge_outdir, "merged.h5"))")

println("\nRendering merged result...")
for (i, render_cfg) in enumerate([
    RenderConfig(zoom=50, colormap=:inferno, clip_percentile=0.999),
    RenderConfig(strategy=HistogramRender(), zoom=20, colormap=:turbo,
                 color_by=:absolute_frame, clip_percentile=nothing),
])
    (_, _) = analyze(smld_merged, render_cfg;
        outdir=merge_outdir, step_number=i+1, verbose=Verbosity.STANDARD)
end

println("\n" * "="^60)
println("COMPLETE - Output: $merge_outdir")
println("="^60)
