"""
Density filter step - removes isolated localizations by neighbor count
"""

using NearestNeighbors

"""
    DensityFilterConfig <: AbstractSMLMConfig

Density-based filtering that removes isolated localizations lacking nearby neighbors.

# Keywords
- `n_sigma`: Search radius in localization uncertainty units (default: 2.0)
- `min_neighbors`: Minimum neighbor count. `:auto` uses valley detection between
  isolated and clustered populations (default: `:auto`)
"""
@kwdef struct DensityFilterConfig <: SMLMData.AbstractSMLMConfig
    n_sigma::Float64 = 2.0
    min_neighbors::Union{Int, Symbol} = :auto  # :auto uses valley detection (see _valley_threshold)
end

"""
    densityfilter_step(smld, cfg; outdir=nothing, step_number=0, verbose=Verbosity.STANDARD)

Filter localizations by neighbor density. Returns `(filtered_smld, DensityFilterInfo)`.

# Arguments
- `smld::BasicSMLD`: Input localizations
- `cfg::DensityFilterConfig`: Density filter parameters

# Keyword Arguments
- `outdir`: Output directory (nothing to skip file output)
- `step_number`: Step number for output directory naming
- `verbose`: Verbosity level

# Returns
`(filtered_smld, DensityFilterInfo)`
"""
function densityfilter_step(smld::BasicSMLD, cfg::DensityFilterConfig;
                            outdir::Union{String,Nothing}=nothing,
                            step_number::Int=0,
                            verbose::Int=Verbosity.STANDARD)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)

    v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg))" n_sigma=cfg.n_sigma min_neighbors=cfg.min_neighbors

    n_before = length(smld.emitters)
    t = @elapsed (filtered, neighbor_counts, threshold) = _filter_by_density(smld, cfg)
    n_after = length(filtered.emitters)
    n_rejected = n_before - n_after

    if dir !== nothing
        _save_densityfilter_outputs!(dir, cfg, v, t, neighbor_counts, threshold, n_before, n_after)
    end

    v >= Verbosity.PROGRESS && @info "  → $n_rejected rejected (threshold=$threshold) ($(round(t, digits=2))s)"
    (filtered, DensityFilterInfo(n_before, n_after, threshold, t))
end

_step_summary(info::DensityFilterInfo) = Dict{Symbol,Any}(
    :n_before => info.n_before,
    :n_after => info.n_after,
    :n_rejected => info.n_before - info.n_after,
    :threshold => info.threshold
)

"""
    analyze(smld, cfg::DensityFilterConfig; kwargs...) -> (filtered_smld, StepInfo)

Filter localizations by neighbor density.
"""
function analyze(smld::BasicSMLD, cfg::DensityFilterConfig;
                 outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD,
                 checkpoint::Int=Checkpoint.EXPENSIVE, kwargs...)
    t = @elapsed (filtered, df_info) = densityfilter_step(smld, cfg;
        outdir=outdir, step_number=step_number, verbose=verbose)

    if checkpoint >= Checkpoint.ALL
        dir = step_outdir(outdir, step_number, cfg)
        _save_step_smld(dir, filtered; filename="smld_density.jld2")
    end

    (filtered, StepInfo(step_number, cfg, t, _step_summary(df_info); info=df_info))
end

function _filter_by_density(smld::BasicSMLD, cfg::DensityFilterConfig)
    emitters = smld.emitters
    n = length(emitters)

    n == 0 && return smld, Int[], 0

    n_sigma = cfg.n_sigma
    σ = [sqrt(e.σ_x^2 + e.σ_y^2) for e in emitters]
    max_σ = maximum(σ)
    max_radius = n_sigma * 2 * max_σ

    coords = zeros(2, n)
    for i in 1:n
        coords[1, i] = emitters[i].x
        coords[2, i] = emitters[i].y
    end
    tree = KDTree(coords)

    neighbor_counts = zeros(Int, n)
    for i in 1:n
        point = [emitters[i].x, emitters[i].y]
        candidates = inrange(tree, point, max_radius)

        for j in candidates
            j == i && continue
            dist = sqrt((emitters[i].x - emitters[j].x)^2 + (emitters[i].y - emitters[j].y)^2)
            σ_combined = sqrt(σ[i]^2 + σ[j]^2)
            if dist < n_sigma * σ_combined
                neighbor_counts[i] += 1
            end
        end
    end

    threshold = if cfg.min_neighbors == :auto
        _valley_threshold(neighbor_counts)
    else
        cfg.min_neighbors
    end

    keep = neighbor_counts .>= threshold
    filtered = emitters[keep]

    BasicSMLD(filtered, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata), neighbor_counts, threshold
end

"""
Find threshold using valley detection between isolated (low neighbors) and clustered (high neighbors).

Algorithm:
1. Build and smooth histogram
2. Find rightmost significant peak (clustered population)
3. Find local minimum (valley) between origin and that peak
4. If no clear valley (unimodal), use heuristics based on peak position
"""
function _valley_threshold(counts::Vector{Int})
    isempty(counts) && return 1

    max_count = maximum(counts)
    max_count == 0 && return 1

    # Build histogram
    hist = zeros(Int, max_count + 1)
    for c in counts
        hist[c + 1] += 1
    end

    # Smooth histogram with simple moving average (window=3)
    smoothed = Float64.(hist)
    for i in 2:(length(hist)-1)
        smoothed[i] = (hist[i-1] + hist[i] + hist[i+1]) / 3
    end

    # Find the rightmost significant peak (clustered population):
    # a local maximum whose value is >= 5% of the max histogram value.
    peak_threshold = 0.05 * maximum(smoothed)
    rightmost_peak_idx = _rightmost_significant_peak(smoothed, peak_threshold)

    # If peak is at very low neighbor counts (< 5), distribution is mostly isolated
    # Use conservative threshold or warn
    if rightmost_peak_idx <= 5
        @warn "Neighbor distribution peaks at low values - most emitters appear isolated"
        return max(1, rightmost_peak_idx)
    end

    # Find valley (local minimum) between origin and the peak
    # Search from bin 1 to peak
    valley_idx = 1
    valley_val = smoothed[1]

    for i in 2:rightmost_peak_idx
        if smoothed[i] < valley_val
            valley_val = smoothed[i]
            valley_idx = i
        end
    end

    # Check if this is a real valley (bimodal) or just monotonic (unimodal)
    # Bimodal: valley should be significantly lower than both sides
    left_max = maximum(smoothed[1:valley_idx])
    right_max = maximum(smoothed[valley_idx:rightmost_peak_idx])

    # Valley is significant if it's below 70% of both adjacent maxima
    is_bimodal = valley_val < 0.7 * left_max && valley_val < 0.7 * right_max

    if is_bimodal
        # Clear valley found - use it as threshold
        return valley_idx - 1  # Convert to neighbor count (0-indexed)
    else
        # Unimodal distribution
        if smoothed[1] > smoothed[rightmost_peak_idx]
            # Peak at zero - mostly isolated, something may be wrong
            @warn "Distribution appears unimodal at low neighbor counts"
            return 3  # Conservative default
        else
            # Peak at high values - mostly clustered, little to filter
            return 1  # Keep almost everything
        end
    end
end

"""
    _rightmost_significant_peak(smoothed, peak_threshold) -> Int

Index of the rightmost local maximum of `smoothed` (a value ≥ both neighbours,
with a missing right edge treated as 0) whose height is at least `peak_threshold`.
Returns 1 if none qualifies. Shared by the valley-threshold search and the
diagnostic figure so both mark the same peak.
"""
function _rightmost_significant_peak(smoothed::AbstractVector{<:Real}, peak_threshold::Real)
    for i in length(smoothed):-1:3
        smoothed[i] >= smoothed[i-1] || continue
        smoothed[i] >= get(smoothed, i+1, 0.0) || continue
        smoothed[i] >= peak_threshold && return i
    end
    return 1
end

function _save_densityfilter_outputs!(dir::String, cfg::DensityFilterConfig, v::Int, t::Float64,
                                      neighbor_counts::Vector{Int}, threshold::Int, n_before::Int, n_after::Int)
    mkpath(dir)
    _save_config!(dir, cfg)
    _save_info!(dir, DensityFilterInfo(n_before, n_after, threshold, t))

    if v >= Verbosity.STANDARD
        _write_densityfilter_stats(dir, cfg, n_before, n_after, threshold, t)
        _save_densityfilter_figures(dir, neighbor_counts, threshold, cfg)
    end
end

function _write_densityfilter_stats(dir, cfg, n_before, n_after, threshold, t)
    n_rejected = n_before - n_after

    filepath = joinpath(dir, "stats.md")
    open(filepath, "w") do io
        println(io, "# Density Filter Statistics\n")
        println(io, "## Summary")
        println(io, "- **Input**: $n_before")
        println(io, "- **Output**: $n_after")
        println(io, "- **Rejected**: $n_rejected ($(round(100*n_rejected/n_before, digits=1))%)")
        println(io, "- **Threshold**: $threshold neighbors")
        println(io, "- **Time**: $(round(t, digits=2))s")
        println(io, "")
        println(io, "## Parameters")
        println(io, "- n_sigma: $(cfg.n_sigma)")
        println(io, "- min_neighbors: $(cfg.min_neighbors)")
    end
end

function _save_densityfilter_figures(dir, neighbor_counts, threshold, cfg)
    isempty(neighbor_counts) && return

    # Build histogram for visualization - use actual max, no artificial cap
    max_neighbor = maximum(neighbor_counts)
    hist_bins = zeros(Int, max_neighbor + 1)
    for c in neighbor_counts
        hist_bins[c + 1] += 1
    end

    # Find peak for visualization
    smoothed = Float64.(hist_bins)
    for i in 2:(length(hist_bins)-1)
        smoothed[i] = (hist_bins[max(1,i-1)] + hist_bins[i] + hist_bins[min(end,i+1)]) / 3
    end

    # Find rightmost significant peak (same routine the threshold uses, so the
    # marked peak matches the threshold decision).
    peak_threshold = 0.05 * maximum(smoothed)
    peak_idx = _rightmost_significant_peak(smoothed, peak_threshold)
    peak_val = hist_bins[peak_idx]

    method_str = cfg.min_neighbors == :auto ? "auto: valley method" : "manual"

    fig = Figure(size=(800, 500))
    ax = Axis(fig[1, 1],
        xlabel = "Number of Neighbors (within $(cfg.n_sigma)σ)",
        ylabel = "Count",
        title = "Neighbor Count Distribution (threshold = $threshold ($method_str))"
    )

    # Draw histogram bars
    barplot!(ax, 0:max_neighbor, hist_bins, color=:steelblue)

    # Draw line from origin to peak (showing valley search region)
    if cfg.min_neighbors == :auto && peak_idx > 1
        lines!(ax, [0, peak_idx - 1], [hist_bins[1], peak_val],
            color=:orange, linewidth=2, linestyle=:solid, label="Search region")
    end

    # Draw threshold
    vlines!(ax, [threshold - 0.5], color=:red, linestyle=:dash, linewidth=2, label="Threshold ($threshold)")

    # Add legend with stats
    n_rejected = sum(hist_bins[1:min(threshold, length(hist_bins))])
    n_kept = sum(hist_bins) - n_rejected
    pct_rejected = round(100 * n_rejected / sum(hist_bins), digits=1)

    Legend(fig[1, 2], ax, "Method",
        framevisible=true,
        padding=(10, 10, 10, 10))

    # Add stats text
    text!(ax, 0.95, 0.95,
        text="Threshold: $threshold neighbors\nRejected: $n_rejected ($pct_rejected%)\nKept: $n_kept",
        align=(:right, :top),
        space=:relative,
        fontsize=12)

    save(joinpath(dir, "neighbor_histogram.png"), fig)
end
