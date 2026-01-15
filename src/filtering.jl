"""
    filtering.jl

Filtering functions for SMLM localizations.
"""

using SMLMData
using NearestNeighbors
using Statistics

"""
    filter_smld(smld::BasicSMLD, config::AnalysisConfig) -> BasicSMLD

Filter SMLD based on config criteria (photons, precision, PSF sigma, pvalue).
"""
function filter_smld(smld::BasicSMLD, config::AnalysisConfig)
    emitters = smld.emitters
    mask = trues(length(emitters))

    if config.min_photons !== nothing
        mask .&= [e.photons > config.min_photons for e in emitters]
    end

    # Precision filter (localization precision from CRLB)
    if config.max_precision !== nothing
        mask .&= [max(e.σ_x, e.σ_y) < config.max_precision for e in emitters]
    end

    # PSF sigma mode filter (for variable sigma fits - filters on fitted PSF width)
    if config.psf_sigma_mode_tolerance !== nothing && length(emitters) > 0
        tol = config.psf_sigma_mode_tolerance

        # Check for isotropic sigma field (GaussianXYNBS)
        if hasfield(typeof(emitters[1]), :σ)
            psf_sigmas = [e.σ for e in emitters]
            psf_sigma_mode = _calculate_mode(psf_sigmas)

            if psf_sigma_mode > 0
                lo = psf_sigma_mode * (1 - tol)
                hi = psf_sigma_mode * (1 + tol)
                mask .&= [lo <= e.σ <= hi for e in emitters]
            end
        # Check for anisotropic sigma fields (GaussianXYNBSXSY)
        elseif hasfield(typeof(emitters[1]), :σx) && hasfield(typeof(emitters[1]), :σy)
            psf_sigmas_x = [e.σx for e in emitters]
            psf_sigmas_y = [e.σy for e in emitters]
            psf_sigma_mode_x = _calculate_mode(psf_sigmas_x)
            psf_sigma_mode_y = _calculate_mode(psf_sigmas_y)

            if psf_sigma_mode_x > 0 && psf_sigma_mode_y > 0
                lo_x = psf_sigma_mode_x * (1 - tol)
                hi_x = psf_sigma_mode_x * (1 + tol)
                lo_y = psf_sigma_mode_y * (1 - tol)
                hi_y = psf_sigma_mode_y * (1 + tol)
                mask .&= [lo_x <= e.σx <= hi_x && lo_y <= e.σy <= hi_y for e in emitters]
            end
        end
    end

    if config.min_pvalue !== nothing
        mask .&= [e.pvalue > config.min_pvalue for e in emitters]
    end

    filtered = emitters[mask]
    return BasicSMLD(filtered, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)
end

"""
    filter_isolated(smld::BasicSMLD, config::AnalysisConfig) -> (filtered_smld, neighbor_counts, threshold)

Filter out isolated emitters based on precision-weighted neighbor counting.
Two emitters are neighbors if their distance < n_sigma × sqrt(σ_i² + σ_j²).

When `isolated_min_neighbors = :auto`, uses triangle method to find optimal threshold.

Returns filtered SMLD, neighbor counts array, and the threshold used.
"""
function filter_isolated(smld::BasicSMLD, config::AnalysisConfig)
    emitters = smld.emitters
    n = length(emitters)

    if n == 0
        return smld, Int[], 0
    end

    n_sigma = config.isolated_n_sigma

    # Get precision for each emitter (localization uncertainty)
    σ = [sqrt(e.σ_x^2 + e.σ_y^2) for e in emitters]
    max_σ = maximum(σ)
    max_radius = n_sigma * 2 * max_σ  # Conservative search radius for KD-tree

    # Build KD-tree for fast neighbor search
    coords = zeros(2, n)
    for i in 1:n
        coords[1, i] = emitters[i].x
        coords[2, i] = emitters[i].y
    end
    tree = KDTree(coords)

    # Count precision-weighted neighbors for each emitter
    neighbor_counts = zeros(Int, n)
    for i in 1:n
        # Find candidates within max radius
        point = [emitters[i].x, emitters[i].y]
        candidates = inrange(tree, point, max_radius)

        for j in candidates
            j == i && continue
            dist = sqrt((emitters[i].x - emitters[j].x)^2 +
                       (emitters[i].y - emitters[j].y)^2)
            σ_combined = sqrt(σ[i]^2 + σ[j]^2)
            if dist < n_sigma * σ_combined
                neighbor_counts[i] += 1
            end
        end
    end

    # Determine threshold
    if config.isolated_min_neighbors == :auto
        min_neighbors = _triangle_threshold(neighbor_counts)
    else
        min_neighbors = config.isolated_min_neighbors
    end

    # Filter - keep emitters with enough neighbors
    keep = neighbor_counts .>= min_neighbors
    filtered = emitters[keep]

    return BasicSMLD(filtered, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata), neighbor_counts, min_neighbors
end

"""
    _triangle_threshold(counts) -> threshold

Triangle method for automatic threshold selection.

Finds the threshold that maximizes the perpendicular distance from a line
connecting the histogram peak to the histogram tail.

Works well for distributions with a peak at low values and a long tail
(like neighbor count histograms where most noise has 0-1 neighbors).
"""
function _triangle_threshold(counts::Vector{Int})
    if isempty(counts)
        return 1
    end

    max_count = maximum(counts)
    if max_count == 0
        return 1
    end

    # Build histogram (bins 0, 1, 2, ..., max_count)
    hist = zeros(Int, max_count + 1)
    for c in counts
        hist[c + 1] += 1  # +1 for 1-based indexing
    end

    # Find peak (mode) - usually at 0 or low values for noise
    peak_idx = argmax(hist)
    peak_val = hist[peak_idx]

    # Find last non-zero bin (tail end)
    tail_idx = findlast(x -> x > 0, hist)
    if tail_idx === nothing || tail_idx <= peak_idx
        return 1
    end

    # Line from peak to tail: (peak_idx, peak_val) to (tail_idx, hist[tail_idx])
    # Normalized line direction
    dx = tail_idx - peak_idx
    dy = hist[tail_idx] - peak_val
    line_len = sqrt(dx^2 + dy^2)
    if line_len == 0
        return 1
    end

    # Find point with maximum perpendicular distance from line
    max_dist = 0.0
    best_idx = peak_idx

    for i in peak_idx:tail_idx
        # Vector from peak to point i
        px = i - peak_idx
        py = hist[i] - peak_val

        # Perpendicular distance = |cross product| / line_length
        # cross = dx*py - dy*px (z-component of 3D cross product)
        cross = dx * py - dy * px
        dist = abs(cross) / line_len

        if dist > max_dist
            max_dist = dist
            best_idx = i
        end
    end

    # Return threshold (convert from 1-based histogram index to neighbor count)
    threshold = best_idx - 1  # Convert back to 0-based count

    # Ensure minimum threshold of 1 (at least 1 neighbor required)
    return max(1, threshold)
end

"""
    _calculate_mode(values; n_bins=100) -> mode_value

Calculate mode of values using histogram binning.
"""
function _calculate_mode(values::Vector{T}; n_bins=100) where T<:Real
    if isempty(values)
        return zero(T)
    end

    # Filter out NaN and Inf
    valid = filter(x -> isfinite(x) && x > 0, values)
    if isempty(valid)
        return zero(T)
    end

    # Create histogram bins
    lo, hi = quantile(valid, [0.01, 0.99])  # Robust range
    if lo >= hi
        return median(valid)  # Fall back to median
    end

    edges = range(lo, hi, length=n_bins+1)
    counts = zeros(Int, n_bins)

    for v in valid
        if lo <= v <= hi
            bin_idx = clamp(floor(Int, (v - lo) / (hi - lo) * n_bins) + 1, 1, n_bins)
            counts[bin_idx] += 1
        end
    end

    # Find mode (bin center with highest count)
    mode_idx = argmax(counts)
    mode_value = (edges[mode_idx] + edges[mode_idx+1]) / 2

    return T(mode_value)
end

"""Adaptive clip percentile based on number of localizations.
Higher percentile for dense data = less clipping = less saturation."""
function adaptive_clip_percentile(n_locs::Int)
    if n_locs < 50_000
        return 0.99
    elseif n_locs < 200_000
        return 0.995
    elseif n_locs < 500_000
        return 0.999
    else
        return 0.9999  # Very dense - minimal clipping
    end
end
