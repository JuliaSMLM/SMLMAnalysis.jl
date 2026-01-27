"""
    calibration.jl

Uncertainty calibration functions for SMLM localizations.
Analyzes frame-to-frame drift from linked emitters to calibrate localization uncertainties.
"""

using SMLMData
using GaussMLE
using Statistics
using LinearAlgebra: det, inv, Diagonal

"""
    analyze_frameconnect_drift(smld_connected) -> NamedTuple

Analyze frame-to-frame drift from linked localizations.

For multi-dataset analysis, drift is calculated per dataset since frame numbers
reset between datasets and datasets may be acquired at different times.

Returns NamedTuple with:
- frame_shifts: Dict{Int, Vector} mapping dataset_id => shifts per frame transition
- chi2_values: Vector of χ² values for each pair (should follow χ²(2) if uncertainties correct)
- mean_chi2: Mean χ² (expected = 2 for correct uncertainties)
- n_pairs_total: Total number of frame-to-frame pairs analyzed
- n_datasets: Number of datasets with tracked emitters
- calibration: Fitted calibration model (A, B, uncertainties, R²)
- stage_motion: Stage motion estimates
"""
function analyze_frameconnect_drift(smld_connected)
    emitters = smld_connected.emitters
    n_datasets = smld_connected.n_datasets

    # Group emitters by track_id
    track_dict = Dict{Int, Vector{eltype(emitters)}}()
    for e in emitters
        if e.track_id > 0
            if !haskey(track_dict, e.track_id)
                track_dict[e.track_id] = eltype(emitters)[]
            end
            push!(track_dict[e.track_id], e)
        end
    end

    # Collect all frame-to-frame pairs, keyed by (dataset_id, frame)
    # Structure: (dataset, frame) => Vector of (Δx, Δy, var_x, var_y)
    frame_pairs = Dict{Tuple{Int, Int}, Vector{NTuple{4, Float64}}}()
    chi2_values = Float64[]

    for (track_id, track_emitters) in track_dict
        # Sort by dataset first, then frame
        sort!(track_emitters, by = e -> (e.dataset, e.frame))

        # Find consecutive frame pairs within same dataset
        for i in 1:(length(track_emitters) - 1)
            e1, e2 = track_emitters[i], track_emitters[i + 1]

            # Only consecutive frames within same dataset
            if e2.dataset == e1.dataset && e2.frame == e1.frame + 1
                Δx = Float64(e2.x - e1.x)
                Δy = Float64(e2.y - e1.y)
                var_x = Float64(e1.σ_x^2 + e2.σ_x^2)
                var_y = Float64(e1.σ_y^2 + e2.σ_y^2)

                key = (e1.dataset, e1.frame)
                if !haskey(frame_pairs, key)
                    frame_pairs[key] = NTuple{4, Float64}[]
                end
                push!(frame_pairs[key], (Δx, Δy, var_x, var_y))

                # Chi-squared for this pair
                if var_x > 0 && var_y > 0
                    χ2 = Δx^2 / var_x + Δy^2 / var_y
                    push!(chi2_values, χ2)
                end
            end
        end
    end

    # Calculate weighted mean shift for each (dataset, frame) transition
    # Output: Dict{dataset_id => Vector of (frame, Δx, Δy, σ_Δx, σ_Δy, n_pairs)}
    frame_shifts = Dict{Int, Vector{Tuple{Int, Float64, Float64, Float64, Float64, Int}}}()

    for (dataset_id, frame) in sort(collect(keys(frame_pairs)))
        pairs = frame_pairs[(dataset_id, frame)]
        n_pairs = length(pairs)

        if n_pairs > 0
            # Weighted average: weight = 1/variance
            sum_wx, sum_wy = 0.0, 0.0
            sum_w_x, sum_w_y = 0.0, 0.0

            for (Δx, Δy, var_x, var_y) in pairs
                if var_x > 0 && var_y > 0
                    w_x, w_y = 1.0 / var_x, 1.0 / var_y
                    sum_wx += w_x * Δx
                    sum_wy += w_y * Δy
                    sum_w_x += w_x
                    sum_w_y += w_y
                end
            end

            if sum_w_x > 0 && sum_w_y > 0
                Δx_mean = sum_wx / sum_w_x
                Δy_mean = sum_wy / sum_w_y
                σ_Δx = 1.0 / sqrt(sum_w_x)
                σ_Δy = 1.0 / sqrt(sum_w_y)

                if !haskey(frame_shifts, dataset_id)
                    frame_shifts[dataset_id] = Tuple{Int, Float64, Float64, Float64, Float64, Int}[]
                end
                push!(frame_shifts[dataset_id], (frame, Δx_mean, Δy_mean, σ_Δx, σ_Δy, n_pairs))
            end
        end
    end

    mean_chi2 = isempty(chi2_values) ? NaN : mean(chi2_values)

    # Collect all pair data for uncertainty calibration analysis
    # Each entry: (Δx², Δy², var_x, var_y)
    pair_data = NTuple{4, Float64}[]
    for pairs in values(frame_pairs)
        for (Δx, Δy, var_x, var_y) in pairs
            if var_x > 0 && var_y > 0
                push!(pair_data, (Δx^2, Δy^2, var_x, var_y))
            end
        end
    end

    # Fit uncertainty calibration model: observed_var = A + B * reported_var
    # Using simple linear regression on binned data
    calibration = _fit_uncertainty_calibration(pair_data)

    # Estimate stage motion from variance of mean shifts
    # If stage vibrates, all emitters move together → mean shift variance > expected from uncertainties
    # Collect all mean shifts across datasets
    all_Δx_mean = Float64[]
    all_Δy_mean = Float64[]
    all_σ_Δx = Float64[]
    all_σ_Δy = Float64[]
    for (dataset_id, shifts) in frame_shifts
        for (frame, Δx, Δy, σx, σy, n) in shifts
            push!(all_Δx_mean, Δx * 1000)  # Convert to nm
            push!(all_Δy_mean, Δy * 1000)
            push!(all_σ_Δx, σx * 1000)
            push!(all_σ_Δy, σy * 1000)
        end
    end

    # Observed variance of mean shifts
    var_Δx_observed = length(all_Δx_mean) > 1 ? var(all_Δx_mean) : 0.0
    var_Δy_observed = length(all_Δy_mean) > 1 ? var(all_Δy_mean) : 0.0

    # Expected variance from localization uncertainty alone
    var_Δx_expected = length(all_σ_Δx) > 0 ? mean(all_σ_Δx.^2) : 0.0
    var_Δy_expected = length(all_σ_Δy) > 0 ? mean(all_σ_Δy.^2) : 0.0

    # Excess variance = stage motion
    excess_var_x = max(0.0, var_Δx_observed - var_Δx_expected)
    excess_var_y = max(0.0, var_Δy_observed - var_Δy_expected)
    σ_stage_x = sqrt(excess_var_x)
    σ_stage_y = sqrt(excess_var_y)
    σ_stage = sqrt((excess_var_x + excess_var_y) / 2)  # Combined estimate

    stage_motion = (
        σ_stage_x = σ_stage_x,
        σ_stage_y = σ_stage_y,
        σ_stage = σ_stage,
        var_observed_x = var_Δx_observed,
        var_observed_y = var_Δy_observed,
        var_expected_x = var_Δx_expected,
        var_expected_y = var_Δy_expected,
        n_frames = length(all_Δx_mean)
    )

    return (
        frame_shifts = frame_shifts,
        chi2_values = chi2_values,
        mean_chi2 = mean_chi2,
        n_pairs_total = length(chi2_values),
        n_datasets = length(frame_shifts),
        pair_data = pair_data,
        calibration = calibration,
        stage_motion = stage_motion
    )
end

"""
    _fit_uncertainty_calibration(pair_data) -> NamedTuple

Fit uncertainty calibration model: observed_variance = A + B * CRLB_variance

Uses weighted least squares on raw data points (weight = 1/reported_var²) since
Var(Δ²) ∝ σ⁴. Binned data is computed separately for plotting only.

Returns NamedTuple with:
- A: additive term (nm²) - represents motion/vibration variance
- B: multiplicative factor - CRLB scale correction
- A_σ, B_σ: uncertainties on A and B
- r_squared: goodness of fit
- bin_centers, bin_observed, bin_expected: binned data for plotting
"""
function _fit_uncertainty_calibration(pair_data)
    if length(pair_data) < 100
        return (A = NaN, B = NaN, A_σ = NaN, B_σ = NaN, r_squared = NaN,
                bin_centers = Float64[], bin_observed = Float64[], bin_expected = Float64[],
                n_per_bin = Int[], bin_var = Float64[],
                n_filtered = 0, n_total = length(pair_data), chi2_threshold = 6.0)
    end

    # Combine x and y data: use average of x and y for each pair
    reported_var = [(p[3] + p[4]) / 2 for p in pair_data]  # Average of var_x and var_y
    observed_var = [(p[1] + p[2]) / 2 for p in pair_data]   # Average of Δx² and Δy²

    # Compute chi² for each pair to filter outliers
    # High chi² pairs are likely double-emitter fits where position shifts between frames
    chi2_per_pair = [p[1]/p[3] + p[2]/p[4] for p in pair_data]

    # Filter out pairs with chi² > threshold (likely double-fits or mismatches)
    chi2_threshold = 6.0  # ~99.7% of chi²(2) is below this
    good_mask = chi2_per_pair .<= chi2_threshold

    n_filtered = sum(.!good_mask)
    n_total = length(pair_data)

    # Apply filter
    reported_var = reported_var[good_mask]
    observed_var = observed_var[good_mask]

    if length(reported_var) < 100
        return (A = NaN, B = NaN, A_σ = NaN, B_σ = NaN, r_squared = NaN,
                bin_centers = Float64[], bin_observed = Float64[], bin_expected = Float64[],
                n_per_bin = Int[], bin_var = Float64[], n_filtered = n_filtered, n_total = n_total)
    end

    # Convert to nm² for interpretability
    reported_var_nm2 = reported_var .* 1e6
    observed_var_nm2 = observed_var .* 1e6

    # =========================================================================
    # First compute bins, then fit the binned data
    # Binning averages out chi-squared noise to reveal the underlying trend
    # Weight by n_per_bin for proper statistical weighting of bin means
    # =========================================================================
    x_min, x_max = extrema(reported_var_nm2)
    x_range = x_max - x_min
    n_bins = 20
    bin_width = x_range / n_bins

    bin_centers = Float64[]
    bin_observed = Float64[]
    bin_expected = Float64[]
    bin_var = Float64[]  # Within-bin variance for weighting
    n_per_bin = Int[]

    for i in 1:n_bins
        bin_lo = x_min + (i - 1) * bin_width
        bin_hi = x_min + i * bin_width

        # Include right edge in last bin
        if i == n_bins
            mask = (reported_var_nm2 .>= bin_lo) .& (reported_var_nm2 .<= bin_hi)
        else
            mask = (reported_var_nm2 .>= bin_lo) .& (reported_var_nm2 .< bin_hi)
        end

        bin_x = reported_var_nm2[mask]
        bin_y = observed_var_nm2[mask]

        if length(bin_x) >= 5  # Minimum points per bin
            push!(bin_centers, mean(bin_x))
            push!(bin_observed, mean(bin_y))
            push!(bin_expected, mean(bin_x))  # Expected if perfectly calibrated
            push!(bin_var, var(bin_y))  # Within-bin variance
            push!(n_per_bin, length(bin_x))
        end
    end

    # Need at least 3 bins for a meaningful fit
    n_valid_bins = length(bin_centers)
    if n_valid_bins < 3
        return (A = NaN, B = NaN, A_σ = NaN, B_σ = NaN, r_squared = NaN,
                bin_centers = bin_centers, bin_observed = bin_observed, bin_expected = bin_expected,
                n_per_bin = n_per_bin, bin_var = bin_var,
                n_filtered = n_filtered, n_total = n_total, chi2_threshold = chi2_threshold)
    end

    # =========================================================================
    # Weighted least squares on binned data
    # Weight = n / var_within_bin (inverse variance of bin mean)
    # This is proper statistical weighting for heteroscedastic data
    # =========================================================================

    # Variance of bin mean = var_within_bin / n
    var_of_bin_mean = bin_var ./ n_per_bin

    # Weights = 1 / variance_of_mean (inverse variance weighting)
    # Avoid division by zero
    weights = [v > 0 ? 1.0 / v : 0.0 for v in var_of_bin_mean]

    # Check we have valid weights
    if sum(weights) ≈ 0
        return (A = NaN, B = NaN, A_σ = NaN, B_σ = NaN, r_squared = NaN,
                bin_centers = bin_centers, bin_observed = bin_observed, bin_expected = bin_expected,
                n_per_bin = n_per_bin, bin_var = bin_var,
                n_filtered = n_filtered, n_total = n_total, chi2_threshold = chi2_threshold)
    end

    # Normalize weights
    weights = weights ./ sum(weights) .* n_valid_bins

    X = hcat(ones(n_valid_bins), bin_centers)
    W = Diagonal(weights)

    XtWX = X' * W * X
    XtWy = X' * W * bin_observed

    if det(XtWX) ≈ 0
        return (A = NaN, B = NaN, A_σ = NaN, B_σ = NaN, r_squared = NaN,
                bin_centers = bin_centers, bin_observed = bin_observed, bin_expected = bin_expected,
                n_per_bin = n_per_bin, bin_var = bin_var,
                n_filtered = n_filtered, n_total = n_total, chi2_threshold = chi2_threshold)
    end

    coeffs = XtWX \ XtWy
    A, B = coeffs[1], coeffs[2]

    # Compute weighted R² and parameter uncertainties
    y_pred = X * coeffs
    residuals = bin_observed .- y_pred

    ss_res = sum(weights .* residuals.^2)
    ss_tot = sum(weights .* (bin_observed .- mean(bin_observed)).^2)
    r_squared = ss_tot > 0 ? 1 - ss_res / ss_tot : NaN

    # Standard errors for weighted regression
    mse = ss_res / max(1, n_valid_bins - 2)
    var_coeffs = mse * inv(XtWX)
    A_σ = sqrt(max(0.0, var_coeffs[1, 1]))
    B_σ = sqrt(max(0.0, var_coeffs[2, 2]))

    return (A = A, B = B, A_σ = A_σ, B_σ = B_σ, r_squared = r_squared,
            bin_centers = bin_centers, bin_observed = bin_observed, bin_expected = bin_expected,
            n_per_bin = n_per_bin, bin_var = bin_var,
            n_filtered = n_filtered, n_total = n_total, chi2_threshold = chi2_threshold)
end

"""
    apply_uncertainty_calibration(smld_connected, σ_motion, k_scale) -> (smld_connected_corrected, smld_combined)

Apply uncertainty calibration to connected localizations and recombine tracks.

The calibration model is: σ²_corrected = σ²_motion + k² × σ²_CRLB

For each localization:
- σ_x_corrected = √(σ_motion² + k² × σ_x²)
- σ_y_corrected = √(σ_motion² + k² × σ_y²)

Then recombines tracks using weighted averaging with corrected uncertainties.

Returns:
- smld_connected_corrected: Connected localizations with corrected uncertainties
- smld_combined: Recombined tracks with corrected uncertainties
"""
function apply_uncertainty_calibration(smld_connected::BasicSMLD, σ_motion::Float64, k_scale::Float64)
    emitters = smld_connected.emitters

    # Apply correction to each localization
    σ_motion_sq = σ_motion^2
    k_sq = k_scale^2

    corrected_emitters = map(emitters) do e
        # Corrected uncertainties: σ_corrected = √(σ_motion² + k² × σ_CRLB²)
        σ_x_corrected = Float32(sqrt(σ_motion_sq + k_sq * e.σ_x^2))
        σ_y_corrected = Float32(sqrt(σ_motion_sq + k_sq * e.σ_y^2))

        # Create new emitter with corrected uncertainties using setproperties pattern
        _copy_emitter_with_uncertainty(e, σ_x_corrected, σ_y_corrected)
    end

    smld_connected_corrected = BasicSMLD(corrected_emitters, smld_connected.camera,
                                         smld_connected.n_frames, smld_connected.n_datasets,
                                         smld_connected.metadata)

    # Recombine tracks with corrected uncertainties
    smld_combined = recombine_tracks(smld_connected_corrected)

    return smld_connected_corrected, smld_combined
end

"""Copy emitter with new σ_x and σ_y values. Works for any GaussMLE emitter type."""
function _copy_emitter_with_uncertainty(e::GaussMLE.Emitter2DFitGaussMLE, σ_x::Float32, σ_y::Float32)
    GaussMLE.Emitter2DFitGaussMLE(
        e.x, e.y, e.photons, e.bg, σ_x, σ_y, e.σ_photons, e.σ_bg,
        e.pvalue, e.frame, e.dataset, e.track_id, e.id
    )
end

function _copy_emitter_with_uncertainty(e::GaussMLE.Emitter2DFitSigma, σ_x::Float32, σ_y::Float32)
    GaussMLE.Emitter2DFitSigma(
        e.x, e.y, e.photons, e.bg, e.σ, σ_x, σ_y, e.σ_photons, e.σ_bg, e.σ_σ,
        e.pvalue, e.frame, e.dataset, e.track_id, e.id
    )
end

function _copy_emitter_with_uncertainty(e::GaussMLE.Emitter2DFitSigmaXY, σ_x::Float32, σ_y::Float32)
    GaussMLE.Emitter2DFitSigmaXY(
        e.x, e.y, e.photons, e.bg, e.σx, e.σy, σ_x, σ_y, e.σ_photons, e.σ_bg, e.σ_σx, e.σ_σy,
        e.pvalue, e.frame, e.dataset, e.track_id, e.id
    )
end

# Fallback for SMLMData Emitter2DFit (no pvalue field)
function _copy_emitter_with_uncertainty(e::SMLMData.Emitter2DFit, σ_x::Float32, σ_y::Float32)
    SMLMData.Emitter2DFit(
        e.x, e.y, e.photons, e.bg, σ_x, σ_y, e.σ_xy, e.σ_photons, e.σ_bg,
        e.frame, e.dataset, e.track_id, e.id
    )
end

"""
    recombine_tracks(smld_connected) -> BasicSMLD

Recombine localizations into tracks using weighted averaging.

Groups localizations by track_id and computes weighted average position
and combined uncertainty for each track.

The combined emitter uses:
- Position: weighted average (weight = 1/σ²)
- Uncertainty: σ_combined = 1/√(Σ 1/σ²)
- Frame: middle frame of track
- Photons: sum of photons
- Background: mean background
- pvalue: geometric mean of pvalues
"""
function recombine_tracks(smld_connected::BasicSMLD)
    emitters = smld_connected.emitters
    EmitterType = eltype(emitters)

    # Group by track_id
    track_dict = Dict{Int, Vector{EmitterType}}()
    for e in emitters
        if e.track_id > 0
            if !haskey(track_dict, e.track_id)
                track_dict[e.track_id] = EmitterType[]
            end
            push!(track_dict[e.track_id], e)
        end
    end

    # Combine each track
    combined_emitters = EmitterType[]

    for (track_id, track_locs) in track_dict
        n = length(track_locs)

        # Weighted average position
        sum_wx, sum_wy = 0.0, 0.0
        sum_w_x, sum_w_y = 0.0, 0.0
        sum_photons = 0.0
        sum_bg = 0.0
        log_pvalue_sum = 0.0
        min_frame, max_frame = typemax(Int), 0
        dataset = track_locs[1].dataset

        for e in track_locs
            # Weights = 1/variance
            w_x = 1.0 / e.σ_x^2
            w_y = 1.0 / e.σ_y^2

            sum_wx += w_x * e.x
            sum_wy += w_y * e.y
            sum_w_x += w_x
            sum_w_y += w_y

            sum_photons += e.photons
            sum_bg += e.bg
            if hasproperty(e, :pvalue)
                log_pvalue_sum += log(max(1e-300, e.pvalue))
            end

            min_frame = min(min_frame, e.frame)
            max_frame = max(max_frame, e.frame)
        end

        # Combined position and uncertainty
        x_combined = Float32(sum_wx / sum_w_x)
        y_combined = Float32(sum_wy / sum_w_y)
        σ_x_combined = Float32(1.0 / sqrt(sum_w_x))
        σ_y_combined = Float32(1.0 / sqrt(sum_w_y))

        # Other combined properties
        frame_combined = div(min_frame + max_frame, 2)
        bg_combined = Float32(sum_bg / n)
        pvalue_combined = Float32(exp(log_pvalue_sum / n))  # geometric mean

        # Create combined emitter using first emitter as template
        template = track_locs[1]
        combined_e = _create_combined_emitter(template, x_combined, y_combined,
                                               Float32(sum_photons), bg_combined,
                                               σ_x_combined, σ_y_combined,
                                               frame_combined, dataset, track_id, pvalue_combined)
        push!(combined_emitters, combined_e)
    end

    return BasicSMLD(combined_emitters, smld_connected.camera,
                     smld_connected.n_frames, smld_connected.n_datasets,
                     smld_connected.metadata)
end

"""Create combined emitter from weighted average. Works for any GaussMLE emitter type."""
function _create_combined_emitter(template::GaussMLE.Emitter2DFitGaussMLE,
                                   x, y, photons, bg, σ_x, σ_y, frame, dataset, track_id, pvalue)
    GaussMLE.Emitter2DFitGaussMLE(x, y, photons, bg, σ_x, σ_y, σ_x, σ_y,  # σ_photons/σ_bg approximated
                                  pvalue, frame, dataset, track_id, 0)
end

function _create_combined_emitter(template::GaussMLE.Emitter2DFitSigma,
                                   x, y, photons, bg, σ_x, σ_y, frame, dataset, track_id, pvalue)
    σ_mean = (σ_x + σ_y) / 2
    GaussMLE.Emitter2DFitSigma(x, y, photons, bg, σ_mean, σ_x, σ_y, σ_x, σ_y, σ_mean,
                               pvalue, frame, dataset, track_id, 0)
end

function _create_combined_emitter(template::GaussMLE.Emitter2DFitSigmaXY,
                                   x, y, photons, bg, σ_x, σ_y, frame, dataset, track_id, pvalue)
    # For anisotropic, use mean of template's PSF sigmas as representative
    σx_psf = (template.σx + σ_x) / 2  # Approximate combined PSF sigma
    σy_psf = (template.σy + σ_y) / 2
    GaussMLE.Emitter2DFitSigmaXY(x, y, photons, bg, σx_psf, σy_psf, σ_x, σ_y,
                                  σ_x, σ_y, σ_x, σ_y,  # Approximations for σ_photons, etc.
                                  pvalue, frame, dataset, track_id, 0)
end

# Fallback for SMLMData Emitter2DFit
function _create_combined_emitter(template::SMLMData.Emitter2DFit,
                                   x, y, photons, bg, σ_x, σ_y, frame, dataset, track_id, pvalue)
    T = eltype(x)
    SMLMData.Emitter2DFit(x, y, photons, bg, σ_x, σ_y, zero(T), σ_x, σ_y,
                          frame, dataset, track_id, 0)
end
