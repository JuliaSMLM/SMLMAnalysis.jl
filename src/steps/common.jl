"""
Common helper functions shared across analysis steps.
"""

# Fallback summary dispatch - overridden per step
_step_summary(::SMLMData.AbstractSMLMInfo) = Dict{Symbol, Any}()

"""
    _calculate_mode(values; n_bins=100)

Calculate the mode of a distribution using histogram binning.
Returns the center of the most populated bin.

Uses a median-centered range (median ± 3×MAD) to avoid outlier peaks
at large fitted PSF sigma pulling the mode away from the true peak.
"""
function _calculate_mode(values::Vector{T}; n_bins=100) where T<:Real
    isempty(values) && return zero(T)

    valid = filter(x -> isfinite(x) && x > 0, values)
    isempty(valid) && return zero(T)

    med = median(valid)
    mad_val = median(abs.(valid .- med))
    mad_val == 0 && return med

    # Median-centered range: captures the primary peak, excludes outlier clusters
    lo = max(med - 3 * mad_val, minimum(valid))
    hi = med + 3 * mad_val
    lo >= hi && return med

    edges = range(lo, hi, length=n_bins+1)
    counts = zeros(Int, n_bins)

    for v in valid
        if lo <= v <= hi
            bin_idx = clamp(floor(Int, (v - lo) / (hi - lo) * n_bins) + 1, 1, n_bins)
            counts[bin_idx] += 1
        end
    end

    mode_idx = argmax(counts)
    T((edges[mode_idx] + edges[mode_idx+1]) / 2)
end

"""
    _grid_figure_size(data; n_cols=4, n_rows=3, panel_height=200)

Calculate figure size for grid overlay plots based on data aspect ratio.
"""
function _grid_figure_size(data; n_cols=4, n_rows=3, panel_height=200)
    data_height, data_width = size(data, 1), size(data, 2)
    data_aspect = data_width / data_height
    panel_width = round(Int, panel_height * data_aspect)
    fig_width = panel_width * n_cols + 100
    fig_height = panel_height * n_rows + 150
    (fig_width, fig_height)
end

"""
    _save_box_overlay(dir, filename, images, x_corners, y_corners, frame_indices, box_size, box_colors; ...)

Core box-overlay renderer: a grid of sample frames (contrast-stretched grayscale) with a
colored `box_size` rectangle drawn at each `(x_corner, y_corner)` on its frame. This is the
shared "boxes on the raw data" diagnostic style — detectfit (`detection_overlay.png`) and filter
(`fit_overlay.png`) render through it so the whole pipeline's overlays are one visual family.

The box-geometry is passed as plain arrays so non-ROIBatch callers use the identical
renderer; see the `roi_batch` convenience method below.
"""
function _save_box_overlay(dir, filename, images, x_corners, y_corners, frame_indices, box_size,
                           box_colors; title_prefix="Frame", frame_labels=nothing, suptitle=nothing)
    n_frames = size(images, 3)
    sample = [round(Int, x) for x in range(1, n_frames, length=min(12, n_frames))]
    display_labels = frame_labels !== nothing ? frame_labels : sample

    # Contrast stretch: dark background, spot cores retain structure
    sample_frames = sample[1:min(4, length(sample))]
    sample_data = vec(images[:, :, sample_frames])
    pmin = Float64(quantile(sample_data, 0.25))    # Below background -> solid black
    pmax = Float64(quantile(sample_data, 0.9995))  # Above most spot peaks -> preserve core detail

    fig = Figure(size=_grid_figure_size(images))
    if suptitle !== nothing
        Label(fig[0, 1:4], suptitle, fontsize=11)
    end

    for (idx, frame_num) in enumerate(sample)
        row = div(idx - 1, 4) + 1
        col = mod(idx - 1, 4) + 1

        ax = Axis(fig[row, col], title="$title_prefix $(display_labels[idx])", aspect=DataAspect(), yreversed=true)
        heatmap!(ax, images[:, :, frame_num]', colormap=:grays, colorrange=(pmin, pmax))

        frame_mask = frame_indices .== frame_num
        if any(frame_mask)
            for (x, y, c) in zip(x_corners[frame_mask], y_corners[frame_mask], box_colors[frame_mask])
                lines!(ax, [x, x+box_size, x+box_size, x, x],
                          [y, y, y+box_size, y+box_size, y],
                    color=c, linewidth=0.5)
            end
        end
        hidedecorations!(ax)
    end

    save(joinpath(dir, filename), fig)
end

"""
    _save_box_overlay(dir, filename, images, roi_batch, box_colors; ...)

ROIBatch convenience (detectfit/filter): unpack `x_corners`/`y_corners`/`frame_indices`/`roi_size`
and delegate to the core renderer above. Behavior unchanged for existing callers.
"""
function _save_box_overlay(dir, filename, images, roi_batch, box_colors; kwargs...)
    _save_box_overlay(dir, filename, images, roi_batch.x_corners, roi_batch.y_corners,
                      roi_batch.frame_indices, roi_batch.roi_size, box_colors; kwargs...)
end

"""
    _estimate_bleaching_rate(frame_counts) -> NamedTuple or nothing

Fit `N(t) = a + b*exp(-k*t)` to a per-frame localization count vector via
linearized initial guess + Nelder-Mead refinement.

This is the OBSERVED decay rate, not the true k_bleach for GenericFluor:
bleaching only occurs from the On state, so `k_observed = k_bleach * P_on`,
where `P_on = k_on/(k_on+k_off)` is the duty cycle. To get true k_bleach:
`k_bleach = k_observed / P_on`.

Returns `(k_bleach, N_0, offset, half_life, r_squared, valid_frames, smoothed)`
or `nothing` if the fit fails.
"""
function _estimate_bleaching_rate(frame_counts::Vector{Int})
    valid_mask = frame_counts .> 0
    valid_frames = findall(valid_mask)
    valid_counts = frame_counts[valid_mask]

    length(valid_counts) < 10 && return nothing

    window = min(50, length(valid_counts) ÷ 10)
    if window > 1
        smoothed = [mean(valid_counts[max(1, i-window):min(end, i+window)]) for i in 1:length(valid_counts)]
    else
        smoothed = Float64.(valid_counts)
    end

    # Initial guess from linearized fit
    tail_start = max(1, round(Int, 0.9 * length(smoothed)))
    a0 = mean(smoothed[tail_start:end])

    shifted = smoothed .- a0
    pos_mask = shifted .> 0
    sum(pos_mask) < 10 && return nothing

    x_lin = Float64.(valid_frames[pos_mask])
    y_lin = log.(shifted[pos_mask])

    n = length(x_lin)
    sum_x = sum(x_lin)
    sum_y = sum(y_lin)
    sum_xy = sum(x_lin .* y_lin)
    sum_x2 = sum(x_lin .^ 2)
    denom = n * sum_x2 - sum_x^2
    abs(denom) < 1e-10 && return nothing

    slope = (n * sum_xy - sum_x * sum_y) / denom
    intercept = (sum_y - slope * sum_x) / n
    k0 = -slope
    b0 = exp(intercept)
    k0 <= 0 && return nothing

    # Nelder-Mead refinement with physical non-negativity constraints on (a, b, k).
    # The model a + b*exp(-k*t) has a flat-data degeneracy where (a ≪ 0, b ≫ 0, k ≈ 0)
    # fits the mean equally well as any physical solution — R² stays high while the fit
    # is meaningless. Inf penalties on negative a or b block that basin.
    t = Float64.(valid_frames)
    function cost(p)
        a, b, k = p
        (k <= 0 || a < 0 || b < 0) && return Inf
        pred = a .+ b .* exp.(-k .* t)
        sum((smoothed .- pred) .^ 2)
    end

    result = Optim.optimize(cost, [a0, b0, k0], Optim.NelderMead(),
                            Optim.Options(iterations=5000, g_tol=1e-8))

    a_fit, b_fit, k_fit = Optim.minimizer(result)
    (k_fit <= 0 || a_fit < 0 || b_fit < 0) && return nothing

    half_life = log(2) / k_fit

    # Degeneracy guards: reject fits that are effectively constant (no decay over window)
    # or where the decaying component is negligible compared to baseline.
    max_val = maximum(smoothed)
    n_window = length(smoothed)
    half_life > 5 * n_window && return nothing
    b_fit < 0.01 * max_val && return nothing

    y_pred = a_fit .+ b_fit .* exp.(-k_fit .* t)
    ss_res = sum((smoothed .- y_pred) .^ 2)
    ss_tot = sum((smoothed .- mean(smoothed)) .^ 2)
    r_squared = ss_tot > 0 ? 1 - ss_res / ss_tot : 0.0

    (k_bleach=k_fit, N_0=b_fit, offset=a_fit, half_life=half_life, r_squared=r_squared,
     valid_frames=valid_frames, smoothed=smoothed)
end

"""
    _save_loc_per_frame(dir, smld; filename, title) -> bleach_result or nothing

Plot localizations per absolute frame across all datasets, with photobleaching
exponential decay fit overlay (`a + b*exp(-k*t)`). Marks dataset boundaries for
multi-dataset SMLDs.

Used by detectfit (raw fits) and filter (post-filter fits) steps.
"""
function _save_loc_per_frame(dir::String, smld::BasicSMLD;
                              filename::String="localizations_per_frame.png",
                              title::String="Localizations per Frame")
    emitters = smld.emitters
    isempty(emitters) && return nothing

    n_frames = smld.n_frames
    n_total = n_frames * smld.n_datasets
    frame_counts = zeros(Int, n_total)
    for e in emitters
        abs_frame = (e.dataset - 1) * n_frames + e.frame
        if abs_frame >= 1 && abs_frame <= n_total
            frame_counts[abs_frame] += 1
        end
    end

    bleach_result = _estimate_bleaching_rate(frame_counts)

    fig = Figure(size=(900, 400))
    ax = Axis(fig[1, 1], xlabel="Absolute Frame", ylabel="Localizations", title=title)
    lines!(ax, 1:n_total, frame_counts, color=(:blue, 0.5), linewidth=0.5, label="Raw")

    if bleach_result !== nothing
        lines!(ax, bleach_result.valid_frames, bleach_result.smoothed,
               color=:blue, linewidth=1.5, label="Smoothed")

        k = round(bleach_result.k_bleach, sigdigits=3)
        tau = round(bleach_result.half_life, digits=0)
        a = round(Int, bleach_result.offset)
        R2 = round(bleach_result.r_squared, digits=3)
        fit_frames = 1:n_total
        fit_counts = bleach_result.N_0 .* exp.(-bleach_result.k_bleach .* fit_frames) .+ bleach_result.offset
        lines!(ax, fit_frames, fit_counts, color=:red, linewidth=2, linestyle=:dash,
               label="Fit: k=$k/frame, t1/2=$(Int(tau)), a=$a, R^2=$R2")

        hlines!(ax, [bleach_result.offset], color=(:red, 0.3), linestyle=:dot, linewidth=1)
    else
        hlines!(ax, [mean(frame_counts)], color=:red, linestyle=:dash,
                label="mean ($(round(mean(frame_counts), digits=1)))")
    end

    axislegend(ax, position=:rt, framevisible=false, labelsize=10)

    if smld.n_datasets > 1
        for ds in 2:smld.n_datasets
            vlines!(ax, [(ds - 1) * n_frames + 0.5], color=(:gray, 0.5), linestyle=:dash)
        end
    end

    save(joinpath(dir, filename), fig)
    return bleach_result
end

# ============================================================
# Sample frame planning (for overlay plots)
# ============================================================

"""
    _plan_sample_frames(ds_frame_counts::Vector{Int}, n_samples::Int=12)

Precompute which frames to sample across all datasets for overlay plots.
Always includes first and last absolute frames, with equal spacing between.

Returns `(plan, abs_frames)`:
- `plan::Dict{Int, Vector{Int}}`: dataset index → local frame indices to capture
- `abs_frames::Vector{Int}`: absolute frame numbers for display labels (in order)
"""
function _plan_sample_frames(ds_frame_counts::Vector{Int}, n_samples::Int=12)
    total = sum(ds_frame_counts)
    total == 0 && return Dict{Int,Vector{Int}}(), Int[]

    n = min(n_samples, total)
    abs_targets = unique([round(Int, x) for x in range(1, total, length=n)])

    cumulative = cumsum(ds_frame_counts)
    plan = Dict{Int,Vector{Int}}()
    abs_frames = Int[]

    for abs_frame in abs_targets
        ds = findfirst(c -> abs_frame <= c, cumulative)
        local_frame = ds == 1 ? abs_frame : abs_frame - cumulative[ds-1]
        if !haskey(plan, ds)
            plan[ds] = Int[]
        end
        push!(plan[ds], local_frame)
        push!(abs_frames, abs_frame)
    end

    (plan, abs_frames)
end

# ============================================================
# Dataset assignment helpers
# ============================================================

"""Update emitter's dataset field using struct reconstruction"""
function _with_dataset(e::Emitter2DFit{T}, ds::Int) where T
    Emitter2DFit{T}(
        e.x, e.y, e.photons, e.bg, e.σ_x, e.σ_y, e.σ_photons, e.σ_bg;
        σ_xy=e.σ_xy, frame=e.frame, dataset=ds, track_id=e.track_id, id=e.id
    )
end

function _with_dataset(e::Emitter2D{T}, ds::Int) where T
    Emitter2D{T}(e.x, e.y, e.photons, e.σ_x, e.σ_y, e.frame, ds, e.track_id)
end

function _with_dataset(e::GaussMLE.Emitter2DFitSigma{T}, ds::Int) where T
    GaussMLE.Emitter2DFitSigma{T}(
        e.x, e.y, e.photons, e.bg, e.σ,
        e.σ_x, e.σ_y, e.σ_xy, e.σ_photons, e.σ_bg, e.σ_σ,
        e.pvalue, e.frame, ds, e.track_id, e.id
    )
end

function _with_dataset(e::GaussMLE.Emitter2DFitSigmaXY{T}, ds::Int) where T
    GaussMLE.Emitter2DFitSigmaXY{T}(
        e.x, e.y, e.photons, e.bg, e.σx, e.σy,
        e.σ_x, e.σ_y, e.σ_xy, e.σ_photons, e.σ_bg, e.σ_σx, e.σ_σy,
        e.pvalue, e.frame, ds, e.track_id, e.id
    )
end

function _with_dataset(e::GaussMLE.Emitter2DFitGaussMLE{T}, ds::Int) where T
    GaussMLE.Emitter2DFitGaussMLE{T}(
        e.x, e.y, e.photons, e.bg,
        e.σ_x, e.σ_y, e.σ_xy, e.σ_photons, e.σ_bg,
        e.pvalue, e.frame, ds, e.track_id, e.id
    )
end

function _with_dataset(e::GaussMLE.Emitter3DFitGaussMLE{T}, ds::Int) where T
    GaussMLE.Emitter3DFitGaussMLE{T}(
        e.x, e.y, e.z, e.photons, e.bg,
        e.σ_x, e.σ_y, e.σ_z, e.σ_xy, e.σ_xz, e.σ_yz, e.σ_photons, e.σ_bg,
        e.pvalue, e.frame, ds, e.track_id, e.id
    )
end

# ============================================================
# Output helpers (used by all step save functions)
# ============================================================

"""
    step_outdir(outdir, step_number, cfg) -> Union{String, Nothing}

Compute output directory for a step: `outdir/02_filter/`.
Returns nothing if outdir is nothing.
"""
function step_outdir(outdir::Union{String,Nothing}, step_number::Int, cfg::SMLMData.AbstractSMLMConfig)
    outdir === nothing && return nothing
    joinpath(outdir, "$(lpad(step_number, 2, '0'))_$(step_name(cfg))")
end

"""Save step config to `config.toml` in the step output directory."""
function _save_config!(dir::String, cfg::SMLMData.AbstractSMLMConfig)
    filepath = joinpath(dir, "config.toml")
    open(filepath, "w") do io
        println(io, "# $(nameof(typeof(cfg)))")
        println(io, "type = \"$(nameof(typeof(cfg)))\"")
        _write_config_fields!(io, cfg)
    end
end

"""Check if a value is a config-like struct (has fields, not a primitive/collection)."""
_is_config_struct(v) = isstructtype(typeof(v)) && !(v isa Union{Number, String, Symbol, AbstractArray, AbstractDict, Tuple, SMLMData.AbstractCamera})

"""Write config fields to TOML. Nested structs become [section] blocks."""
function _write_config_fields!(io::IO, cfg; section::String="")
    for f in fieldnames(typeof(cfg))
        v = getfield(cfg, f)
        v isa SMLMData.AbstractCamera && continue
        v === nothing && continue
        key = section == "" ? string(f) : "$(section).$(f)"
        if _is_config_struct(v)
            # Nested config -> TOML section
            println(io, "\n[$f]")
            println(io, "type = \"$(nameof(typeof(v)))\"")
            _write_config_fields!(io, v; section=string(f))
        elseif v isa String
            println(io, "$f = \"$v\"")
        elseif v isa Symbol
            println(io, "$f = \"$v\"")
        else
            println(io, "$f = $v")
        end
    end
end

"""
    _save_info!(dir::String, info; section::String="")

Write upstream Info struct fields to `info.toml` in TOML format.

Writes scalar fields (numbers, bools, strings, symbols, tuples of scalars).
Skips complex fields (arrays, dicts, structs like BasicSMLD, models).

When `section` is empty, writes a fresh file with type header.
When `section` is provided, appends a `[section]` block.
"""
function _save_info!(dir::String, info; section::String="")
    filepath = joinpath(dir, "info.toml")
    open(filepath, section == "" ? "w" : "a") do io
        if section == ""
            println(io, "# Upstream package info")
            println(io, "type = \"$(nameof(typeof(info)))\"")
        else
            println(io, "\n[$section]")
        end
        for f in fieldnames(typeof(info))
            v = getfield(info, f)
            _write_info_field!(io, f, v)
        end
    end
end

"""Write a single field to info.toml, skipping complex types."""
function _write_info_field!(io::IO, name::Symbol, v::Number)
    # Bool <: Number, so this method also handles true/false (TOML-valid as-is).
    println(io, "$name = $v")
end
function _write_info_field!(io::IO, name::Symbol, v::String)
    println(io, "$name = \"$v\"")
end
function _write_info_field!(io::IO, name::Symbol, v::Symbol)
    println(io, "$name = \"$v\"")
end
function _write_info_field!(io::IO, name::Symbol, v::Nothing)
    println(io, "$name = \"nothing\"")
end
function _write_info_field!(io::IO, name::Symbol, v::Tuple)
    # Only write tuples of scalars
    if all(x -> x isa Union{Number, Bool, String, Symbol}, v)
        vals = join([x isa String || x isa Symbol ? "\"$x\"" : "$x" for x in v], ", ")
        println(io, "$name = [$vals]")
    end
    # Skip tuples containing complex types
end
function _write_info_field!(io::IO, ::Symbol, ::Any)
    # Skip: AbstractVector, AbstractArray, AbstractDict, complex structs
end

# ============================================================
# Step SMLD checkpointing (opt-in via checkpoint=true on analyze())
# ============================================================

"""
    _save_step_smld(dir, smld; filename, kwargs...)

Persist a step's output SMLD via JLD2 so downstream iteration (e.g., diagnostic
plots, parameter sweeps, BaGoL re-runs) can resume without re-running the
upstream pipeline. The SMLD is stored under the `smld` key:

```julia
data = JLD2.load("path/smld_corrected.jld2")
smld = data["smld"]   # full BasicSMLD with camera, n_frames, n_datasets
```

Extra named values are stored as additional top-level keys (e.g., pass
`drift_model=drift_model` to embed the drift model alongside the SMLD).

No-op if `dir` is nothing.
"""
function _save_step_smld(dir::Union{String,Nothing}, smld::BasicSMLD;
                          filename::String="smld.jld2",
                          kwargs...)
    dir === nothing && return nothing
    mkpath(dir)   # ensure the step dir exists (some steps gate their own mkpath behind verbosity)
    path = joinpath(dir, filename)
    JLD2.jldsave(path; smld=smld, kwargs...)
    return path
end

# ============================================================
# Pipeline cache helpers (inter-step data passing via filesystem)
# ============================================================

"""
    cache_dir(outdir) -> Union{String, Nothing}

Returns `joinpath(outdir, ".cache")` or nothing if outdir is nothing.
"""
cache_dir(outdir::Union{String,Nothing}) = outdir === nothing ? nothing : joinpath(outdir, ".cache")

"""
    save_cache(outdir, filename; kwargs...)

Save data to `outdir/.cache/filename` via JLD2. No-op if outdir is nothing.
"""
function save_cache(outdir::Union{String,Nothing}, filename::String; kwargs...)
    outdir === nothing && return nothing
    dir = cache_dir(outdir)
    mkpath(dir)
    path = joinpath(dir, filename)
    JLD2.jldsave(path; kwargs...)
    return path
end

"""
    load_cache(outdir, filename) -> Union{Dict, Nothing}

Load data from `outdir/.cache/filename` via JLD2. Returns nothing if missing or outdir is nothing.
"""
function load_cache(outdir::Union{String,Nothing}, filename::String)
    outdir === nothing && return nothing
    dir = cache_dir(outdir)
    dir === nothing && return nothing
    path = joinpath(dir, filename)
    isfile(path) || return nothing
    JLD2.load(path)
end
