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
    _save_box_overlay(dir, filename, images, roi_batch, box_colors; title_prefix="Frame")

Save frame overlay with colored boxes around ROIs.

# Arguments
- `dir`: Output directory
- `filename`: Output filename (e.g., "detection_overlay.png")
- `images`: 3D array of frame images
- `roi_batch`: ROIBatch with x_corners, y_corners, frame_indices
- `box_colors`: Vector of colors, one per ROI (same order as roi_batch)
- `title_prefix`: Prefix for frame titles (default "Frame")
"""
function _save_box_overlay(dir, filename, images, roi_batch, box_colors; title_prefix="Frame", frame_labels=nothing, suptitle=nothing)
    n_frames = size(images, 3)
    frame_indices = [round(Int, x) for x in range(1, n_frames, length=min(12, n_frames))]
    # Use provided frame_labels for display, or fall back to frame_indices
    display_labels = frame_labels !== nothing ? frame_labels : frame_indices

    # For SMLM data: high-contrast stretch with dark background
    # Clip at background level for dark base, aggressive upper clip for bright spots
    sample_frames = frame_indices[1:min(4, length(frame_indices))]
    sample_data = vec(images[:, :, sample_frames])
    bg_level = median(sample_data)
    pmin = Float64(bg_level)  # Clip at background for dark base
    pmax = Float64(quantile(sample_data, 0.995))  # Aggressive clip to brighten spots

    fig = Figure(size=_grid_figure_size(images))
    box_size = roi_batch.roi_size

    # Optional supertitle
    if suptitle !== nothing
        Label(fig[0, 1:4], suptitle, fontsize=11)
    end

    for (idx, frame_num) in enumerate(frame_indices)
        row = div(idx - 1, 4) + 1
        col = mod(idx - 1, 4) + 1
        display_frame = display_labels[idx]

        ax = Axis(fig[row, col], title="$title_prefix $display_frame", aspect=DataAspect(), yreversed=true)
        frame_data = images[:, :, frame_num]'
        heatmap!(ax, frame_data, colormap=:grays, colorrange=(pmin, pmax))

        frame_mask = roi_batch.frame_indices .== frame_num
        if any(frame_mask)
            det_x = roi_batch.x_corners[frame_mask]
            det_y = roi_batch.y_corners[frame_mask]
            colors = box_colors[frame_mask]
            for (x, y, c) in zip(det_x, det_y, colors)
                lines!(ax, [x, x+box_size, x+box_size, x, x],
                          [y, y, y+box_size, y+box_size, y],
                    color=c, linewidth=0.5)
            end
        end
        hidedecorations!(ax)
    end

    save(joinpath(dir, filename), fig)
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
    println(io, "$name = $v")
end
function _write_info_field!(io::IO, name::Symbol, v::Bool)
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
