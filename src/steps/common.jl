"""
Common helper functions shared across analysis steps.
"""

"""
    _calculate_mode(values; n_bins=100)

Calculate the mode of a distribution using histogram binning.
Returns the center of the most populated bin.
"""
function _calculate_mode(values::Vector{T}; n_bins=100) where T<:Real
    isempty(values) && return zero(T)

    valid = filter(x -> isfinite(x) && x > 0, values)
    isempty(valid) && return zero(T)

    lo, hi = quantile(valid, [0.01, 0.99])
    lo >= hi && return median(valid)

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
function _save_box_overlay(dir, filename, images, roi_batch, box_colors; title_prefix="Frame")
    n_frames = size(images, 3)
    frame_indices = [round(Int, x) for x in range(1, n_frames, length=min(12, n_frames))]

    pmin = Float64(quantile(vec(images[:,:,1]), 0.01))
    pmax = Float64(quantile(vec(images[:,:,1]), 0.99))

    fig = Figure(size=_grid_figure_size(images))
    box_size = roi_batch.roi_size

    for (idx, frame_num) in enumerate(frame_indices)
        row = div(idx - 1, 4) + 1
        col = mod(idx - 1, 4) + 1

        ax = Axis(fig[row, col], title="$title_prefix $frame_num", aspect=DataAspect(), yreversed=true)
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
