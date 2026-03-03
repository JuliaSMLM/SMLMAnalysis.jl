# MIC (MATLAB Instrument Control) H5 format loader
# Structure varies:
#   Old format: Channel01/Zposition001/DataXXXX (Dataset directly)
#   New format: Channel01/Zposition001/DataXXXX/DataXXXX (Group with nested Dataset)
# With per-pixel calibration in Calibration/ group

using HDF5

# Helper to check if path exists in HDF5 file (works with nested paths)
_h5_exists(f, path) = try; f[path]; true; catch; false; end

# Helper to resolve actual dataset path (handles both old and new H5 formats)
function _resolve_data_path(f, dk::String)
    base_path = "Channel01/Zposition001/$dk"
    obj = f[base_path]
    if obj isa HDF5.Dataset
        # Old format: DataXXXX is the dataset directly
        return base_path
    elseif obj isa HDF5.Group
        # New format: DataXXXX is a group containing DataXXXX dataset
        nested_path = "$base_path/$dk"
        if _h5_exists(f, nested_path) && f[nested_path] isa HDF5.Dataset
            return nested_path
        end
    end
    error("Cannot find dataset for $dk")
end

"""
    load_mic_h5_info(filepath) -> NamedTuple

Get info about MIC H5 file without loading all data.

Returns NamedTuple with fields:
- height, width: image dimensions
- n_frames: total frames across all blocks
- n_blocks: number of data blocks (datasets)
- frames_per_block: Vector of frame counts per block
- has_calibration: whether calibration data exists
"""
function load_mic_h5_info(filepath::String)
    h5open(filepath, "r") do f
        # Find data blocks
        zpos = f["Channel01/Zposition001"]
        data_keys = sort([k for k in keys(zpos) if startswith(k, "Data")])

        # Count frames and get size
        n_frames = 0
        h, w = 0, 0
        frames_per_block = Int[]

        for dk in data_keys
            try
                data_path = _resolve_data_path(f, dk)
                sz = size(f[data_path])
                if h == 0
                    h, w = sz[1], sz[2]
                end
                n_frames += sz[3]
                push!(frames_per_block, sz[3])
            catch
                # Skip blocks that don't have valid data
                continue
            end
        end

        return (
            height = h,
            width = w,
            n_frames = n_frames,
            n_blocks = length(data_keys),
            frames_per_block = frames_per_block,
            has_calibration = _h5_exists(f, "Calibration"),
            file_size_gb = filesize(filepath) / 1e9
        )
    end
end

"""
    load_mic_h5_calibration(filepath) -> NamedTuple

Load calibration data from MIC H5 file.

Returns NamedTuple with (offset, variance, gain) as 2D arrays.

NOTE: The gain in MIC H5 files is INVERTED compared to our convention.
Our convention: ADU = photons * gain, so gain ≈ 0.24 e-/ADU
MIC H5: gain stored as ~4, which is 1/gain in our convention
Use load_mic_h5_calibration_for_scmos() to get corrected values.
"""
function load_mic_h5_calibration(filepath::String)
    h5open(filepath, "r") do f
        offset = read(f["Calibration/CCDOffset"])
        variance = read(f["Calibration/CCDVar"])
        gain = read(f["Calibration/Gain"])
        return (offset=offset, variance=variance, gain=gain)
    end
end

"""
    load_mic_h5_calibration_for_scmos(filepath) -> NamedTuple

Load calibration data and convert to SCMOSCamera convention.

Returns NamedTuple with:
- offset: per-pixel offset (unchanged)
- readnoise: per-pixel readnoise (sqrt of variance)
- gain: per-pixel gain (INVERTED from stored value)
"""
function load_mic_h5_calibration_for_scmos(filepath::String)
    cal = load_mic_h5_calibration(filepath)
    return (
        offset = Float32.(cal.offset),
        readnoise = Float32.(sqrt.(cal.variance)),
        gain = Float32.(1.0 ./ cal.gain)  # Invert gain to our convention
    )
end

"""
    build_camera_from_mic_h5(filepath; pixel_size, qe=1.0) -> SCMOSCamera

Build an SCMOSCamera from MIC H5 per-pixel calibration data (offset, readnoise, gain).

Pixel size and QE are not stored in MIC H5 files and must be provided.

# Arguments
- `filepath`: Path to MIC H5 file with Calibration/ group
- `pixel_size`: Pixel size in μm (required)
- `qe`: Quantum efficiency 0-1 (default: 1.0)
"""
function build_camera_from_mic_h5(filepath::String; pixel_size::Real, qe::Real=1.0)
    cal = load_mic_h5_calibration_for_scmos(filepath)
    ny, nx = size(cal.readnoise)
    SCMOSCamera(nx, ny, Float32(pixel_size), cal.readnoise;
                offset=cal.offset, gain=cal.gain, qe=Float32(qe))
end

"""
    load_mic_h5_block(filepath, block_num::Int) -> Array{Float32,3}

Load a single data block from MIC H5 file.
block_num is 1-indexed.
"""
function load_mic_h5_block(filepath::String, block_num::Int)
    h5open(filepath, "r") do f
        zpos = f["Channel01/Zposition001"]
        data_keys = sort([k for k in keys(zpos) if startswith(k, "Data")])

        if block_num < 1 || block_num > length(data_keys)
            error("Block $block_num out of range (1:$(length(data_keys)))")
        end

        dk = data_keys[block_num]
        data_path = _resolve_data_path(f, dk)
        Float32.(read(f[data_path]))
    end
end

"""
    load_mic_h5(filepath; max_frames=nothing, max_blocks=nothing) -> images, dataset_indices

Load MIC dSTORM H5 file. Each block becomes a "dataset".

Returns:
- images: 3D array (height × width × n_frames)
- dataset_indices: Vector{Int} mapping each frame to its block (1-indexed)
"""
function load_mic_h5(filepath::String;
                          max_frames::Union{Int,Nothing}=nothing,
                          max_blocks::Union{Int,Nothing}=nothing)
    h5open(filepath, "r") do f
        # Find all data blocks
        zpos = f["Channel01/Zposition001"]
        data_keys = sort([k for k in keys(zpos) if startswith(k, "Data")])

        # Limit blocks if requested
        if max_blocks !== nothing
            data_keys = data_keys[1:min(max_blocks, length(data_keys))]
        end

        # Count total frames (using resolved paths)
        frames_per_block = Int[]
        valid_keys = String[]
        for dk in data_keys
            try
                data_path = _resolve_data_path(f, dk)
                push!(frames_per_block, size(f[data_path], 3))
                push!(valid_keys, dk)
            catch
                continue
            end
        end
        data_keys = valid_keys
        n_frames_total = sum(frames_per_block)
        n_blocks = length(frames_per_block)

        # Limit frames if requested
        n_to_load = max_frames === nothing ? n_frames_total : min(max_frames, n_frames_total)

        # Get image size from first block
        first_data_path = _resolve_data_path(f, data_keys[1])
        first_data = f[first_data_path]
        h, w = size(first_data)[1:2]

        # Pre-allocate output
        images = Array{Float32}(undef, h, w, n_to_load)
        dataset_indices = Vector{Int}(undef, n_to_load)

        # Load data blocks
        frame_idx = 1
        for (block_num, dk) in enumerate(data_keys)
            data_path = _resolve_data_path(f, dk)
            block_data = read(f[data_path])
            n_block = size(block_data, 3)

            # How many frames to copy from this block
            frames_remaining = n_to_load - frame_idx + 1
            n_copy = min(n_block, frames_remaining)

            if n_copy <= 0
                break
            end

            images[:, :, frame_idx:frame_idx+n_copy-1] = Float32.(block_data[:, :, 1:n_copy])
            dataset_indices[frame_idx:frame_idx+n_copy-1] .= block_num
            frame_idx += n_copy

            if frame_idx > n_to_load
                break
            end
        end

        return images, dataset_indices
    end
end
