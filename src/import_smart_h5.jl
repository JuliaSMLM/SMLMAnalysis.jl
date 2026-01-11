"""
Import functions for SMART microscope HDF5 data format.

The SMART microscope stores data in HDF5 files with structure:
- /Main/data: 3D array (width, height, frames) of UInt16 camera data
- /Main/camera: Empty group (metadata elsewhere)
- /Main/laser_XXX: Laser control data
- /Main/stage_XXX: Stage position data
"""

using HDF5
using SMLMData

export load_smart_h5, load_smart_h5_info

"""
    load_smart_h5_info(filepath::String)

Load metadata about a SMART microscope HDF5 file without reading the full dataset.

Returns a NamedTuple with:
- filepath: Full path to the file
- width, height, nframes: Image dimensions
- dtype: Data type of the images
- file_size_gb: Approximate file size in GB

# Example
```julia
info = load_smart_h5_info("data/gatta_ruler/2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5")
println("File contains ", info.nframes, " frames of size ", info.width, "x", info.height)
```
"""
function load_smart_h5_info(filepath::String)
    HDF5.h5open(filepath, "r") do file
        data = file["Main/data"]
        dims = size(data)

        return (
            filepath = filepath,
            width = dims[1],
            height = dims[2],
            nframes = dims[3],
            dtype = eltype(data),
            file_size_gb = filesize(filepath) / 1e9
        )
    end
end

"""
    load_smart_h5(filepath::String; frame_range=nothing)

Load image data from a SMART microscope HDF5 file.

# Arguments
- `filepath::String`: Path to the HDF5 file
- `frame_range`: Optional range of frames to load (e.g., 1:1000). If nothing, loads all frames.

# Returns
- `data::Array{UInt16, 3}`: Image data (width, height, frames)

# Example
```julia
# Load all frames
data = load_smart_h5("data/experiment.h5")

# Load first 1000 frames
data = load_smart_h5("data/experiment.h5", frame_range=1:1000)
```
"""
function load_smart_h5(filepath::String; frame_range=nothing)
    HDF5.h5open(filepath, "r") do file
        data = file["Main/data"]

        if frame_range === nothing
            return read(data)
        else
            return data[:, :, frame_range]
        end
    end
end

"""
    load_smart_h5_frame(filepath::String, frame::Int)

Load a single frame from a SMART microscope HDF5 file.

# Arguments
- `filepath::String`: Path to the HDF5 file
- `frame::Int`: Frame number to load (1-indexed)

# Returns
- `frame_data::Matrix{UInt16}`: Single frame image

# Example
```julia
frame1 = load_smart_h5_frame("data/experiment.h5", 1)
```
"""
function load_smart_h5_frame(filepath::String, frame::Int)
    HDF5.h5open(filepath, "r") do file
        data = file["Main/data"]
        return data[:, :, frame]
    end
end

"""
    smart_h5_to_array(filepath::String; max_frames=nothing)

Load SMART microscope data as a properly formatted array for SMLM processing.

This function transposes the data from (width, height, frames) to (height, width, frames)
to match standard image conventions where the first dimension is rows (y) and second is columns (x).

# Arguments
- `filepath::String`: Path to the HDF5 file
- `max_frames`: Optional maximum number of frames to load

# Returns
- `data::Array{UInt16, 3}`: Image data in (rows, cols, frames) format
- `info::NamedTuple`: File metadata

# Example
```julia
data, info = smart_h5_to_array("data/experiment.h5", max_frames=1000)
println("Loaded ", size(data, 3), " frames of size ", size(data, 1), "x", size(data, 2))
```
"""
function smart_h5_to_array(filepath::String; max_frames=nothing, verbose=false)
    info = load_smart_h5_info(filepath)

    frame_range = if max_frames === nothing
        nothing
    else
        1:min(max_frames, info.nframes)
    end

    t1 = @elapsed data = load_smart_h5(filepath, frame_range=frame_range)
    verbose && println("    HDF5 read: $(round(t1, digits=2))s")

    # Transpose from (width, height, frames) to (height, width, frames)
    # to match standard image conventions
    t2 = @elapsed data = permutedims(data, (2, 1, 3))
    verbose && println("    permutedims: $(round(t2, digits=2))s")

    return data, info
end
