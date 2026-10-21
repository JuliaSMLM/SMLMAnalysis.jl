"""
Import functions for SMART microscope HDF5 data format.

The SMART microscope stores data in HDF5 files with structure:
- /Main/data: 3D array of UInt16 camera data, `(ny, nx, frames)` when the dataset carries
  a `dimension_order` attribute, legacy `(nx, ny, frames)` otherwise
- /Main/camera: Empty group (metadata elsewhere)
- /Main/laser_XXX: Laser control data
- /Main/stage_XXX: Stage position data

Image layout (lab decision 0040, SMLMDocs `pixels.md`
<https://JuliaSMLM.github.io/SMLMDocs.jl/dev/pixels/>): every loader returns images
as `[y, x]` then frame. A dataset's stored layout comes from its `dimension_order`
attribute; files without it (SMART files written before 2026-01-29) are legacy `[x, y]`
and are converted once, at load.
"""

# Stored layout of an image dataset from its `dimension_order` attribute, exact spellings
# only: `:yx` for "y,x,frame", "y,x", "y,x,z", "HWN", "HW"; `:xy` for "x,y,frame", "x,y",
# "x,y,z". With no attribute, `default` (the legacy rule of the caller's file format)
# applies.
function _h5_layout(obj, default::Symbol)
    a = HDF5.attrs(obj)
    haskey(a, "dimension_order") || return default
    tag = string(a["dimension_order"])
    tag in ("y,x,frame", "y,x", "y,x,z", "HWN", "HW") && return :yx
    tag in ("x,y,frame", "x,y", "x,y,z") && return :xy
    throw(ArgumentError("unrecognized dimension_order \"$tag\" on $(HDF5.name(obj))"))
end

# Legacy [x, y(, frame)] -> [y, x(, frame)]; the one place an image array is permuted.
_to_yx(a::AbstractArray) = ndims(a) == 2 ? permutedims(a) : permutedims(a, (2, 1, 3))

"""
    load_smart_h5_info(filepath::String)

Load metadata about a SMART microscope HDF5 file without reading the full dataset.

Returns a NamedTuple with:
- filepath: Full path to the file
- width, height, nframes: Image dimensions (`height` = rows = y, `width` = columns = x,
  whichever layout the file stores; see [`load_smart_h5`](@ref))
- dtype: Data type of the images
- file_size_gb: Approximate file size in GB

# Example
```julia
info = load_smart_h5_info("path/to/smart_data.h5")
println("File contains ", info.nframes, " frames of size ", info.width, "x", info.height)
```
"""
function load_smart_h5_info(filepath::String)
    return HDF5.h5open(filepath, "r") do file
        data = file["Main/data"]
        dims = size(data)
        layout = _h5_layout(data, :xy)
        ny, nx = layout == :yx ? (dims[1], dims[2]) : (dims[2], dims[1])

        return (
            filepath = filepath,
            width = nx,
            height = ny,
            nframes = dims[3],
            dtype = eltype(data),
            file_size_gb = filesize(filepath) / 1.0e9,
        )
    end
end

"""
    load_smart_h5(filepath::String; frame_range=nothing)

Load image data from a SMART microscope HDF5 file, as `[y, x]` then frame.

The stored layout is read from the dataset's `dimension_order` attribute (`"y,x,frame"`,
`"y,x"`, `"y,x,z"`, `"HWN"`, `"HW"` load as stored; `"x,y,frame"`, `"x,y"`, `"x,y,z"` are
converted; any other value throws `ArgumentError`). A file without the attribute is legacy
`[x, y]` (written before 2026-01-29) and is converted once here. See lab decision 0040 and
[SMLMDocs `pixels.md`](https://JuliaSMLM.github.io/SMLMDocs.jl/dev/pixels/).

# Arguments
- `filepath::String`: Path to the HDF5 file
- `frame_range`: Optional range of frames to load (e.g., 1:1000), or a single frame
  index (e.g., 5), which returns that one frame as a `Matrix`. If nothing, loads
  all frames.

# Returns
- `data::Array{UInt16, 3}`: Image data `(ny, nx, frames)` for a range or `nothing`;
  `Matrix{UInt16}` `(ny, nx)` for a single frame index.

# Example
```julia
# Load all frames
data = load_smart_h5("data/experiment.h5")

# Load first 1000 frames
data = load_smart_h5("data/experiment.h5", frame_range=1:1000)

# Load a single frame
frame5 = load_smart_h5("data/experiment.h5", frame_range=5)
```
"""
function load_smart_h5(filepath::String; frame_range = nothing)
    return HDF5.h5open(filepath, "r") do file
        data = file["Main/data"]
        layout = _h5_layout(data, :xy)

        raw = frame_range === nothing ? read(data) : data[:, :, frame_range]
        return layout == :yx ? raw : _to_yx(raw)
    end
end

"""
    smart_h5_to_array(filepath::String; max_frames=nothing)

Load SMART microscope data as an `[y, x]` array for SMLM processing.

Same layout rule as [`load_smart_h5`](@ref): `[y, x]` tags load unchanged, `[x, y]` tags
and untagged legacy files are converted once. Nothing is permuted twice.

# Arguments
- `filepath::String`: Path to the HDF5 file
- `max_frames`: Optional maximum number of frames to load

# Returns
- `data::Array{UInt16, 3}`: Image data `(ny, nx, frames)`
- `info::NamedTuple`: File metadata

# Example
```julia
data, info = smart_h5_to_array("data/experiment.h5", max_frames=1000)
println("Loaded ", size(data, 3), " frames of size ", size(data, 1), "x", size(data, 2))
```
"""
function smart_h5_to_array(filepath::String; max_frames = nothing, verbose = false)
    info = load_smart_h5_info(filepath)

    frame_range = if max_frames === nothing
        nothing
    else
        1:min(max_frames, info.nframes)
    end

    t1 = @elapsed data = load_smart_h5(filepath, frame_range = frame_range)
    verbose && println("    HDF5 read: $(round(t1, digits = 2))s")

    return data, info
end
