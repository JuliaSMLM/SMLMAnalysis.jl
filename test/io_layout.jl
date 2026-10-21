using SMLMAnalysis
using Test

# ============================================================
# Image layout (lab decision 0040): readers return [y, x] then frame, maps (ny, nx).
# All tests use non-square arrays (ny = 5, nx = 7) with per-pixel content 100y + x, so a
# swapped axis shows up in the shape and in the values.
# ============================================================

const _LAYOUT_NY, _LAYOUT_NX = 5, 7
function _layout_yx(nf = 3)
    return Float32[100y + x + 1000f for y in 1:_LAYOUT_NY, x in 1:_LAYOUT_NX, f in 1:nf]
end
_layout_xy(a) = permutedims(a, (2, 1, 3))   # what a legacy [x, y] writer stores

function _write_smart(path, data; tag = nothing)
    return SMLMAnalysis.HDF5.h5open(path, "w") do f
        g = SMLMAnalysis.HDF5.create_group(f, "Main")
        g["data"] = UInt16.(data)
        tag === nothing || (SMLMAnalysis.HDF5.attrs(g["data"])["dimension_order"] = tag)
    end
end

# blocks: vector of [y, x, frame] arrays; maps: (ny, nx) calibration maps or nothing.
# stored_xy: store the data transposed (a legacy MATLAB file); maps_xy: same for the maps
# (default: like the data). tag / map_tag: dimension_order of the data / of each map.
function _write_mic(
        path, blocks; maps = nothing, stored_xy = false, maps_xy = stored_xy,
        writer = nothing, tag = nothing, map_tag = nothing
    )
    H5 = SMLMAnalysis.HDF5
    st(a, xy) = xy ? (ndims(a) == 2 ? permutedims(a) : _layout_xy(a)) : a
    return H5.h5open(path, "w") do f
        md = H5.create_group(f, "Metadata")
        H5.attrs(md)["Version"] = "SEQv2"
        writer === nothing || (H5.attrs(md)["Writer"] = writer)
        for (i, b) in enumerate(blocks)
            dk = "Data" * lpad(i, 4, '0')
            g = H5.create_group(f, "Channel01/Zposition001/$dk")
            g[dk] = round.(UInt16, st(b, stored_xy))
            btag = tag isa AbstractVector ? tag[i] : tag
            btag === nothing || (H5.attrs(g[dk])["dimension_order"] = btag)
        end
        if maps !== nothing
            cg = H5.create_group(f, "Calibration")
            for (k, m) in (
                    ("CCDOffset", maps.offset), ("CCDVar", maps.variance),
                    ("Gain", maps.gain),
                )
                cg[k] = st(m, maps_xy)
                map_tag === nothing || (H5.attrs(cg[k])["dimension_order"] = map_tag)
            end
        end
    end
end

function _layout_maps()
    return (
        offset = Float32[100y + x for y in 1:_LAYOUT_NY, x in 1:_LAYOUT_NX],
        variance = fill(4.0f0, _LAYOUT_NY, _LAYOUT_NX),
        gain = fill(2.0f0, _LAYOUT_NY, _LAYOUT_NX),
    )
end

@testset "SMART legacy untagged [x, y] file loads as [y, x]" begin
    mktempdir() do dir
        path = joinpath(dir, "legacy.h5")
        yx = _layout_yx()
        _write_smart(path, _layout_xy(yx))             # stored (7, 5, 3), untagged
        @test SMLMAnalysis.HDF5.h5open(f -> size(f["Main/data"]), path) == (7, 5, 3)
        @test load_smart_h5(path) == UInt16.(yx)
        @test size(load_smart_h5(path)) == (5, 7, 3)
        data, info = smart_h5_to_array(path)
        @test data == UInt16.(yx)
        @test (info.height, info.width, info.nframes) == (5, 7, 3)
        @test load_smart_h5_info(path).height == 5
        @test load_smart_h5_info(path).width == 7
        @test load_smart_h5(path; frame_range = 2) == UInt16.(yx[:, :, 2])
        @test load_smart_h5(path; frame_range = 2:3) == UInt16.(yx[:, :, 2:3])
        @test smart_h5_to_array(path; max_frames = 2)[1] == UInt16.(yx[:, :, 1:2])
    end
end

@testset "SMART tagged files load as stored" begin
    mktempdir() do dir
        yx = _layout_yx()
        for tag in ("HWN", "y,x,frame")
            path = joinpath(dir, "tagged.h5")
            _write_smart(path, yx; tag = tag)         # stored (5, 7, 3)
            @test load_smart_h5(path) == UInt16.(yx)
            @test smart_h5_to_array(path)[1] == UInt16.(yx)
            info = load_smart_h5_info(path)
            @test (info.height, info.width, info.nframes) == (5, 7, 3)
            @test load_smart_h5(path; frame_range = 3) == UInt16.(yx[:, :, 3])
        end
    end
end

@testset "MIC legacy MATLAB [x, y] file loads as [y, x]" begin
    mktempdir() do dir
        path = joinpath(dir, "matlab.h5")
        b1, b2 = _layout_yx(3), _layout_yx(2) .+ 5000.0f0
        maps = _layout_maps()
        _write_mic(path, [b1, b2]; maps = maps, stored_xy = true)   # no Writer
        other = joinpath(dir, "other.h5")   # a Writer that is not MicroscopeSeqSR: legacy
        _write_mic(
            other, [b1, b2]; maps = maps, stored_xy = true, writer = "SomeOtherWriter 1.0"
        )
        @test load_mic_h5(other)[1] == load_mic_h5(path)[1]
        @test load_mic_h5_calibration(other).offset == maps.offset
        @test load_mic_h5_block(path, 1) == b1
        @test load_mic_h5_block(path, 2) == b2
        images, idx = load_mic_h5(path)
        @test size(images) == (5, 7, 5)
        @test images == cat(b1, b2; dims = 3)
        @test idx == [1, 1, 1, 2, 2]
        @test load_mic_h5(path; max_frames = 4)[1] == images[:, :, 1:4]
        info = load_mic_h5_info(path)
        @test (info.height, info.width, info.n_frames) == (5, 7, 5)
        cal = load_mic_h5_calibration(path)
        @test size(cal.offset) == (5, 7)
        @test cal.offset == maps.offset
        cam = build_camera_from_mic_h5(path; pixel_size = 0.1)
        @test size(cam.offset) == (5, 7)
        @test (length(cam.pixel_edges_x), length(cam.pixel_edges_y)) == (8, 6)
        # the map pixel [y, x] belongs to the data pixel [y, x] it was written with
        @test cam.offset[4, 6] == 100 * 4 + 6
        @test cam.offset[2, 7] == 100 * 2 + 7
    end
end

@testset "MIC from the SEQv2 writer and tagged MIC load as stored" begin
    mktempdir() do dir
        b = _layout_yx(3)
        maps = _layout_maps()
        cases = (
            (writer = "MicroscopeSeqSR 0.1.0-DEV", tag = nothing),
            (writer = nothing, tag = "y,x,frame"),
        )
        for c in cases
            path = joinpath(dir, "julia.h5")
            _write_mic(path, [b]; maps = maps, writer = c.writer, tag = c.tag)
            @test load_mic_h5_block(path, 1) == b
            @test load_mic_h5(path)[1] == b
            info = load_mic_h5_info(path)
            @test (info.height, info.width) == (5, 7)
            @test load_mic_h5_calibration(path).offset == maps.offset
            cam = build_camera_from_mic_h5(path; pixel_size = 0.1)
            @test size(cam.offset) == (5, 7)
            @test cam.offset[4, 6] == 406
        end
    end
end

@testset "unknown dimension_order throws ArgumentError" begin
    mktempdir() do dir
        smart = joinpath(dir, "s.h5")
        _write_smart(smart, _layout_yx(); tag = "XYN (dim1 = SLM x)")
        @test_throws ArgumentError load_smart_h5(smart)
        @test_throws ArgumentError load_smart_h5_info(smart)
        @test_throws ArgumentError smart_h5_to_array(smart)
        err = try
            load_smart_h5(smart)
        catch e
            e
        end
        @test occursin("XYN (dim1 = SLM x)", sprint(showerror, err))
        mic = joinpath(dir, "m.h5")
        _write_mic(mic, [_layout_yx()]; tag = "frame,y,x")
        @test_throws ArgumentError load_mic_h5(mic)
        @test_throws ArgumentError load_mic_h5_block(mic, 1)
    end
end

@testset "exact dimension_order spellings only" begin
    mktempdir() do dir
        yx, bad = _layout_yx(), ["x,yframe", "x,y,garbage", "x,y,z,frame", "y,x,", "hwn"]
        for tag in bad
            smart = joinpath(dir, "s.h5")
            _write_smart(smart, yx; tag = tag)
            @test_throws ArgumentError load_smart_h5(smart)
            @test_throws ArgumentError smart_h5_to_array(smart)
            mic = joinpath(dir, "m.h5")
            _write_mic(mic, [yx]; tag = tag)
            @test_throws ArgumentError load_mic_h5(mic)
        end
        for tag in ("x,y,frame", "x,y,z")   # exact [x, y] spellings convert
            smart = joinpath(dir, "sx.h5")
            _write_smart(smart, _layout_xy(yx); tag = tag)
            @test load_smart_h5(smart) == UInt16.(yx)
            mic = joinpath(dir, "mx.h5")
            _write_mic(mic, [yx]; stored_xy = true, tag = tag)
            @test load_mic_h5(mic)[1] == yx
        end
        for tag in ("y,x", "y,x,z", "HW")   # exact [y, x] spellings load as stored
            smart = joinpath(dir, "sy.h5")
            _write_smart(smart, yx; tag = tag)
            @test load_smart_h5(smart) == UInt16.(yx)
        end
    end
end

@testset "unknown tag in any MIC block throws, never a zero-size info" begin
    mktempdir() do dir
        b = _layout_yx()
        for tags in (["frame,y,x", "y,x,frame"], ["y,x,frame", "frame,y,x"])
            path = joinpath(dir, "blocks.h5")
            _write_mic(path, [b, b]; tag = tags)
            @test_throws ArgumentError load_mic_h5_info(path)
            @test_throws ArgumentError load_mic_h5(path)
            @test_throws ArgumentError load_mic_h5_block(path, 1)
            @test_throws ArgumentError load_mic_h5_block(path, 2)
        end
    end
end

@testset "MIC maps follow the data layout, own tag wins, mismatch throws" begin
    mktempdir() do dir
        b, maps = _layout_yx(), _layout_maps()
        # tagged [y, x] data, untagged maps, no Writer: maps take the data's layout
        p1 = joinpath(dir, "tagged_data.h5")
        _write_mic(p1, [b]; maps = maps, tag = "y,x,frame")
        @test load_mic_h5(p1)[1] == b
        @test load_mic_h5_calibration(p1).offset == maps.offset
        @test size(build_camera_from_mic_h5(p1; pixel_size = 0.1).offset) == (5, 7)
        # tagged [y, x] maps with untagged legacy [x, y] data
        p2 = joinpath(dir, "tagged_maps.h5")
        _write_mic(p2, [b]; maps = maps, stored_xy = true, maps_xy = false, map_tag = "y,x")
        @test load_mic_h5(p2)[1] == b
        @test load_mic_h5_calibration(p2).offset == maps.offset
        # untagged maps in a Writer file follow the Writer rule too
        p3 = joinpath(dir, "writer.h5")
        _write_mic(p3, [b]; maps = maps, writer = "MicroscopeSeqSR 0.1.0-DEV")
        @test load_mic_h5_calibration(p3).offset == maps.offset
        # legacy data with maps whose size does not match: ArgumentError naming both
        p4 = joinpath(dir, "mismatch.h5")
        bad = (
            offset = fill(1.0f0, 6, 7), variance = fill(1.0f0, 6, 7),
            gain = fill(1.0f0, 6, 7),
        )
        _write_mic(
            p4, [b]; maps = bad, stored_xy = true, maps_xy = false, map_tag = "y,x"
        )
        err = try
            load_mic_h5_calibration(p4)
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("(6, 7)", err.msg) && occursin("(5, 7)", err.msg)
        @test_throws ArgumentError build_camera_from_mic_h5(p4; pixel_size = 0.1)
    end
end

@testset "max_blocks does not skip the layout check of later blocks" begin
    mktempdir() do dir
        b = _layout_yx()
        path = joinpath(dir, "later.h5")
        _write_mic(path, [b, b]; tag = ["y,x,frame", "frame,y,x"])
        @test_throws ArgumentError load_mic_h5(path; max_blocks = 1)
        @test_throws ArgumentError load_mic_h5(path)
    end
end
