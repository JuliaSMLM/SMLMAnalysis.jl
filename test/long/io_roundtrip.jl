using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Random
using TOML
using Statistics

@testset "SMLD HDF5 round-trip" begin
    # Locks the σ_xy regression: save_smld/load_smld must preserve every
    # emitter field, including the position covariance σ_xy that the
    # GaussianXYNBS → Emitter2DFitSigma path carries. This bug survived
    # because the prior tests only constructed types, never round-tripped.
    cam = IdealCamera(8, 8, 0.1)
    T = Float64

    mktempdir() do dir
        # Emitter2DFitSigma (16 fields) — the primary GaussianXYNBS output.
        es = [
            GaussMLE.Emitter2DFitSigma{T}(
                0.1i, 0.2i, 1000.0 + i, 5.0, 0.13,     # x, y, photons, bg, σ
                0.01, 0.012, 0.003, 20.0, 0.5, 0.002,  # σ_x, σ_y, σ_xy, σ_photons, σ_bg, σ_σ
                0.4, i, 1, 0, i
            )                        # pvalue, frame, dataset, track_id, id
                for i in 1:5
        ]
        smld = SMLMAnalysis.BasicSMLD(es, cam, 10, 1, Dict{String, Any}())
        path = joinpath(dir, "sigma.h5")
        save_smld(path, smld)
        loaded = load_smld(path)

        @test length(loaded.emitters) == 5
        @test eltype(loaded.emitters) <: GaussMLE.Emitter2DFitSigma
        for (a, b) in zip(smld.emitters, loaded.emitters)
            for f in fieldnames(GaussMLE.Emitter2DFitSigma)
                @test getfield(a, f) ≈ getfield(b, f)
            end
        end
        @test loaded.emitters[3].σ_xy ≈ 0.003   # the field that used to vanish

        # Emitter2DFitSigmaXY (18 fields) — GaussianXYNBSXSY output.
        exy = [
            GaussMLE.Emitter2DFitSigmaXY{T}(
                0.1i, 0.2i, 1000.0 + i, 5.0, 0.13, 0.14, # x, y, photons, bg, σx, σy
                0.01, 0.012, 0.003, 20.0, 0.5,           # σ_x, σ_y, σ_xy, σ_photons, σ_bg
                0.002, 0.0021, 0.4, i, 1, 0, i
            )          # σ_σx, σ_σy, pvalue, frame, dataset, track_id, id
                for i in 1:4
        ]
        smld_xy = SMLMAnalysis.BasicSMLD(exy, cam, 10, 1, Dict{String, Any}())
        pxy = joinpath(dir, "sigmaxy.h5")
        save_smld(pxy, smld_xy)
        loaded_xy = load_smld(pxy)

        @test length(loaded_xy.emitters) == 4
        @test eltype(loaded_xy.emitters) <: GaussMLE.Emitter2DFitSigmaXY
        for (a, b) in zip(smld_xy.emitters, loaded_xy.emitters)
            for f in fieldnames(GaussMLE.Emitter2DFitSigmaXY)
                @test getfield(a, f) ≈ getfield(b, f)
            end
        end
        @test loaded_xy.emitters[2].σ_xy ≈ 0.003

        # Schema validation: a valid HDF5 file that isn't an SMLD fails with a
        # friendly ErrorException, not a raw KeyError deep in the read.
        bogus = joinpath(dir, "bogus.h5")
        SMLMAnalysis.HDF5.h5open(bogus, "w") do f
            f["junk"] = 1
        end
        @test_throws ErrorException load_smld(bogus)
    end
end

@testset "step checkpoint is versioned HDF5" begin
    # _save_step_smld writes through save_smld; load_smld must read it back unchanged.
    cam = IdealCamera(16, 16, 0.1)
    em = [
        SMLMAnalysis.Emitter2DFit{Float64}(
            0.1i, 0.2i, 1000.0 + i, 5.0, 0.01, 0.012, 20.0, 0.5;
            frame = i, dataset = 1 + (i % 2)
        ) for i in 1:6
    ]
    smld = SMLMAnalysis.BasicSMLD(em, cam, 6, 2, Dict{String, Any}())
    dm = SMLMDriftCorrection.LegendrePolynomial(smld; degree = 2)
    mktempdir() do dir
        p = SMLMAnalysis._save_step_smld(
            joinpath(dir, "03_driftcorrect"), smld;
            filename = "smld_corrected.h5", drift_model = dm
        )
        @test p == joinpath(dir, "03_driftcorrect", "smld_corrected.h5") && isfile(p)
        s2 = load_smld(p)
        @test s2.emitters isa Vector{SMLMAnalysis.Emitter2DFit{Float64}}
        @test all(
            getfield(a, f) == getfield(b, f) for (a, b) in zip(em, s2.emitters)
                for f in fieldnames(SMLMAnalysis.Emitter2DFit{Float64})
        )
        @test (s2.n_frames, s2.n_datasets) == (6, 2)
        @test s2.camera.pixel_edges_x == cam.pixel_edges_x
        @test s2.metadata["drift_correction"]["model_type"] == "LegendrePolynomial"
        @test SMLMAnalysis._save_step_smld(nothing, smld; filename = "x.h5") === nothing
    end
end

@testset "TOML provenance is valid" begin
    # Upstream field names outside TOML's bare-key set (DriftConfig.σ_loc) must be
    # written as quoted keys, or config.toml fails to parse.
    let io = IOBuffer()
        SMLMAnalysis._write_config_fields!(io, DriftConfig())
        parsed = TOML.parse(String(take!(io)))
        @test haskey(parsed, "σ_loc")
    end
    # Provenance files are named .toml and must parse back. The hand-rolled
    # serializer used to emit invalid TOML for tuples ((500.0, Inf)), ranges
    # (1:19), symbol vectors ([:red, :blue]), and unescaped strings. _toml_value
    # fixes every scalar-emit point; assert the round-trip actually works.
    mktempdir() do dir
        # FilterConfig: photons/precision are Tuple{Float64,Float64}; Inf must
        # serialize as TOML `inf`. FilterConfig has no nested-config fields, so
        # every value sits at the root table.
        cfg = FilterConfig(photons = (500.0, Inf), precision = (0.0, 0.02))
        SMLMAnalysis._save_config!(dir, cfg)
        parsed = TOML.parsefile(joinpath(dir, "config.toml"))   # throws if invalid TOML
        @test parsed["type"] == "FilterConfig"
        @test parsed["photons"][1] == 500.0
        @test isinf(parsed["photons"][2])        # `inf` parses back to Inf::Float64
        @test parsed["precision"] == [0.0, 0.02]

        # DetectFitConfig.datasets is an AbstractVector{Int}; a UnitRange (1:19)
        # used to emit the bare, invalid `datasets = 1:19`. It must now parse.
        cfg2 = DetectFitConfig(datasets = 1:19)
        SMLMAnalysis._save_config!(dir, cfg2)
        @test TOML.parsefile(joinpath(dir, "config.toml")) isa AbstractDict  # no parse error

        # _toml_value unit behavior for the hard scalar types.
        @test SMLMAnalysis._toml_value(1:19) == "\"1:19\""       # range -> quoted, not materialized
        @test SMLMAnalysis._toml_value((500.0, Inf)) == "[500.0, inf]"
        @test SMLMAnalysis._toml_value([:red, :blue]) == "[\"red\", \"blue\"]"
        @test SMLMAnalysis._toml_value(true) == "true"

        # String escaping: embedded quote, backslash, and newline survive a
        # write -> TOML.parse round-trip (the old String branch escaped nothing).
        s = "he said \"hi\"\\\n done"
        open(joinpath(dir, "esc.toml"), "w") do io
            println(io, "s = $(SMLMAnalysis._toml_value(s))")
        end
        @test TOML.parsefile(joinpath(dir, "esc.toml"))["s"] == s

        # info.toml: a `nothing` field is OMITTED (as in config.toml), never written
        # as the string "nothing" — otherwise a key's TOML type flips between runs.
        info = (; n_after = 3, p2_estimate = nothing, field_mode = :gaussian)
        SMLMAnalysis._save_info!(dir, info)
        parsed_info = TOML.parsefile(joinpath(dir, "info.toml"))
        @test parsed_info["n_after"] == 3
        @test parsed_info["field_mode"] == "gaussian"
        @test !haskey(parsed_info, "p2_estimate")
    end
end

@testset "MIC H5 calibration units" begin
    # SMITE convention: RawData = Gain_stored*photons + Offset (ADU); CCDVar is the
    # dark variance in ADU². Regression coverage for the ADU→e⁻ readnoise conversion
    # (load_mic_h5_calibration_for_scmos must divide by the stored ADU/e⁻ gain, not
    # just take sqrt(variance)).
    mktempdir() do dir
        path = joinpath(dir, "cal.h5")
        offset_adu = fill(100.0, 4, 4)
        gain_stored = fill(2.0, 4, 4)  # ADU/e⁻
        var_adu2 = fill(4.0, 4, 4)     # ADU² → readnoise = sqrt(4)/2 = 1.0 e⁻ rms
        SMLMAnalysis.HDF5.h5open(path, "w") do f
            g = SMLMAnalysis.HDF5.create_group(f, "Calibration")
            g["CCDOffset"] = offset_adu
            g["CCDVar"] = var_adu2
            g["Gain"] = gain_stored
        end

        cal = SMLMAnalysis.load_mic_h5_calibration_for_scmos(path)
        @test all(cal.readnoise .≈ 1.0f0)
        @test all(cal.gain .≈ 0.5f0)     # e⁻/ADU = 1/gain_stored
        @test all(cal.offset .≈ 100.0f0)

        cam = build_camera_from_mic_h5(path; pixel_size = 0.1)
        @test all(cam.readnoise .≈ 1.0f0)
        @test all(cam.gain .≈ 0.5f0)
        @test all(cam.offset .≈ 100.0f0)
    end
end

@testset "load_mic_h5 warns on an unreadable data block" begin
    # Data002 is a group without its nested dataset, so _resolve_data_path
    # throws for it — load_mic_h5 must @warn (naming the block and the error)
    # and skip it, rather than the previous bare `catch; continue`.
    mktempdir() do dir
        path = joinpath(dir, "corrupt.h5")
        SMLMAnalysis.HDF5.h5open(path, "w") do f
            g = SMLMAnalysis.HDF5.create_group(f, "Channel01/Zposition001")
            g["Data001"] = rand(Float32, 4, 4, 3)
            SMLMAnalysis.HDF5.create_group(g, "Data002")
        end
        @test_logs (:warn, r"load_mic_h5: skipping unreadable data block \"Data002\"") match_mode = :any load_mic_h5(path)
        images, dataset_indices = load_mic_h5(path)
        @test size(images, 3) == 3
        @test all(==(1), dataset_indices)
    end
end

@testset "atomic saves refuse a directory destination" begin
    # save_smld goes through _replace_atomically: temp file + rename(2). rename
    # fails on a directory destination rather than deleting it (and everything
    # inside it), so it must throw and leave the directory intact. The exception
    # type differs by Julia version, so only a throw is asserted.
    cam = IdealCamera(8, 8, 0.1)
    T = Float64
    es = [
        SMLMAnalysis.Emitter2DFit{T}(0.1i, 0.2i, 1000.0 + i, 5.0, 0.01, 0.012, 20.0, 0.5, 0.4, i, 1, 0, i)
            for i in 1:3
    ]
    smld = SMLMAnalysis.BasicSMLD(es, cam, 10, 1, Dict{String, Any}())

    mktempdir() do dir
        # save_smld: destination is a directory, not a file.
        target = joinpath(dir, "out.h5")
        mkdir(target)
        inner = joinpath(target, "keepme.txt")
        write(inner, "do not delete me")

        @test_throws Exception save_smld(target, smld)
        @test isdir(target)
        @test isfile(inner)
        @test read(inner, String) == "do not delete me"
        @test isempty(filter(f -> f != "keepme.txt", readdir(target)))

        # A normal save (no pre-existing directory) still works and leaves no
        # stray temp files behind.
        good = joinpath(dir, "good.h5")
        save_smld(good, smld)
        @test isfile(good)
        @test isempty(filter(f -> f ∉ ("out.h5", "good.h5"), readdir(dir)))
    end

    # Race: a directory appears at the destination after the isdir pre-check.
    # rename(2) itself must refuse it (never a recursive delete), and the temp
    # must be cleaned up. This is the property the pre-check cannot prove.
    mktempdir() do dir
        p = joinpath(dir, "out.h5")
        @test_throws Exception SMLMAnalysis._replace_atomically(p) do tmp
            write(tmp, "x")
            mkdir(p)
            write(joinpath(p, "keep"), "k")
        end
        @test read(joinpath(p, "keep"), String) == "k"
        @test readdir(dir) == ["out.h5"]
    end

    # The writer's temp path lives in a directory only we can enter (POSIX), so
    # no other user sharing the destination directory can pre-place a symlink at it.
    if !Sys.iswindows()
        mktempdir() do dir
            p = joinpath(dir, "out.txt")
            SMLMAnalysis._replace_atomically(p) do tmp
                @test dirname(tmp) != dir
                @test dirname(dirname(tmp)) == dir
                @test filemode(dirname(tmp)) & 0o777 == 0o700
                @test !ispath(tmp)
                write(tmp, "x")
            end
            @test read(p, String) == "x"
            @test readdir(dir) == ["out.txt"]
        end
    end

    # A failing writer leaves no temp behind and never touches the destination.
    mktempdir() do dir
        p = joinpath(dir, "out.h5")
        @test_throws ErrorException SMLMAnalysis._replace_atomically(tmp -> error("boom"), p)
        @test !ispath(p)
        @test isempty(readdir(dir))

        write(p, "old content")
        @test_throws ErrorException SMLMAnalysis._replace_atomically(p) do tmp
            write(tmp, "partial")
            error("boom")
        end
        @test read(p, String) == "old content"
        @test readdir(dir) == ["out.h5"]
    end

    if !Sys.iswindows()
        # A symlink destination is replaced as a directory entry, not written through.
        mktempdir() do dir
            target = joinpath(dir, "target.h5")
            write(target, "target content")
            link = joinpath(dir, "link.h5")
            symlink(target, link)
            save_smld(link, smld)
            @test !islink(link)
            @test isfile(link)
            @test length(load_smld(link).emitters) == 3
            @test read(target, String) == "target content"
        end

        # The saved file's mode follows the umask, as a direct create would
        # (not the 0600 a mktemp-created temp would carry over).
        mktempdir() do dir
            p = joinpath(dir, "mode.h5")
            mask = UInt32(0o022)
            old = ccall(:umask, UInt32, (UInt32,), mask)
            try
                save_smld(p, smld)
            finally
                ccall(:umask, UInt32, (UInt32,), old)
            end
            @test filemode(p) & 0o777 == 0o666 & ~mask
        end
    end
end

@testset "load_smart_h5 with an integer frame_range returns one frame" begin
    # `data[:, :, frame_range]` returns a Matrix (not a 3-D array) for an Int
    # frame_range -- a single frame index, not just a range, is a supported
    # (and documented) call shape.
    mktempdir() do dir
        path = joinpath(dir, "smart.h5")
        data = rand(UInt16, 8, 6, 3)
        SMLMAnalysis.HDF5.h5open(path, "w") do f
            g = SMLMAnalysis.HDF5.create_group(f, "Main")
            g["data"] = data
            # tagged [y, x]: an untagged file would be legacy [x, y] and convert
            SMLMAnalysis.HDF5.attrs(g["data"])["dimension_order"] = "HWN"
        end

        frame = load_smart_h5(path; frame_range = 2)
        @test frame isa Matrix
        @test frame == data[:, :, 2]

        frames = load_smart_h5(path; frame_range = 2:3)
        @test size(frames) == (8, 6, 2)
    end
end

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
# stored_xy: store data and maps transposed (a legacy MATLAB file).
function _write_mic(
        path, blocks; maps = nothing, stored_xy = false, writer = nothing, tag = nothing
    )
    H5 = SMLMAnalysis.HDF5
    st(a) = stored_xy ? (ndims(a) == 2 ? permutedims(a) : _layout_xy(a)) : a
    return H5.h5open(path, "w") do f
        md = H5.create_group(f, "Metadata")
        H5.attrs(md)["Version"] = "SEQv2"
        writer === nothing || (H5.attrs(md)["Writer"] = writer)
        for (i, b) in enumerate(blocks)
            dk = "Data" * lpad(i, 4, '0')
            g = H5.create_group(f, "Channel01/Zposition001/$dk")
            g[dk] = round.(UInt16, st(b))
            tag === nothing || (H5.attrs(g[dk])["dimension_order"] = tag)
        end
        if maps !== nothing
            cg = H5.create_group(f, "Calibration")
            cg["CCDOffset"] = st(maps.offset)
            cg["CCDVar"] = st(maps.variance)
            cg["Gain"] = st(maps.gain)
            if tag !== nothing
                for k in ("CCDOffset", "CCDVar", "Gain")
                    H5.attrs(cg[k])["dimension_order"] = "y,x"
                end
            end
        end
    end
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
        maps = (
            offset = Float32[100y + x for y in 1:5, x in 1:7],
            variance = fill(4.0f0, 5, 7),
            gain = fill(2.0f0, 5, 7),
        )
        _write_mic(path, [b1, b2]; maps = maps, stored_xy = true)   # no Writer
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
        @test (cam.pixel_edges_x |> length, cam.pixel_edges_y |> length) == (8, 6)
        # the map pixel [y, x] belongs to the data pixel [y, x] it was written with
        @test cam.offset[4, 6] == 100 * 4 + 6
        @test cam.offset[2, 7] == 100 * 2 + 7
    end
end

@testset "MIC from the SEQv2 writer and tagged MIC load as stored" begin
    mktempdir() do dir
        b = _layout_yx(3)
        maps = (
            offset = Float32[100y + x for y in 1:5, x in 1:7],
            variance = fill(4.0f0, 5, 7),
            gain = fill(2.0f0, 5, 7),
        )
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

@testset "detectfit from a legacy MIC file and a SEQv2-writer file agree" begin
    # One bright Gaussian spot at pixel (row r, col c) of a non-square frame, written once
    # as legacy MATLAB [x, y] and once as the Julia writer's [y, x]. Both must localize at
    # x = (c - 0.5) Δ, y = (r - 0.5) Δ (pixel [1, 1] center at (0.5Δ, 0.5Δ), pixels.md).
    ny, nx, nf, Δ = 24, 40, 6, 0.1
    r, c = 8, 29
    σpx = 0.13 / Δ
    frame = Float32[
        100 + 3000 * exp(-((y - r)^2 + (x - c)^2) / (2σpx^2)) for y in 1:ny, x in 1:nx
    ]
    block = repeat(frame; outer = (1, 1, nf))
    maps = (
        offset = fill(100.0f0, ny, nx),
        variance = fill(4.0f0, ny, nx),
        gain = fill(1.0f0, ny, nx),
    )
    mktempdir() do dir
        locs = map(
            (
                (stored_xy = true, writer = nothing),
                (stored_xy = false, writer = "MicroscopeSeqSR 0.1.0-DEV"),
            )
        ) do kw
            path = joinpath(dir, "spot_$(kw.stored_xy).h5")
            _write_mic(path, [block]; maps = maps, kw...)
            cfg = DetectFitConfig(
                boxer = BoxerConfig(boxsize = 7, psf_sigma = 0.13, backend = :cpu),
                fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu),
                path = path, pixel_size = Δ,
            )
            smld, _ = analyze(nothing, cfg; verbose = Verbosity.SILENT)
            @test length(smld.emitters) >= nf
            smld.emitters
        end
        for es in locs
            @test all(abs(e.x - (c - 0.5) * Δ) < 0.02 for e in es)
            @test all(abs(e.y - (r - 0.5) * Δ) < 0.02 for e in es)
        end
        @test length(locs[1]) == length(locs[2])
        @test locs[1][1].x ≈ locs[2][1].x atol = 1.0e-4
        @test locs[1][1].y ≈ locs[2][1].y atol = 1.0e-4
    end
end
