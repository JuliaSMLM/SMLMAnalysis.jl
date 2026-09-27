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

@testset "GaussMLE emitter types + abstract-eltype round-trip" begin
    # Regression coverage for the AbstractEmitter type-erasure + serialization fix:
    #  - detectfit sets each emitter's dataset field by mutating it in place
    #    (all advertised emitter types are mutable structs, so this needs no
    #    per-type dispatch — unlike the old positional-reconstruction
    #    `_with_dataset`, removed once every type was confirmed mutable), for
    #    EVERY advertised GaussMLE emitter type (GaussianXYNB → Emitter2DFitGaussMLE,
    #    AstigmaticXYZNB → Emitter3DFitGaussMLE) as well as the standard ones.
    #  - save_smld/load_smld must round-trip those types and the standard-3D
    #    off-diagonal covariances, AND must key off the concrete emitter even when
    #    the SMLD is typed SMLMAnalysis.BasicSMLD{T,AbstractEmitter} (as the pre-narrowing
    #    pipeline produced) rather than degrading to SMLMAnalysis.Emitter2DFit and dropping σ.
    cam = IdealCamera(8, 8, 0.1)
    T = Float64

    # Setting .dataset in place must work for the fixed-width / astigmatic types
    # and leave every other field untouched.
    e2g = GaussMLE.Emitter2DFitGaussMLE{T}(
        0.1, 0.2, 1000.0, 5.0,
        0.01, 0.012, 0.003, 20.0, 0.5, 0.4, 1, 1, 0, 1
    )
    e3g = GaussMLE.Emitter3DFitGaussMLE{T}(
        0.1, 0.2, 0.3, 1000.0, 5.0,
        0.01, 0.012, 0.02, 0.003, 0.001, 0.002, 20.0, 0.5, 0.4, 1, 1, 0, 1
    )
    e2g.dataset = 7
    @test e2g.dataset == 7
    @test e2g.σ_xy ≈ 0.003
    e3g.dataset = 7
    @test e3g.dataset == 7
    @test e3g.σ_yz ≈ 0.002

    mktempdir() do dir
        # Emitter2DFitGaussMLE (GaussianXYNB) round-trip.
        g2 = [
            GaussMLE.Emitter2DFitGaussMLE{T}(
                0.1i, 0.2i, 1000.0 + i, 5.0,
                0.01, 0.012, 0.003, 20.0, 0.5, 0.4, i, 1, 0, i
            ) for i in 1:4
        ]
        s2 = SMLMAnalysis.BasicSMLD(g2, cam, 10, 1, Dict{String, Any}())
        p2 = joinpath(dir, "g2.h5"); save_smld(p2, s2); l2 = load_smld(p2)
        @test eltype(l2.emitters) <: GaussMLE.Emitter2DFitGaussMLE
        for (a, b) in zip(s2.emitters, l2.emitters), f in fieldnames(GaussMLE.Emitter2DFitGaussMLE)
            @test getfield(a, f) ≈ getfield(b, f)
        end

        # Emitter3DFitGaussMLE (AstigmaticXYZNB) round-trip: z + full covariance.
        g3 = [
            GaussMLE.Emitter3DFitGaussMLE{T}(
                0.1i, 0.2i, 0.3i, 1000.0 + i, 5.0,
                0.01, 0.012, 0.02, 0.003, 0.001, 0.002, 20.0, 0.5, 0.4, i, 1, 0, i
            ) for i in 1:4
        ]
        s3 = SMLMAnalysis.BasicSMLD(g3, cam, 10, 1, Dict{String, Any}())
        p3 = joinpath(dir, "g3.h5"); save_smld(p3, s3); l3 = load_smld(p3)
        @test eltype(l3.emitters) <: GaussMLE.Emitter3DFitGaussMLE
        @test l3.emitters[2].z ≈ 0.6
        @test l3.emitters[2].σ_xz ≈ 0.001
        @test l3.emitters[2].σ_yz ≈ 0.002

        # Standard SMLMAnalysis.Emitter3DFit: off-diagonal covariances σ_xz/σ_yz survive
        # (they were never written before this fix).
        e3 = [
            SMLMAnalysis.Emitter3DFit{T}(
                0.1i, 0.2i, 0.3i, 1000.0 + i, 5.0,
                0.01, 0.012, 0.02, 20.0, 0.5;
                σ_xy = 0.003, σ_xz = 0.001, σ_yz = 0.002, frame = i, dataset = 1, id = i
            ) for i in 1:4
        ]
        s3s = SMLMAnalysis.BasicSMLD(e3, cam, 10, 1, Dict{String, Any}())
        p3s = joinpath(dir, "e3.h5"); save_smld(p3s, s3s); l3s = load_smld(p3s)
        @test eltype(l3s.emitters) <: SMLMAnalysis.Emitter3DFit
        @test l3s.emitters[2].σ_xz ≈ 0.001
        @test l3s.emitters[2].σ_yz ≈ 0.002

        # Abstract-eltype SMLD (as the pipeline produced before narrowing):
        # save_smld must reload it as concrete Emitter2DFitSigma, NOT degrade to
        # SMLMAnalysis.Emitter2DFit and drop the PSF-width σ.
        abs_v = SMLMAnalysis.AbstractEmitter[
            GaussMLE.Emitter2DFitSigma{T}(
                0.1i, 0.2i, 1000.0 + i, 5.0, 0.13,
                0.01, 0.012, 0.003, 20.0, 0.5, 0.002,
                0.4, i, 1, 0, i
            ) for i in 1:4
        ]
        @test eltype(abs_v) == SMLMAnalysis.AbstractEmitter
        s_abs = SMLMAnalysis.BasicSMLD(abs_v, cam, 10, 1, Dict{String, Any}())
        pabs = joinpath(dir, "abs.h5"); save_smld(pabs, s_abs); labs = load_smld(pabs)
        @test eltype(labs.emitters) <: GaussMLE.Emitter2DFitSigma   # NOT SMLMAnalysis.Emitter2DFit
        @test labs.emitters[2].σ ≈ 0.13                             # PSF-width σ preserved
        @test labs.emitters[2].σ_xy ≈ 0.003
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

@testset "edge geometry metadata round-trip" begin
    # Edge classification mirrors its cell mask into metadata; save_smld must keep it.
    CP = SMLMAnalysis.SMLMClustering.CellPolygon
    sq(x0, s) = NTuple{2, Float64}[(x0, x0), (x0 + s, x0), (x0 + s, x0 + s), (x0, x0 + s)]
    cells = [
        CP(sq(0.0, 4.0), [sq(0.5, 1.0), sq(2.0, 0.5)]),   # two holes
        CP(sq(5.0, 1.0)),                                # no holes
        CP(sq(7.0, 2.0), [sq(7.5, 0.2)]),
    ]
    outer = sq(0.0, 4.0)
    cam = IdealCamera(16, 16, 0.1)
    em = [SMLMAnalysis.Emitter2DFit{Float64}(0.1i, 0.1i, 1000.0, 5.0, 0.01, 0.01, 20.0, 0.5; frame = i) for i in 1:3]
    cellkey(cs) = [(c.outer, c.holes) for c in cs]   # CellPolygon has no ==; compare its fields
    mktempdir() do dir
        md = Dict{String, Any}(
            "edge_cells" => cells, "edge_outer_polygon" => outer,
            "empty_cells" => CP[], "empty_polygon" => NTuple{2, Float64}[]
        )
        p = joinpath(dir, "geom.h5")
        save_smld(p, SMLMAnalysis.BasicSMLD(em, cam, 3, 1, md))
        m2 = load_smld(p).metadata
        @test m2["edge_outer_polygon"] == outer
        @test m2["edge_outer_polygon"] isa Vector{NTuple{2, Float64}}
        @test m2["edge_cells"] isa Vector{CP}
        @test cellkey(m2["edge_cells"]) == cellkey(cells)
        @test m2["empty_cells"] isa Vector{CP} && isempty(m2["empty_cells"])
        @test m2["empty_polygon"] isa Vector{NTuple{2, Float64}} && isempty(m2["empty_polygon"])
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
