using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Random
using TOML
using Statistics

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

@testset "config.toml keeps each key in its own table" begin
    # A TOML table owns every key after its header. DetectFitConfig's plain fields
    # (h5_format, qe) come after its nested boxer and fitter configs, and GaussMLEConfig's
    # iterations after its psf_model; each must parse back where it belongs.
    io = IOBuffer()
    SMLMAnalysis._write_config_fields!(io, DetectFitConfig())
    parsed = TOML.parse(String(take!(io)))
    @test parsed["h5_format"] == "auto"
    @test parsed["qe"] == 1.0
    @test !haskey(parsed["boxer"], "qe")
    @test !haskey(parsed["fitter"], "qe")
    @test parsed["fitter"]["iterations"] == 20
    @test !haskey(parsed["fitter"]["psf_model"], "iterations")
    @test parsed["fitter"]["psf_model"]["type"] == "GaussianXYNBS"
end
