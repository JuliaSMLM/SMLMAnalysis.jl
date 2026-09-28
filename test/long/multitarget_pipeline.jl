using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Random
using TOML
using Statistics

@testset "multi-target saves and returns aligned channels" begin
    # Phase-2 dispatch + result assembly on hand-made SMLDs (no detectfit): channel
    # B is channel A shifted by a known offset; after CrossAlign the saved file,
    # result[:B].smld and result.smlds[2] must all be the aligned data.
    rng = MersenneTwister(3)
    cam = IdealCamera(64, 64, 0.1)
    N = 1500
    xs = 1.0 .+ 4.4 .* rand(rng, N); ys = 1.0 .+ 4.4 .* rand(rng, N)
    mk(dx, dy) = SMLMAnalysis.BasicSMLD(
        [
            SMLMAnalysis.Emitter2DFit{Float64}(
                xs[i] + dx, ys[i] + dy, 1000.0, 10.0,
                0.01, 0.01, 50.0, 2.0; frame = 1 + (i % 10)
            ) for i in 1:N
        ],
        cam, 10, 1, Dict{String, Any}()
    )
    a, b = mk(0.0, 0.0), mk(0.08, -0.05)
    labels = [:A, :B]
    (state, si) = analyze(
        [a, b], CrossAlignConfig();
        outdir = nothing, step_number = 1, verbose = 0
    )
    meanxy(s) = (sum(e.x for e in s.emitters) / N, sum(e.y for e in s.emitters) / N)
    off(s1, s2) = hypot((meanxy(s2) .- meanxy(s1))...)
    @test off(state[1], state[2]) < 0.01          # known ~94 nm offset removed to < 10 nm

    channels = Dict{Symbol, SMLMAnalysis.AnalysisResult}(
        :A => SMLMAnalysis.AnalysisResult(a, a, nothing),
        :B => SMLMAnalysis.AnalysisResult(b, b, nothing)
    )
    mktempdir() do dir
        SMLMAnalysis._finalize_channels!(channels, state, labels, dir; verbose = 0)
        mtr = SMLMAnalysis.MultiTargetResult(labels, state, channels, [si], dir)
        @test mtr[:B].smld === mtr.smlds[2] === state[2]
        @test mtr[:B].smld_connected === b                # pre-alignment data kept
        saved = load_smld(joinpath(dir, "smld_B.h5"))
        @test [e.x for e in saved.emitters] == [e.x for e in state[2].emitters]
        @test off(saved, b) > 0.05                        # file holds aligned, not raw, B
    end

    # A state that is not one SMLD per label is a contract violation: throw
    # rather than silently leave the channel results alone (a caller that
    # continued on to _write_composite_readme! would otherwise BoundsError).
    ch2 = Dict{Symbol, SMLMAnalysis.AnalysisResult}(
        :A => SMLMAnalysis.AnalysisResult(a, nothing, nothing),
        :B => SMLMAnalysis.AnalysisResult(b, nothing, nothing)
    )
    mktempdir() do dir
        @test_throws ArgumentError SMLMAnalysis._finalize_channels!(ch2, state[1:1], labels, dir; verbose = 0)
        @test ch2[:B].smld === b
        @test isempty(readdir(dir))

        # Right length, but an untyped Vector{Any} container: rejected even
        # though its actual elements are BasicSMLDs -- _write_composite_readme!
        # requires Vector{<:SMLMAnalysis.BasicSMLD}, so this would otherwise MethodError
        # there instead of failing loudly here.
        @test_throws ArgumentError SMLMAnalysis._finalize_channels!(ch2, Any[state[1], state[2]], labels, dir; verbose = 0)

        # A view is an AbstractVector{<:SMLMAnalysis.BasicSMLD} but not a Vector -- same
        # rejection, for the same reason (_write_composite_readme! requires
        # exactly Vector{<:SMLMAnalysis.BasicSMLD}).
        @test_throws ArgumentError SMLMAnalysis._finalize_channels!(ch2, view(state, 1:2), labels, dir; verbose = 0)
    end
end

@testset "crosscorr g(r)" begin
    # Exercises the internal _compute_crosscorr for the three fixed defects:
    # (a) bin-width consistency when r_max is NOT an integer multiple of dr,
    # (b) exactly-coincident cross-channel pairs (dist==0) are counted, not
    #     dropped, and (c) CSR normalization. CPU-only; no GPU needed.
    T = Float64

    # Minimal emitter builder — only x/y matter to _compute_crosscorr; every
    # other field is a constant dummy. Layout mirrors the Emitter2DFitSigma
    # constructor used in the HDF5 round-trip testset above.
    mkem(x, y, i) = GaussMLE.Emitter2DFitSigma{T}(
        x, y, 1000.0, 5.0, 0.13,
        0.01, 0.012, 0.003, 20.0, 0.5, 0.002,
        0.4, 1, 1, 0, i
    )
    mksmld(xs, ys, cam) = SMLMAnalysis.BasicSMLD(
        [mkem(xs[i], ys[i], i) for i in eachindex(xs)],
        cam, 1, 1, Dict{String, Any}()
    )

    cam = IdealCamera(64, 64, 0.1)   # 6.4 × 6.4 μm FOV
    xlo, xhi = first(cam.pixel_edges_x), last(cam.pixel_edges_x)
    ylo, yhi = first(cam.pixel_edges_y), last(cam.pixel_edges_y)
    nx() = xlo + rand() * (xhi - xlo)
    ny() = ylo + rand() * (yhi - ylo)

    # (a) Non-divisor consistency: r_max=1.0 is not a multiple of dr=0.03.
    Random.seed!(7)
    cfg_nd = CrossCorrConfig(r_max = 1.0, dr = 0.03)
    sa = mksmld([nx() for _ in 1:50], [ny() for _ in 1:50], cam)
    sb = mksmld([nx() for _ in 1:50], [ny() for _ in 1:50], cam)
    (r_nd, g_nd, area_nd) = SMLMAnalysis._compute_crosscorr(sa, sb, cfg_nd)
    # r_centers are spaced by exactly dr — one shared width for counting,
    # annulus areas, and the r-axis.
    @test all(isapprox.(diff(r_nd), cfg_nd.dr; atol = 1.0e-9))
    # Last bin's outer edge (= last center + dr/2) reaches at least r_max.
    @test (last(r_nd) + cfg_nd.dr / 2) >= cfg_nd.r_max - 1.0e-9
    @test length(r_nd) == ceil(Int, cfg_nd.r_max / cfg_nd.dr)

    # (b) Zero-distance kept: A and B share exactly-coincident points (plus
    # random background). The coincident pairs must elevate the first bin.
    Random.seed!(11)
    ncoinc = 30
    # Coincident locations on a diagonal, spacing ≈ 0.27 μm ≫ dr, so only the
    # exact same-location A/B pairs (dist==0) fall into bin 1.
    cxs = [0.4 + 0.19k for k in 0:(ncoinc - 1)]
    cys = [0.4 + 0.19k for k in 0:(ncoinc - 1)]
    bg_ax = [nx() for _ in 1:100]; bg_ay = [ny() for _ in 1:100]
    bg_bx = [nx() for _ in 1:100]; bg_by = [ny() for _ in 1:100]
    sa_z = mksmld(vcat(cxs, bg_ax), vcat(cys, bg_ay), cam)
    sb_z = mksmld(vcat(cxs, bg_bx), vcat(cys, bg_by), cam)
    (r_z, g_z, _) = SMLMAnalysis._compute_crosscorr(sa_z, sb_z, CrossCorrConfig())
    @test g_z[1] > 1   # coincident cross-channel pairs are counted, not dropped
    @test g_z[1] > 2   # and clearly elevated above the CSR background

    # (c) CSR sanity: two independent uniform-random channels give g(r) ≈ 1
    # across mid-range bins (away from r→0 noise and the far edge).
    Random.seed!(1234)
    npts = 3000
    ax = [nx() for _ in 1:npts]; ay = [ny() for _ in 1:npts]
    bx = [nx() for _ in 1:npts]; by = [ny() for _ in 1:npts]
    sa_c = mksmld(ax, ay, cam)
    sb_c = mksmld(bx, by, cam)
    (r_c, g_c, _) = SMLMAnalysis._compute_crosscorr(sa_c, sb_c, CrossCorrConfig())
    mid = 15:35        # r ≈ 0.145–0.345 μm
    @test isapprox(sum(g_c[mid]) / length(mid), 1.0; atol = 0.1)
    @test all(0.7 .< g_c[mid] .< 1.3)
end
