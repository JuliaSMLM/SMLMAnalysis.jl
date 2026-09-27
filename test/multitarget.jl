using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Random
using TOML
using Statistics

# Test-local fixtures for the "multi-target steps receive labels and colors"
# testset below. `struct` and `analyze` method definitions must be at module
# top level (not inside a `@testset` body), so they live here.

# A user `AbstractMultiTargetStep` that only records the keywords the
# orchestrator called it with -- stands in for a real cross-channel step.
const _recorded_labels_colors = Ref{Any}(nothing)

struct _RecordingMultiTargetStep <: SMLMAnalysis.AbstractMultiTargetStep end

function SMLMAnalysis.analyze(
        smlds::Vector{<:SMLMAnalysis.BasicSMLD}, cfg::_RecordingMultiTargetStep;
        outdir = nothing, step_number::Int = 0, verbose::Int = 0,
        labels::Vector{Symbol} = Symbol[], colors::Vector{Symbol} = Symbol[]
    )
    _recorded_labels_colors[] = (labels, colors)
    return (smlds, SMLMAnalysis.StepInfo(step_number, cfg, 0.0, Dict{Symbol, Any}()))
end

# A cheap per-channel step standing in for DetectFitConfig: hands back a
# pre-built BasicSMLD without running real detection/fitting, so phase 1 of
# `analyze(channels, MultiTargetConfig)` stays fast in this test.
struct _FakeChannelStep <: SMLMAnalysis.AbstractSMLMConfig
    smld::SMLMAnalysis.BasicSMLD
end

function SMLMAnalysis.analyze(
        ::Any, cfg::_FakeChannelStep;
        outdir = nothing, step_number::Int = 0, verbose::Int = 0,
        checkpoint::Int = Checkpoint.EXPENSIVE
    )
    return (cfg.smld, SMLMAnalysis.StepInfo(step_number, cfg, 0.0, Dict{Symbol, Any}()))
end

@testset "Multi-target step types" begin
    # Type hierarchy
    @test SMLMAnalysis.AbstractMultiTargetStep <: SMLMAnalysis.AbstractSMLMConfig
    @test CompositeRenderConfig <: SMLMAnalysis.AbstractMultiTargetStep
    @test CrossAlignConfig <: SMLMAnalysis.AbstractMultiTargetStep

    # CompositeRenderConfig defaults
    cr = CompositeRenderConfig()
    @test cr.strategy isa GaussianRender
    @test cr.zoom == 20.0
    @test cr.colors === nothing
    @test cr.clip_percentile === :auto
    @test cr.normalize_each === nothing
    @test cr.scalebar == true
    @test cr.scalebar_position == :br

    # CompositeRenderConfig with custom fields
    cr2 = CompositeRenderConfig(strategy = HistogramRender(), zoom = 10.0, colors = [:red, :blue])
    @test cr2.strategy isa HistogramRender
    @test cr2.zoom == 10.0
    @test cr2.colors == [:red, :blue]

    # CrossAlignConfig defaults
    ca = CrossAlignConfig()
    @test ca.align isa AlignConfig
    @test ca.align.method == :entropy
    @test ca.align.maxn == 100
    @test ca.align.histbinsize == 0.05

    # CrossAlignConfig custom: upstream AlignConfig passed through as-is
    ca2 = CrossAlignConfig(align = AlignConfig(method = :fft, maxn = 50, verbose = 1))
    @test ca2.align.method == :fft
    @test ca2.align.maxn == 50
    @test ca2.align.verbose == 1

    # step_name dispatch
    @test SMLMAnalysis.step_name(cr) == "compositerender"
    @test SMLMAnalysis.step_name(ca) == "crossalign"

    # analyze dispatch methods exist for multi-target steps
    @test hasmethod(analyze, Tuple{Vector{<:SMLMAnalysis.BasicSMLD}, CompositeRenderConfig})
    @test hasmethod(analyze, Tuple{Vector{<:SMLMAnalysis.BasicSMLD}, CrossAlignConfig})

    # MultiTargetConfig with steps vector
    mt = MultiTargetConfig(
        labels = [:A, :B],
        steps = [
            CompositeRenderConfig(zoom = 20.0),
            CrossAlignConfig(),
            CompositeRenderConfig(zoom = 10.0, strategy = HistogramRender()),
        ],
        outdir = "/tmp/test_mt",
    )
    @test length(mt.steps) == 3
    @test mt.steps[1] isa CompositeRenderConfig
    @test mt.steps[2] isa CrossAlignConfig
    @test mt.steps[3] isa CompositeRenderConfig
    @test mt.colors == [:cyan, :magenta]

    # AlignConfig/AlignInfo re-exports
    @test AlignConfig === SMLMDriftCorrection.AlignConfig
    @test AlignInfo === SMLMDriftCorrection.AlignInfo
end

@testset "multi-target config.toml is valid with two nested-config steps" begin
    # CompositeRenderConfig (strategy) and CrossAlignConfig (align) each carry a
    # nested config field. Written into the same multi_target_config.toml under
    # separate [[steps]] entries, a bare `[strategy]` / `[align]` table floats to
    # the document root and collides across steps -- TOML.parsefile then fails
    # with "key already defined". table_prefix="steps." scopes each nested table
    # to its own [[steps]] array element instead.
    mktempdir() do dir
        mt = MultiTargetConfig(
            labels = [:A, :B],
            steps = [
                CompositeRenderConfig(),
                CrossAlignConfig(),
                CompositeRenderConfig(),
            ],
            outdir = dir,
        )
        SMLMAnalysis._save_multitarget_config!(mt)
        parsed = TOML.parsefile(joinpath(dir, "multi_target_config.toml"))
        @test length(parsed["steps"]) == 3
        @test haskey(parsed["steps"][1], "strategy")
        @test haskey(parsed["steps"][2], "align")
        @test haskey(parsed["steps"][3], "strategy")
    end
end

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

@testset "multi-target steps receive labels and colors" begin
    # Every multi-target step's analyze() is called with both `labels` and
    # `colors`, whether or not it reads either -- this must fail before the
    # fix (the old `_multitarget_extra_kwargs` dispatch only forwarded the one
    # keyword each built-in step declared, so a user step reading neither, or
    # both, never saw them).
    cam = IdealCamera(8, 8, 0.1)
    mk(dx) = SMLMAnalysis.BasicSMLD(
        [
            SMLMAnalysis.Emitter2DFit{Float64}(
                1.0 + dx, 1.0, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0; frame = 1
            ),
        ],
        cam, 1, 1, Dict{String, Any}()
    )
    channel_cfg(smld) = AnalysisConfig(camera = cam, steps = [_FakeChannelStep(smld)])

    mt = MultiTargetConfig(
        labels = [:A, :B], colors = [:red, :green],
        steps = [_RecordingMultiTargetStep()],
        outdir = mktempdir(),
    )
    _recorded_labels_colors[] = nothing
    analyze(
        [
            (zeros(Float32, 1, 1, 1), channel_cfg(mk(0.0))),
            (zeros(Float32, 1, 1, 1), channel_cfg(mk(0.1))),
        ], mt
    )
    @test _recorded_labels_colors[] == (mt.labels, mt.colors)
end

@testset "cross-align carries edge geometry with the emitters" begin
    # align_smld only moves emitters; edge-classification metadata
    # (edge_outer_polygon / edge_cells) must move with them or a saved
    # post-alignment SMLD pairs aligned emitters with a stale mask (Codex #51).
    CP = SMLMAnalysis.SMLMClustering.CellPolygon
    rng = MersenneTwister(11)
    cam = IdealCamera(64, 64, 0.1)
    N = 1500
    xs = 1.0 .+ 4.4 .* rand(rng, N); ys = 1.0 .+ 4.4 .* rand(rng, N)
    sq(x0, s) = NTuple{2, Float64}[(x0, x0), (x0 + s, x0), (x0 + s, x0 + s), (x0, x0 + s)]
    shiftpts(pts, dx, dy) = NTuple{2, Float64}[(p[1] + dx, p[2] + dy) for p in pts]

    outer_a = sq(1.0, 4.0)
    cells_a = [CP(sq(1.0, 4.0), [sq(1.5, 1.0)])]
    dx, dy = 0.094, 0.0
    outer_b = shiftpts(outer_a, dx, dy)
    cells_b = [
        CP(
            shiftpts(cells_a[1].outer, dx, dy),
            [shiftpts(h, dx, dy) for h in cells_a[1].holes]
        ),
    ]

    mkgeom(ddx, ddy, md) = SMLMAnalysis.BasicSMLD(
        [
            SMLMAnalysis.Emitter2DFit{Float64}(
                xs[i] + ddx, ys[i] + ddy, 1000.0,
                10.0, 0.01, 0.01, 50.0, 2.0; frame = 1 + (i % 10)
            ) for i in 1:N
        ],
        cam, 10, 1, md
    )
    a = mkgeom(0.0, 0.0, Dict{String, Any}("edge_outer_polygon" => outer_a, "edge_cells" => cells_a))
    b = mkgeom(dx, dy, Dict{String, Any}("edge_outer_polygon" => outer_b, "edge_cells" => cells_b))

    # Default CrossAlignConfig() uses AlignConfig's default transform=:shift, so
    # every vertex and every emitter is offset by the exact same [dx, dy] —
    # emitter shift and polygon-vertex shift must agree exactly (not just close).
    (aligned, si) = analyze([a, b], CrossAlignConfig(); verbose = 0)
    @test si.info isa SMLMAnalysis.CrossAlignInfo

    mean_x(s) = sum(e.x for e in s.emitters) / length(s.emitters)
    mean_y(s) = sum(e.y for e in s.emitters) / length(s.emitters)
    emitter_dx = mean_x(b) - mean_x(aligned[2])
    emitter_dy = mean_y(b) - mean_y(aligned[2])

    # Every vertex (outer ring, cell outer ring, and the cell's hole) must have
    # moved by the same [dx, dy] the emitters did, not merely a nearby amount.
    vertex_shift_matches(before, after) = all(
        isapprox(pb[1] - pa[1], emitter_dx; atol = 1.0e-9) &&
            isapprox(pb[2] - pa[2], emitter_dy; atol = 1.0e-9)
            for (pb, pa) in zip(before, after)
    )

    @test vertex_shift_matches(outer_b, aligned[2].metadata["edge_outer_polygon"])
    aligned_cell = aligned[2].metadata["edge_cells"][1]
    @test vertex_shift_matches(cells_b[1].outer, aligned_cell.outer)
    @test vertex_shift_matches(cells_b[1].holes[1], aligned_cell.holes[1])

    # The reference channel (aligned[1] === a) is never touched.
    @test aligned[1] === a
    @test aligned[1].metadata["edge_outer_polygon"] == outer_a
    @test aligned[1].metadata["edge_cells"][1].outer == cells_a[1].outer
end

@testset "cross-align drops edge geometry under :affine alignment" begin
    # align_smld's :affine transform composes two sequential affine passes,
    # but info.diagnostic only records each pass's own coefficients summed
    # together, which is NOT the composed map (Codex measured a 50nm
    # mismatch at x=50μm with a=0.1, a2=0.01). Recovering the exact composed
    # map from the aligned emitter positions is unsound in general
    # (degenerate/collinear layouts give a spurious exact fit) and Float32
    # data can legitimately fail a tight residual tolerance, so the ruling is
    # to not carry geometry through :affine at all: both metadata keys are
    # dropped, with one @warn per channel.
    CP = SMLMAnalysis.SMLMClustering.CellPolygon
    cam = IdealCamera(16, 16, 0.1)
    em = [SMLMAnalysis.Emitter2DFit{Float64}(0.1i, 0.1i, 1000.0, 5.0, 0.01, 0.01, 20.0, 0.5; frame = i) for i in 1:3]
    outer = NTuple{2, Float64}[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)]
    hole = NTuple{2, Float64}[(1.5, 1.5), (2.0, 1.5), (2.0, 2.0)]
    md = Dict{String, Any}("edge_outer_polygon" => outer, "edge_cells" => [CP(outer, [hole])])
    ref = SMLMAnalysis.BasicSMLD(em, cam, 3, 1, Dict{String, Any}())
    chan = SMLMAnalysis.BasicSMLD(em, cam, 3, 1, deepcopy(md))
    aligned = [ref, chan]

    info = SMLMDriftCorrection.AlignInfo(
        [zeros(2), zeros(2)],
        SMLMDriftCorrection.AbstractAlignTransform[
            SMLMDriftCorrection.AffineTransform2D(0.0, 1.0, 0.0, 0.0),
            SMLMDriftCorrection.AffineTransform2D(0.0, 1.0, 0.0, 0.0),
        ],
        0.01, :entropy, :affine, :cpu, nothing
    )

    @test_logs (:warn, r"not carried through :affine") SMLMAnalysis._align_edge_geometry!(aligned, info)
    @test !haskey(aligned[2].metadata, "edge_outer_polygon")
    @test !haskey(aligned[2].metadata, "edge_cells")
    @test aligned[1].metadata == Dict{String, Any}()   # reference untouched
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
