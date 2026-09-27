using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Random
using TOML
using Statistics

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
