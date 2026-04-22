using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using Test
using Random

@testset "SMLMAnalysis.jl" begin
    @testset "Types" begin
        # Test AnalysisInfo constructor
        info = AnalysisInfo()
        @test info.elapsed_s == 0.0
        @test isempty(info.steps)
        @test isempty(info.step_infos)

        # Test AnalysisInfo with data
        steps = Dict{Symbol, Any}(:test => (a=1, b=2))
        info = AnalysisInfo(1.5, steps, StepInfo[])
        @test info.elapsed_s == 1.5
        @test info.steps[:test].a == 1

        # Test StepInfo with typed info
        cfg = FilterConfig()
        filter_info = FilterInfo(100, 80, 0.5)
        step_info = StepInfo(1, cfg, 0.5, Dict{Symbol,Any}(:n_before => 100); info=filter_info)
        @test step_info.info !== nothing
        @test step_info.info isa FilterInfo
        @test step_info.info.n_before == 100
        @test step_info.info.n_after == 80
        @test step_info.elapsed_s == 0.5

        # Test StepInfo without info
        step_info2 = StepInfo(2, cfg, 0.3, Dict{Symbol,Any}())
        @test step_info2.info === nothing

        # Test native info structs
        # Back-compat 7-arg constructor (defaults selected_source_indices to nothing)
        di = DetectFitInfo([], [], 2, 1000, 950, 5000, 1.5)
        @test di.n_datasets == 2
        @test di.n_rois == 1000
        @test di.n_fits == 950
        @test di.selected_source_indices === nothing

        # Full 8-arg constructor with provenance
        di_sel = DetectFitInfo([], [], 3, 500, 450, 1000, 0.5, [1, 3, 5])
        @test di_sel.selected_source_indices == [1, 3, 5]
        @test di_sel.n_datasets == 3

        dfi = DensityFilterInfo(1000, 800, 5, 0.3)
        @test dfi.n_before == 1000
        @test dfi.threshold == 5
    end

    @testset "analyze dispatch" begin
        # Verify analyze() dispatch methods exist for each step config type
        @test hasmethod(analyze, Tuple{Vector{<:AbstractArray{<:Real,3}}, DetectFitConfig})
        @test hasmethod(analyze, Tuple{AbstractArray{<:Real,3}, DetectFitConfig})
        @test hasmethod(analyze, Tuple{DetectFitConfig})
        @test hasmethod(analyze, Tuple{BasicSMLD, FilterConfig})
        @test hasmethod(analyze, Tuple{BasicSMLD, FrameConnectConfig})
        @test hasmethod(analyze, Tuple{BasicSMLD, DriftConfig})
        @test hasmethod(analyze, Tuple{BasicSMLD, DensityFilterConfig})
        @test hasmethod(analyze, Tuple{BasicSMLD, RenderConfig})

        # Verify old step function names are not exported
        @test !isdefined(Main, :detectfit)
        @test !isdefined(Main, :filter_step)
        @test !isdefined(Main, :frameconnect_step)
        @test !isdefined(Main, :driftcorrect_step)
        @test !isdefined(Main, :densityfilter_step)
        @test !isdefined(Main, :render_step)

        # DetectFitConfig camera field
        cfg = DetectFitConfig()
        @test cfg.camera === nothing
        cam = IdealCamera(64, 64, 0.1)
        cfg2 = DetectFitConfig(camera=cam, boxer=BoxerConfig(boxsize=7))
        @test cfg2.camera === cam
        @test cfg2.boxer.boxsize == 7
        @test cfg2.fitter.psf_model isa GaussianXYNBS

        # DetectFitConfig.datasets selection field
        @test cfg.datasets === nothing                          # default is no selection
        cfg_range = DetectFitConfig(datasets=1:19)
        @test cfg_range.datasets == 1:19
        @test cfg_range.datasets isa UnitRange{Int}
        cfg_sparse = DetectFitConfig(datasets=[1, 2, 3, 5, 7])
        @test cfg_sparse.datasets == [1, 2, 3, 5, 7]
        @test cfg_sparse.datasets isa Vector{Int}

        # _select_sources: pass-through when nothing, bounds-checked otherwise
        src = [(i=j,) for j in 1:5]
        @test SMLMAnalysis._select_sources(src, nothing) === src
        @test SMLMAnalysis._select_sources(src, [1, 3, 5]) == [src[1], src[3], src[5]]
        @test SMLMAnalysis._select_sources(src, 2:4) == src[2:4]
        @test_throws ErrorException SMLMAnalysis._select_sources(src, [1, 6])
        @test_throws ErrorException SMLMAnalysis._select_sources(src, [0, 1])
    end

    @testset "CalibrationConfig re-export" begin
        # CalibrationConfig is re-exported from SMLMFrameConnection
        @test CalibrationConfig === SMLMFrameConnection.CalibrationConfig
        @test CalibrationResult === SMLMFrameConnection.CalibrationResult

        # CalibrationConfig can be nested in FrameConnectConfig
        cal_cfg = CalibrationConfig(clamp_k_to_one=true)
        fc_cfg = FrameConnectConfig(max_frame_gap=5, calibration=cal_cfg)
        @test fc_cfg.calibration !== nothing
        @test fc_cfg.calibration.clamp_k_to_one == true

        # Default is nothing (no calibration)
        fc_default = FrameConnectConfig()
        @test fc_default.calibration === nothing
    end

    @testset "Info struct subtypes" begin
        # All info structs should be AbstractSMLMInfo subtypes
        @test StepInfo <: AbstractSMLMInfo
        @test DetectFitInfo <: AbstractSMLMInfo
        @test FilterInfo <: AbstractSMLMInfo
        @test DensityFilterInfo <: AbstractSMLMInfo
        @test CompositeRenderInfo <: AbstractSMLMInfo
        @test CrossAlignInfo <: AbstractSMLMInfo
        @test AnalysisInfo <: AbstractSMLMInfo
    end

    @testset "Multi-target step types" begin
        # Type hierarchy
        @test AbstractMultiTargetStep <: AbstractSMLMConfig
        @test CompositeRenderConfig <: AbstractMultiTargetStep
        @test CrossAlignConfig <: AbstractMultiTargetStep

        # CompositeRenderConfig defaults
        cr = CompositeRenderConfig()
        @test cr.strategy isa GaussianRender
        @test cr.zoom == 20.0
        @test cr.colors === nothing
        @test cr.clip_percentile == 0.99
        @test cr.normalize_each === nothing
        @test cr.scalebar == true
        @test cr.scalebar_position == :br

        # CompositeRenderConfig with custom fields
        cr2 = CompositeRenderConfig(strategy=HistogramRender(), zoom=10.0, colors=[:red, :blue])
        @test cr2.strategy isa HistogramRender
        @test cr2.zoom == 10.0
        @test cr2.colors == [:red, :blue]

        # CrossAlignConfig defaults
        ca = CrossAlignConfig()
        @test ca.method == :entropy
        @test ca.maxn == 100
        @test ca.histbinsize == 0.05

        # CrossAlignConfig custom
        ca2 = CrossAlignConfig(method=:fft, maxn=50)
        @test ca2.method == :fft
        @test ca2.maxn == 50

        # step_name dispatch
        @test SMLMAnalysis.step_name(cr) == "compositerender"
        @test SMLMAnalysis.step_name(ca) == "crossalign"

        # analyze dispatch methods exist for multi-target steps
        @test hasmethod(analyze, Tuple{Vector{<:BasicSMLD}, CompositeRenderConfig})
        @test hasmethod(analyze, Tuple{Vector{<:BasicSMLD}, CrossAlignConfig})

        # MultiTargetConfig with steps vector
        mt = MultiTargetConfig(
            labels=[:A, :B],
            steps=[
                CompositeRenderConfig(zoom=20.0),
                CrossAlignConfig(),
                CompositeRenderConfig(zoom=10.0, strategy=HistogramRender()),
            ],
            outdir="/tmp/test_mt",
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

    @testset "Bleaching fit degeneracy guard" begin
        Random.seed!(42)

        # Flat data: constant ~30 locs/frame with Poisson-like noise, 19k frames.
        # Reproduces the @hstirf case that produced a=-42873, k≈0 with unbounded NelderMead.
        # Expected: reject as degenerate, return nothing.
        flat = [max(0, round(Int, 30 + randn() * sqrt(30))) for _ in 1:19_000]
        @test SMLMAnalysis._estimate_bleaching_rate(flat) === nothing

        # Pure exponential decay: should recover parameters accurately (non-regression).
        t_exp = 1:5000
        exp_data = [max(0, round(Int, 10 + 50 * exp(-0.0005 * i) + randn() * 2)) for i in t_exp]
        res_exp = SMLMAnalysis._estimate_bleaching_rate(exp_data)
        @test res_exp !== nothing
        @test isapprox(res_exp.k_bleach, 5e-4, rtol=0.1)
        @test isapprox(res_exp.offset, 10.0, atol=2.0)
        @test isapprox(res_exp.N_0, 50.0, atol=3.0)
        @test res_exp.r_squared > 0.9

        # Bleach-then-flat (realistic DNA-PAINT): should still fit with physical params.
        t_bf = 1:10_000
        bf_data = [max(0, round(Int, 20 + 30 * exp(-0.001 * i) + randn() * 1.5)) for i in t_bf]
        res_bf = SMLMAnalysis._estimate_bleaching_rate(bf_data)
        @test res_bf !== nothing
        @test res_bf.offset >= 0     # physical bound held
        @test res_bf.N_0 >= 0        # physical bound held
        @test res_bf.k_bleach > 0
    end
end
