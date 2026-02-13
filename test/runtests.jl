using SMLMAnalysis
using SMLMFrameConnection
using Test

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
        di = DetectFitInfo([], [], 2, 1000, 950, 5000, 1.5)
        @test di.n_datasets == 2
        @test di.n_rois == 1000
        @test di.n_fits == 950

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
        @test AnalysisInfo <: AbstractSMLMInfo
    end
end
