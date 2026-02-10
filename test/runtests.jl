using SMLMAnalysis
using Test

@testset "SMLMAnalysis.jl" begin
    @testset "Types" begin
        # Test AnalysisInfo constructor
        info = AnalysisInfo()
        @test info.elapsed_s == 0.0
        @test isempty(info.steps)

        # Test AnalysisInfo with data
        steps = Dict{Symbol, Any}(:test => (a=1, b=2))
        info = AnalysisInfo(1.5, steps, StepRecord[])
        @test info.elapsed_s == 1.5
        @test info.steps[:test].a == 1

        # Test StepRecord with info
        cfg = FilterConfig()
        record = StepRecord(1, cfg, 0.5, Dict{Symbol,Any}(); info=(test=true,))
        @test record.info !== nothing
        @test record.info.test == true

        # Test StepRecord without info
        record2 = StepRecord(2, cfg, 0.3, Dict{Symbol,Any}())
        @test record2.info === nothing
    end

    @testset "analyze dispatch" begin
        # Verify analyze() dispatch methods exist for each step config type
        @test hasmethod(analyze, Tuple{Vector{<:AbstractArray{<:Real,3}}, DetectFitConfig})
        @test hasmethod(analyze, Tuple{AbstractArray{<:Real,3}, DetectFitConfig})
        @test hasmethod(analyze, Tuple{DetectFitConfig})
        @test hasmethod(analyze, Tuple{BasicSMLD, FilterConfig})
        @test hasmethod(analyze, Tuple{BasicSMLD, FrameConnectConfig})
        @test hasmethod(analyze, Tuple{BasicSMLD, DriftCorrectConfig})
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
        cfg2 = DetectFitConfig(camera=cam, boxsize=7)
        @test cfg2.camera === cam
        @test cfg2.boxsize == 7
    end
end
