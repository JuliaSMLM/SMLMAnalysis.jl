using SMLMAnalysis
using Test

@testset "SMLMAnalysis.jl" begin
    @testset "Types" begin
        # Test AnalysisInfo constructor
        info = AnalysisInfo()
        @test info.elapsed_ns == UInt64(0)
        @test isempty(info.steps)

        # Test AnalysisInfo with data
        steps = Dict{Symbol, Any}(:test => (a=1, b=2))
        info = AnalysisInfo(UInt64(1000), steps)
        @test info.elapsed_ns == UInt64(1000)
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
end
