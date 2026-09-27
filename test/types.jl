using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Random
using TOML
using Statistics

@testset "Verbosity/Checkpoint validation" begin
    cam = IdealCamera(64, 64, 0.1)

    # Valid levels round-trip through AnalysisConfig / MultiTargetConfig construction.
    @test AnalysisConfig(camera = cam, verbose = Verbosity.SILENT, checkpoint = Checkpoint.NONE).verbose == Verbosity.SILENT
    @test AnalysisConfig(camera = cam, verbose = Verbosity.DEBUG, checkpoint = Checkpoint.ALL).checkpoint == Checkpoint.ALL
    @test MultiTargetConfig(labels = [:A], outdir = "x", verbose = Verbosity.DEBUG).verbose == Verbosity.DEBUG

    # Out-of-range verbose/checkpoint must raise ArgumentError at construction.
    @test_throws ArgumentError AnalysisConfig(camera = cam, verbose = -1)
    @test_throws ArgumentError AnalysisConfig(camera = cam, verbose = Verbosity.DEBUG + 1)
    @test_throws ArgumentError AnalysisConfig(camera = cam, checkpoint = -1)
    @test_throws ArgumentError AnalysisConfig(camera = cam, checkpoint = Checkpoint.ALL + 1)
    @test_throws ArgumentError MultiTargetConfig(labels = [:A], outdir = "x", verbose = -1)
    @test_throws ArgumentError MultiTargetConfig(labels = [:A], outdir = "x", verbose = Verbosity.DEBUG + 1)

    # Same validation at the _run_pipeline entry point, for direct calls that
    # bypass AnalysisConfig entirely.
    steps = SMLMAnalysis.AbstractSMLMConfig[]
    @test_throws ArgumentError SMLMAnalysis._run_pipeline(nothing, steps, cam, nothing, -1)
    @test_throws ArgumentError SMLMAnalysis._run_pipeline(nothing, steps, cam, nothing, Verbosity.STANDARD, -1)
    @test_throws ArgumentError SMLMAnalysis._run_pipeline(nothing, steps, cam, nothing, Verbosity.STANDARD, Checkpoint.ALL + 1)
end

@testset "CalibrationConfig re-export" begin
    # CalibrationConfig is re-exported from SMLMFrameConnection
    @test CalibrationConfig === SMLMFrameConnection.CalibrationConfig
    @test CalibrationResult === SMLMFrameConnection.CalibrationResult

    # CalibrationConfig can be nested in FrameConnectConfig
    cal_cfg = CalibrationConfig(clamp_k_to_one = true)
    fc_cfg = FrameConnectConfig(max_frame_gap = 5, calibration = cal_cfg)
    @test fc_cfg.calibration !== nothing
    @test fc_cfg.calibration.clamp_k_to_one == true

    # Default is nothing (no calibration)
    fc_default = FrameConnectConfig()
    @test fc_default.calibration === nothing
end

@testset "Info struct subtypes" begin
    # All info structs should be SMLMAnalysis.AbstractSMLMInfo subtypes
    @test SMLMAnalysis.StepInfo <: SMLMAnalysis.AbstractSMLMInfo
    @test SMLMAnalysis.DetectFitInfo <: SMLMAnalysis.AbstractSMLMInfo
    @test SMLMAnalysis.FilterInfo <: SMLMAnalysis.AbstractSMLMInfo
    @test SMLMAnalysis.DensityFilterInfo <: SMLMAnalysis.AbstractSMLMInfo
    @test SMLMAnalysis.CompositeRenderInfo <: SMLMAnalysis.AbstractSMLMInfo
    @test SMLMAnalysis.CrossAlignInfo <: SMLMAnalysis.AbstractSMLMInfo
    @test SMLMAnalysis.AnalysisInfo <: SMLMAnalysis.AbstractSMLMInfo
end

@testset "stepinfo/stepinfos accessors" begin
    # Build StepInfos directly; FilterConfig → name "filter" (repeated),
    # DensityFilterConfig → name "densityfilter" (unique).
    si_a = SMLMAnalysis.StepInfo(1, FilterConfig(), 0.1, Dict{Symbol, Any}(); info = SMLMAnalysis.FilterInfo(100, 90, 0.1))
    si_b = SMLMAnalysis.StepInfo(2, DensityFilterConfig(), 0.2, Dict{Symbol, Any}(); info = SMLMAnalysis.DensityFilterInfo(90, 80, 5, 0.2))
    si_c = SMLMAnalysis.StepInfo(3, FilterConfig(), 0.3, Dict{Symbol, Any}(); info = SMLMAnalysis.FilterInfo(200, 150, 0.3))
    steps = SMLMAnalysis.StepInfo[si_a, si_b, si_c]

    for info in (
            SMLMAnalysis.AnalysisInfo(1.0, steps),
            SMLMAnalysis.MultiTargetInfo(1.0, Dict{Symbol, SMLMAnalysis.AnalysisInfo}(), steps),
        )
        # Symbol and String name lookup return the same SMLMAnalysis.StepInfo
        @test stepinfo(info, :densityfilter) === si_b
        @test stepinfo(info, "densityfilter") === si_b
        @test stepinfo(info, :densityfilter).info isa SMLMAnalysis.DensityFilterInfo

        # Repeated name: stepinfo returns the FIRST, stepinfos returns ALL in order
        @test stepinfo(info, :filter) === si_a
        @test stepinfo(info, "filter") === si_a
        @test stepinfos(info, :filter) == [si_a, si_c]
        @test [si.number for si in stepinfos(info, "filter")] == [1, 3]

        # Missing name: stepinfo throws ArgumentError listing available names,
        # stepinfos returns empty
        @test_throws ArgumentError stepinfo(info, :nope)
        @test isempty(stepinfos(info, :nope))
    end
end
