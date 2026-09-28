using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Random
using TOML
using Statistics

@testset "Types" begin
    # Test SMLMAnalysis.AnalysisInfo constructor
    info = SMLMAnalysis.AnalysisInfo()
    @test info.elapsed_s == 0.0
    @test isempty(info.step_infos)

    # Test SMLMAnalysis.AnalysisInfo with data
    cfg0 = FilterConfig()
    si0 = SMLMAnalysis.StepInfo(1, cfg0, 0.2, Dict{Symbol, Any}(); info = SMLMAnalysis.FilterInfo(10, 8, 0.2))
    info = SMLMAnalysis.AnalysisInfo(1.5, SMLMAnalysis.StepInfo[si0])
    @test info.elapsed_s == 1.5
    @test length(info.step_infos) == 1
    @test info.step_infos[1].info isa SMLMAnalysis.FilterInfo

    # Test SMLMAnalysis.StepInfo with typed info
    cfg = FilterConfig()
    filter_info = SMLMAnalysis.FilterInfo(100, 80, 0.5)
    step_info = SMLMAnalysis.StepInfo(1, cfg, 0.5, Dict{Symbol, Any}(:n_before => 100); info = filter_info)
    @test step_info.info !== nothing
    @test step_info.info isa SMLMAnalysis.FilterInfo
    @test step_info.info.n_before == 100
    @test step_info.info.n_after == 80
    @test step_info.elapsed_s == 0.5

    # Test SMLMAnalysis.StepInfo without info
    step_info2 = SMLMAnalysis.StepInfo(2, cfg, 0.3, Dict{Symbol, Any}())
    @test step_info2.info === nothing

    # Test native info structs
    # Back-compat 7-arg constructor (defaults selected_source_indices to nothing)
    di = SMLMAnalysis.DetectFitInfo([], [], 2, 1000, 950, 5000, 1.5)
    @test di.n_datasets == 2
    @test di.n_rois == 1000
    @test di.n_fits == 950
    @test di.selected_source_indices === nothing

    @test di.n_nonfinite == 0

    # Full 8-arg constructor with provenance
    di_sel = SMLMAnalysis.DetectFitInfo([], [], 3, 500, 450, 1000, 0.5, [1, 3, 5])
    @test di_sel.selected_source_indices == [1, 3, 5]
    @test di_sel.n_datasets == 3
    @test di_sel.n_nonfinite == 0

    # _drop_nonfinite_emitters: drops emitters with a NaN/Inf in any AbstractFloat field.
    e_ok = SMLMAnalysis.Emitter2DFit(1.0, 1.0, 100.0, 5.0, 0.01, 0.01, 0.0, 1.0, 1.0, 1, 1, 0, 1)
    e_nan_y = SMLMAnalysis.Emitter2DFit(1.0, NaN, 100.0, 5.0, 0.01, 0.01, 0.0, 1.0, 1.0, 1, 1, 0, 2)
    e_inf_photons = SMLMAnalysis.Emitter2DFit(1.0, 1.0, Inf, 5.0, 0.01, 0.01, 0.0, 1.0, 1.0, 1, 1, 0, 3)
    kept, n_dropped = SMLMAnalysis._drop_nonfinite_emitters([e_ok, e_nan_y, e_inf_photons])
    @test kept == [e_ok]
    @test n_dropped == 2

    # A NaN in a diagnostic field (pvalue) is NOT a required-finite field --
    # the emitter must be kept, not dropped.
    e_nan_pvalue = GaussMLE.Emitter2DFitSigma{Float64}(
        1.0, 1.0, 100.0, 5.0, 0.03, 0.01, 0.01, 0.0, 20.0, 0.5, 0.002, NaN, 1, 1, 0, 4
    )
    kept_pv, n_dropped_pv = SMLMAnalysis._drop_nonfinite_emitters([e_nan_pvalue])
    @test kept_pv == [e_nan_pvalue]
    @test n_dropped_pv == 0

    # Once kept, the fit overlay must colour it the way the pvalue filter treats
    # it (rejected), not green because both NaN comparisons are false.
    cfg_pv = FilterConfig(pvalue = (1.0e-3, 1.0))
    @test SMLMAnalysis._fit_box_color(e_nan_pvalue, cfg_pv, nothing) === :purple
    smld_pv = SMLMAnalysis.BasicSMLD([e_nan_pvalue], IdealCamera(16, 16, 0.1), 1, 1, Dict{String, Any}())
    @test isempty(SMLMAnalysis._filter_smld(smld_pv, cfg_pv).emitters)

    dfi = SMLMAnalysis.DensityFilterInfo(1000, 800, 5, 0.3)
    @test dfi.n_before == 1000
    @test dfi.threshold == 5

    # densityfilter's per-point query radius n_sigma*sqrt(σ_i^2+max_σ^2) is an
    # exact upper bound for the pair test dist < n_sigma*sqrt(σ_i^2+σ_j^2) it
    # gates on, so switching the KD-tree query from a single global radius to a
    # per-point one must not change which emitters survive. Cross-check against
    # an independent brute-force count of that same pair test.
    cam = IdealCamera(64, 64, 0.1)
    xs = [0.0, 0.05, 0.0, 0.05, 5.0]
    ys = [0.0, 0.0, 0.05, 0.05, 5.0]
    σs = [0.01, 0.01, 0.01, 0.05, 0.2]
    emitters = [SMLMAnalysis.Emitter2DFit(xs[i], ys[i], 1000.0, 10.0, σs[i], σs[i], 0.0, 1.0, 1.0, i, 1, 0, i) for i in 1:5]
    smld_df = SMLMAnalysis.BasicSMLD(emitters, cam, 1, 1, Dict{String, Any}())
    cfg_df = DensityFilterConfig(n_sigma = 3.0, min_neighbors = 1)
    filtered_df, _ = SMLMAnalysis.densityfilter_step(smld_df, cfg_df)

    σ = [sqrt(e.σ_x^2 + e.σ_y^2) for e in emitters]
    expected_counts = zeros(Int, 5)
    for i in 1:5, j in 1:5
        i == j && continue
        dist = sqrt((emitters[i].x - emitters[j].x)^2 + (emitters[i].y - emitters[j].y)^2)
        dist < cfg_df.n_sigma * sqrt(σ[i]^2 + σ[j]^2) && (expected_counts[i] += 1)
    end
    @test length(filtered_df.emitters) == count(>=(cfg_df.min_neighbors), expected_counts)

    # Non-finite σ is rejected up front rather than silently propagating.
    bad_emitters = [SMLMAnalysis.Emitter2DFit(0.0, 0.0, 1000.0, 10.0, NaN, 0.01, 0.0, 1.0, 1.0, 1, 1, 0, 1)]
    smld_bad = SMLMAnalysis.BasicSMLD(bad_emitters, cam, 1, 1, Dict{String, Any}())
    @test_throws ArgumentError SMLMAnalysis.densityfilter_step(smld_bad, cfg_df)

    # A finite σ_x that overflows once squared must throw too, not silently
    # drop both emitters as having no neighbours.
    huge_σ = [
        SMLMAnalysis.Emitter2DFit(x, 0.0, 1000.0, 10.0, 1.0e200, 0.01, 0.0, 1.0, 1.0, 1, 1, 0, i)
            for (i, x) in enumerate((0.0, 1.0e199))
    ]
    smld_huge_σ = SMLMAnalysis.BasicSMLD(huge_σ, cam, 1, 1, Dict{String, Any}())
    @test_throws ArgumentError SMLMAnalysis.densityfilter_step(smld_huge_σ, cfg_df)

    # IntensityFilter: non-finite photons rejected up front, at 100 emitters.
    # The finite check runs before the step's "too few emitters" shortcut, so
    # a 99-emitter input is rejected too (tested for a NaN x below).
    many = [SMLMAnalysis.Emitter2DFit(0.01i, 0.01i, 1000.0, 10.0, 0.01, 0.01, 0.0, 1.0, 1.0, i, 1, 0, i) for i in 2:100]
    bad_photon = SMLMAnalysis.Emitter2DFit(0.5, 0.5, NaN, 10.0, 0.01, 0.01, 0.0, 1.0, 1.0, 1, 1, 0, 1)
    smld_ifbad = SMLMAnalysis.BasicSMLD(vcat([bad_photon], many), cam, 1, 1, Dict{String, Any}())
    @test_throws ArgumentError SMLMAnalysis.intensityfilter_step(smld_ifbad, IntensityFilterConfig())

    # Regression: a NaN x (not just a NaN σ/photons) must throw ArgumentError
    # up front, not InexactError from floor(Int, NaN) deep in _bin_index /
    # the KD-tree query.
    bad_x = SMLMAnalysis.Emitter2DFit(NaN, 0.0, 1000.0, 10.0, 0.01, 0.01, 0.0, 1.0, 1.0, 1, 1, 0, 1)
    smld_df_bad_x = SMLMAnalysis.BasicSMLD([bad_x], cam, 1, 1, Dict{String, Any}())
    @test_throws ArgumentError SMLMAnalysis.densityfilter_step(smld_df_bad_x, cfg_df)

    bad_x_if = SMLMAnalysis.Emitter2DFit(NaN, 0.5, 1000.0, 10.0, 0.01, 0.01, 0.0, 1.0, 1.0, 1, 1, 0, 1)
    smld_if_bad_x = SMLMAnalysis.BasicSMLD(vcat([bad_x_if], many), cam, 1, 1, Dict{String, Any}())
    @test_throws ArgumentError SMLMAnalysis.intensityfilter_step(smld_if_bad_x, IntensityFilterConfig())
    # Below the 100-emitter "too few to estimate the field" shortcut as well.
    smld_if_bad_x_99 = SMLMAnalysis.BasicSMLD(vcat([bad_x_if], many[1:98]), cam, 1, 1, Dict{String, Any}())
    @test_throws ArgumentError SMLMAnalysis.intensityfilter_step(smld_if_bad_x_99, IntensityFilterConfig())

    # psf_sigma: an explicit (lo, hi) — including a 0.0 lower bound — is
    # always applied, unlike :auto's "skip when degenerate" behavior.
    mkem_sigma(σ, i) = GaussMLE.Emitter2DFitSigma{Float64}(
        0.0, 0.0, 1000.0, 5.0, σ, 0.01, 0.01, 0.0, 20.0, 0.5, 0.002, 1.0, 1, 1, 0, i
    )
    smld_ps = SMLMAnalysis.BasicSMLD([mkem_sigma(0.05, 1), mkem_sigma(0.3, 2)], cam, 1, 1, Dict{String, Any}())
    filtered_ps, _ = SMLMAnalysis.filter_step(smld_ps, FilterConfig(psf_sigma = (0.0, 0.1)))
    @test length(filtered_ps.emitters) == 1
    @test filtered_ps.emitters[1].σ == 0.05
end
