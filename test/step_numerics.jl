using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Random
using TOML
using Statistics

@testset "Bleaching fit degeneracy guard" begin
    Random.seed!(42)

    # Flat data: constant ~30 locs/frame with Poisson-like noise, 19k frames.
    # Reproduces a real-data case that produced a=-42873, k≈0 with unbounded NelderMead.
    # Expected: reject as degenerate, return nothing.
    flat = [max(0, round(Int, 30 + randn() * sqrt(30))) for _ in 1:19_000]
    @test SMLMAnalysis._estimate_bleaching_rate(flat) === nothing

    # Pure exponential decay: should recover parameters accurately (non-regression).
    t_exp = 1:5000
    exp_data = [max(0, round(Int, 10 + 50 * exp(-0.0005 * i) + randn() * 2)) for i in t_exp]
    res_exp = SMLMAnalysis._estimate_bleaching_rate(exp_data)
    @test res_exp !== nothing
    @test isapprox(res_exp.k_bleach, 5.0e-4, rtol = 0.1)
    @test isapprox(res_exp.offset, 10.0, atol = 2.0)
    @test isapprox(res_exp.N_0, 50.0, atol = 3.0)
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

@testset "Driftcorrect max drift is measured within one dataset" begin
    # Dataset 2 sits 5 um off (inter shift) and drifts with a P1 intra term running -0.1 to
    # +0.1 um over its frames. max_drift_nm measures from dataset 2's first frame, so the
    # excursion is 0.2 um = 200 nm and the 5 um inter shift never enters it; _calc_inter_shifts
    # reports the shift separately.
    m = SMLMDriftCorrection.LegendrePolynomial(2, 2, 100; degree = 1)
    m.inter[2].dm[1] = 5.0
    m.intra[2].dm[1].coefficients = [0.1]
    @test SMLMAnalysis._calc_max_drift(m, 100) ≈ 200.0 atol = 1.0e-6
    @test SMLMAnalysis._calc_inter_shifts(m)[2] ≈ 5000.0 atol = 1.0e-6

    # One dataset with traj(frame 1) = 0, as continuous-mode normalization pins it: the
    # since-first-frame value equals the old drift_trajectory excursion.
    m1 = SMLMDriftCorrection.LegendrePolynomial(2, 1, 100; degree = 1)
    m1.inter[1].dm[1] = 0.1
    m1.intra[1].dm[1].coefficients = [0.1]
    old = maximum(abs.(SMLMDriftCorrection.drift_trajectory(m1; dataset = 1).x)) * 1000
    @test old ≈ 200.0 atol = 1.0e-6
    @test SMLMAnalysis._calc_max_drift(m1, 100) ≈ old atol = 1.0e-6
end
