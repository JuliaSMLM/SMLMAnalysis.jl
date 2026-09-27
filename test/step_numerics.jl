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
