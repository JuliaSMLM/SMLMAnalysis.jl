using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Random
using TOML
using Statistics

@testset "upstream API smoke (tiny data)" begin
    # Runs the full detect/fit -> filter -> frame-connect -> drift -> render
    # pipeline through analyze(), at default verbosity so the figure/stats
    # writers (CairoMakie) run too, on data tiny enough to stay in the fast
    # tier. Its purpose is compat detection: each upstream API gets exercised
    # cheaply on every CI run, not just in the local thorough tier.
    Random.seed!(1)
    cam = IdealCamera(32, 32, 0.1)
    sim = SMLMAnalysis.StaticSMLMConfig(density = 5.0, σ_psf = 0.13, nframes = 50, ndatasets = 2)
    (_, si) = SMLMAnalysis.simulate(
        sim;
        pattern = SMLMAnalysis.Nmer2D(n = 8, d = 0.05),
        molecule = SMLMAnalysis.GenericFluor(photons = 5.0e4, k_off = 20.0, k_on = 0.04),
        camera = cam
    )
    images = [
        SMLMAnalysis.gen_images(
            si.smld_model, SMLMAnalysis.MicroscopePSFs.GaussianPSF(0.13);
            dataset = d, bg = 20.0, poisson_noise = true
        )[1] for d in 1:2
    ]

    cfg = AnalysisConfig(
        DetectFitConfig(
            boxer = BoxerConfig(boxsize = 7, psf_sigma = 0.13, backend = :cpu),
            fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu)
        ),
        FilterConfig(photons = (100.0, Inf)),
        FrameConnectConfig(max_frame_gap = 2),
        DriftConfig(degree = 1),
        RenderConfig(zoom = 5);
        camera = cam,
        outdir = mktempdir(),
    )
    t = @elapsed (result, info) = analyze(images, cfg)
    @info "upstream API smoke (tiny data) wall time" seconds = t

    @test result isa SMLMAnalysis.AnalysisResult
    @test length(result.smld.emitters) >= 1

    # stepinfo: by name, by config type, and the not-found ArgumentError.
    @test stepinfo(info, :driftcorrect).name == "driftcorrect"
    @test stepinfo(info, DriftConfig).name == "driftcorrect"
    @test_throws ArgumentError stepinfo(info, :nosuchstep)
    @test_throws ArgumentError stepinfo(info, SMLMAnalysis.CrossAlignConfig)
end
