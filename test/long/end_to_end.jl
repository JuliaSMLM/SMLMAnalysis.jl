using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Random
using TOML
using Statistics

@testset "end-to-end analyze() tuple contract" begin
    # Runs the full detect/fit → filter → render pipeline on synthetic
    # CPU-only data and checks the (result, info) contract the whole
    # package rests on. Mirrors the precompile workload (known-good),
    # but as an assertable test rather than a build-time smoke run.
    Random.seed!(1)
    cam = IdealCamera(32, 32, 0.1)
    sim = SMLMAnalysis.StaticSMLMConfig(density = 5.0, σ_psf = 0.13, nframes = 50, ndatasets = 1)
    (_, si) = SMLMAnalysis.simulate(sim;
        pattern  = SMLMAnalysis.Nmer2D(n = 8, d = 0.05),
        molecule = SMLMAnalysis.GenericFluor(photons = 5.0e4, k_off = 20.0, k_on = 0.04),
        camera   = cam)
    (imgs, _) = SMLMAnalysis.gen_images(si.smld_model, SMLMAnalysis.MicroscopePSFs.GaussianPSF(0.13);
        dataset = 1, bg = 20.0, poisson_noise = true)

    cfg = AnalysisConfig(
        DetectFitConfig(boxer  = BoxerConfig(boxsize = 7, psf_sigma = 0.13, backend = :cpu),
                        fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu)),
        FilterConfig(photons = (100.0, Inf)),
        RenderConfig(zoom = 10);
        camera  = cam,
        verbose = Verbosity.SILENT,
    )
    (result, info) = analyze([imgs], cfg)

    @test result isa SMLMAnalysis.AnalysisResult
    @test info isa SMLMAnalysis.AnalysisInfo
    @test result.smld isa SMLMAnalysis.BasicSMLD
    @test length(result.smld.emitters) > 0
    @test length(info.step_infos) == 3
    @test info.elapsed_s >= 0

    # detectfit narrows the accumulated AbstractEmitter[] to a concrete type;
    # the default GaussianXYNBS fitter yields Emitter2DFitSigma.
    @test isconcretetype(eltype(result.smld.emitters))
    @test eltype(result.smld.emitters) <: GaussMLE.Emitter2DFitSigma

    # The fitted output round-trips through HDF5 losslessly — concrete type and
    # PSF-width σ preserved (σ was silently lost when the SMLD was AbstractEmitter-typed).
    mktempdir() do dir
        p = joinpath(dir, "pipeline.h5")
        save_smld(p, result.smld)
        reloaded = load_smld(p)
        @test length(reloaded.emitters) == length(result.smld.emitters)
        @test eltype(reloaded.emitters) <: GaussMLE.Emitter2DFitSigma
        @test reloaded.emitters[1].σ ≈ result.smld.emitters[1].σ
        @test reloaded.emitters[1].σ_xy ≈ result.smld.emitters[1].σ_xy
    end

    # Multi-dataset detectfit: exercises the per-dataset loop (dataset field,
    # per-dataset frame numbering, n_datasets tracking) that the in-memory and
    # file-based paths share. Reuse the same stack as two datasets.
    (smld2, si2) = analyze([imgs, imgs],
        DetectFitConfig(camera = cam,
                        boxer  = BoxerConfig(boxsize = 7, psf_sigma = 0.13, backend = :cpu),
                        fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu));
        verbose = Verbosity.SILENT)
    @test si2.info.n_datasets == 2
    @test smld2.n_datasets == 2
    @test smld2.n_frames == size(imgs, 3)          # equal-length → per-dataset count
    @test Set(e.dataset for e in smld2.emitters) == Set([1, 2])
    @test all(1 <= e.frame <= size(imgs, 3) for e in smld2.emitters)  # frames are per-dataset
end

@testset "end-to-end with outdir" begin
    # Moderate simulated data through the full pipeline with disk output,
    # checking the on-disk contract: numbered step directories, valid TOML
    # provenance, every saved SMLD loadable and non-empty, and a render PNG.
    Random.seed!(1)
    cam = IdealCamera(64, 64, 0.1)
    sim = SMLMAnalysis.StaticSMLMConfig(density = 5.0, σ_psf = 0.13, nframes = 500, ndatasets = 2)
    (_, si) = SMLMAnalysis.simulate(sim;
        pattern  = SMLMAnalysis.Nmer2D(n = 8, d = 0.05),
        molecule = SMLMAnalysis.GenericFluor(photons = 5.0e4, k_off = 20.0, k_on = 0.04),
        camera   = cam)
    images = [SMLMAnalysis.gen_images(si.smld_model, SMLMAnalysis.MicroscopePSFs.GaussianPSF(0.13);
                          dataset = d, bg = 20.0, poisson_noise = true)[1] for d in 1:2]

    outdir = mktempdir()
    cfg = AnalysisConfig(
        DetectFitConfig(boxer  = BoxerConfig(boxsize = 7, psf_sigma = 0.13, backend = :cpu),
                        fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu)),
        FilterConfig(photons = (100.0, Inf)),
        FrameConnectConfig(max_frame_gap = 2),
        DriftConfig(degree = 1),
        RenderConfig(zoom = 5);
        camera = cam,
        outdir = outdir,
    )
    (result, info) = analyze(images, cfg)

    # Numbered step directories only -- outdir also holds a top-level
    # .cache/ (checkpoint cache, see common.jl cache_dir) that isn't a step.
    step_dirs = filter(f -> isdir(f) && occursin(r"^\d\d_", basename(f)),
                        readdir(outdir; join=true))
    @test length(step_dirs) == length(cfg.steps)

    tomls = String[]
    h5s = String[]
    pngs = String[]
    for (root, _, files) in walkdir(outdir)
        for f in files
            endswith(f, ".toml") && push!(tomls, joinpath(root, f))
            endswith(f, ".h5")   && push!(h5s, joinpath(root, f))
            endswith(f, ".png")  && push!(pngs, joinpath(root, f))
        end
    end
    @test !isempty(tomls)
    for p in tomls
        @test TOML.parsefile(p) isa Dict
    end

    smld_h5s = filter(p -> occursin("smld_", basename(p)), h5s)
    @test !isempty(smld_h5s)
    last_n = 0
    for p in smld_h5s
        loaded = load_smld(p)
        @test length(loaded.emitters) > 0
        if occursin("smld_corrected.h5", p)
            last_n = length(loaded.emitters)
        end
    end
    @test last_n > 0
    @test length(result.smld.emitters) == last_n

    @test !isempty(pngs)
end

@testset "multi-target orchestrator writes smld_<label>.h5 through _finalize_channels!" begin
    # A multi-target channel is always raw images (or a file path) that goes
    # through its own DetectFitConfig -- analyze(channels, MultiTargetConfig)
    # has no path for already-localized data (steps=[] leaves _run_pipeline's
    # state a raw image Vector, never a SMLMAnalysis.BasicSMLD, so it errors "Pipeline
    # produced no SMLD"). So this drives a tiny simulated 2-channel run
    # end-to-end -- the multi-target analogue of "upstream API smoke" above --
    # to catch a regression where the orchestrator stops calling
    # _finalize_channels! after the phase-2 multi-target steps (Codex #51).
    # Channel B is channel A's own pattern (same seed) offset by a known
    # 0.1 μm, so the saved file must show the *aligned* B, not raw B.
    cam = IdealCamera(32, 32, 0.1)
    dx_true = 0.1
    gen_channel(seed, dx, dy) = begin
        Random.seed!(seed)
        sim = SMLMAnalysis.StaticSMLMConfig(density = 5.0, σ_psf = 0.13, nframes = 50, ndatasets = 1)
        (_, si) = SMLMAnalysis.simulate(sim;
            pattern  = SMLMAnalysis.Nmer2D(n = 8, d = 0.05),
            molecule = SMLMAnalysis.GenericFluor(photons = 5.0e4, k_off = 20.0, k_on = 0.04),
            camera   = cam)
        model = deepcopy(si.smld_model)
        for e in model.emitters
            e.x += dx
            e.y += dy
        end
        (imgs, _) = SMLMAnalysis.gen_images(model, SMLMAnalysis.MicroscopePSFs.GaussianPSF(0.13);
            dataset = 1, bg = 20.0, poisson_noise = true)
        [imgs]
    end
    images_a = gen_channel(11, 0.0, 0.0)
    images_b = gen_channel(11, dx_true, 0.0)

    chan_cfg() = AnalysisConfig(
        DetectFitConfig(boxer  = BoxerConfig(boxsize = 7, psf_sigma = 0.13, backend = :cpu),
                        fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu)),
        FilterConfig(photons = (100.0, Inf)),
        FrameConnectConfig(max_frame_gap = 2);
        camera = cam,
    )

    outdir = mktempdir()
    mt = MultiTargetConfig(
        labels = [:A, :B],
        steps  = [CrossAlignConfig(align = AlignConfig(method = :fft))],
        outdir = outdir,
    )
    (result, info) = analyze([(images_a, chan_cfg()), (images_b, chan_cfg())], mt)

    @test result isa SMLMAnalysis.MultiTargetResult
    for label in (:A, :B)
        p = joinpath(outdir, "smld_$(label).h5")
        @test isfile(p)
        loaded = load_smld(p)
        @test [e.x for e in loaded.emitters] == [e.x for e in result[label].smld.emitters]
        @test [e.y for e in loaded.emitters] == [e.y for e in result[label].smld.emitters]
    end

    # smld_B.h5 == result.smlds[2] == result[:B].smld, all three agreeing on
    # the SAME (aligned) coordinates.
    loaded_b = load_smld(joinpath(outdir, "smld_B.h5"))
    @test [e.x for e in loaded_b.emitters] == [e.x for e in result.smlds[2].emitters]
    @test [e.y for e in loaded_b.emitters] == [e.y for e in result.smlds[2].emitters]
    @test result.smlds[2] === result[:B].smld

    # And the aligned B must actually differ from B's pre-alignment data
    # -- catches a regression where the "aligned" save silently falls
    # back to writing unaligned data. `result[:B].smld_connected` is
    # NOT a safe baseline for this: it's FrameConnectInfo's pre-COMBINE
    # linked data, and combining alone (independent of cross-align)
    # shifts the mean by its own ~0.2 μm, which would make this
    # assertion pass even with a broken cross-align. Instead, compare
    # against the channel's own last-saved pre-alignment SMLD -- the
    # frame-connect step's `smld_combined.h5`, i.e. B's phase-1 result
    # exactly as it entered CrossAlignConfig. Found by directory pattern
    # (not a hard-coded step number), since the channel's own step
    # numbering is an implementation detail of chan_cfg() above.
    b_dir = joinpath(outdir, "B")
    fc_dirs = filter(d -> occursin(r"^\d+_frameconnect$", d), readdir(b_dir))
    @test length(fc_dirs) == 1
    pre = load_smld(joinpath(b_dir, only(fc_dirs), "smld_combined.h5"))
    mean_x(s) = sum(e.x for e in s.emitters) / length(s.emitters)
    # dx_true = 0.1 μm was removed by alignment; require most of it back.
    @test abs(mean_x(result[:B].smld) - mean_x(pre)) >= 0.05
end
