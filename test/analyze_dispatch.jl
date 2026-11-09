using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Random
using TOML
using Statistics
using Logging

struct _SelectiveLogger <: Logging.AbstractLogger
    io::IOBuffer
end
Logging.min_enabled_level(::_SelectiveLogger) = Logging.Info
Logging.shouldlog(::_SelectiveLogger, level, _module, group, id) = id !== :hidden
Logging.catch_exceptions(::_SelectiveLogger) = false
Logging.handle_message(logger::_SelectiveLogger, level, message, args...; kwargs...) = println(logger.io, message)

@testset "TeeLogger honors each sink's shouldlog" begin
    buffer = IOBuffer()
    mktempdir() do outdir
        with_logger(_SelectiveLogger(buffer)) do
            SMLMAnalysis._with_log_file(outdir) do
                @info "visible message" _id = :visible
                @info "hidden message" _id = :hidden
            end
        end
        @test String(take!(buffer)) == "visible message\n"
        log = read(joinpath(outdir, "log.txt"), String)
        @test occursin("visible message", log)
        @test occursin("hidden message", log)
    end
end

@testset "TeeLogger file sink records all admitted messages" begin
    for level in (Logging.Debug, Logging.Error)
        buffer = IOBuffer()
        mktempdir() do outdir
            with_logger(ConsoleLogger(buffer, level)) do
                SMLMAnalysis._with_log_file(outdir) do
                    @debug "admitted debug message"
                    @info "admitted info message"
                end
            end
            log = read(joinpath(outdir, "log.txt"), String)
            caller_log = String(take!(buffer))
            @test occursin("admitted info message", log)
            @test occursin("admitted debug message", log) == (level == Logging.Debug)
            @test occursin("admitted info message", caller_log) == (level == Logging.Debug)
            @test occursin("admitted debug message", caller_log) == (level == Logging.Debug)
        end
    end
end

@testset "analyze dispatch" begin
    # Verify analyze() dispatch methods exist for each step config type
    @test hasmethod(analyze, Tuple{Vector{<:AbstractArray{<:Real, 3}}, DetectFitConfig})
    @test hasmethod(analyze, Tuple{AbstractArray{<:Real, 3}, DetectFitConfig})
    @test hasmethod(analyze, Tuple{DetectFitConfig})
    @test hasmethod(analyze, Tuple{SMLMAnalysis.BasicSMLD, FilterConfig})
    @test hasmethod(analyze, Tuple{SMLMAnalysis.BasicSMLD, FrameConnectConfig})
    @test hasmethod(analyze, Tuple{SMLMAnalysis.BasicSMLD, DriftConfig})
    @test hasmethod(analyze, Tuple{SMLMAnalysis.BasicSMLD, DensityFilterConfig})
    @test hasmethod(analyze, Tuple{SMLMAnalysis.BasicSMLD, RenderConfig})

    # analyze(data, config::AnalysisConfig) with a data type _normalize_data
    # doesn't recognize (e.g. a String — data is never itself a file path;
    # use `nothing` with a file-based DetectFitConfig instead).
    @test_throws ArgumentError analyze("some/path.h5", AnalysisConfig(camera = IdealCamera(8, 8, 0.1)))

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
    @test DetectFitConfig().boxer.boxsize == 9    # default box is 9 pixels
    cam = IdealCamera(64, 64, 0.1)
    cfg2 = DetectFitConfig(camera = cam, boxer = BoxerConfig(boxsize = 7))
    @test cfg2.camera === cam
    @test cfg2.boxer.boxsize == 7
    @test cfg2.fitter.psf_model isa GaussianXYNBS

    # _inject_camera (AnalysisConfig pipeline path): pixel_size/qe on the
    # DetectFitConfig would be silently ignored (the pipeline camera always
    # wins), so it's rejected instead of accepted-and-dropped.
    @test SMLMAnalysis._inject_camera(cfg, cam).camera === cam    # no pixel_size/qe: fine
    @test_throws ArgumentError SMLMAnalysis._inject_camera(DetectFitConfig(pixel_size = 0.1), cam)
    @test_throws ArgumentError SMLMAnalysis._inject_camera(DetectFitConfig(qe = 0.9), cam)
    @test_throws ArgumentError DetectFitConfig(camera = cam, pixel_size = 0.1)
    @test SMLMAnalysis._inject_camera(cfg2, cam) === cfg2  # camera already set: unchanged

    # DetectFitConfig.datasets selection field
    @test cfg.datasets === nothing                          # default is no selection
    cfg_range = DetectFitConfig(datasets = 1:19)
    @test cfg_range.datasets == 1:19
    @test cfg_range.datasets isa Vector{Int}   # concrete field type: any AbstractVector{Int} is accepted but stored as Vector{Int}
    cfg_sparse = DetectFitConfig(datasets = [1, 2, 3, 5, 7])
    @test cfg_sparse.datasets == [1, 2, 3, 5, 7]
    @test cfg_sparse.datasets isa Vector{Int}

    # _select_sources: pass-through when nothing, bounds-checked otherwise
    src = [(i = j,) for j in 1:5]
    @test SMLMAnalysis._select_sources(src, nothing) === src
    @test SMLMAnalysis._select_sources(src, [1, 3, 5]) == [src[1], src[3], src[5]]
    @test SMLMAnalysis._select_sources(src, 2:4) == src[2:4]
    @test_throws ArgumentError SMLMAnalysis._select_sources(src, [1, 6])
    @test_throws ArgumentError SMLMAnalysis._select_sources(src, [0, 1])

    # No kwargs... catch-all on step analyze() methods: a misspelled keyword
    # must raise MethodError, not silently vanish into a kwargs sink.
    smld_empty = SMLMAnalysis.BasicSMLD(SMLMAnalysis.Emitter2DFit{Float64}[], cam, 1, 1, Dict{String, Any}())
    @test_throws MethodError analyze(smld_empty, FilterConfig(); bogus_kwarg = 1)
end

@testset "DetectFit blank datasets and terminal statistics checkpoint" begin
    Random.seed!(4)
    cam = IdealCamera(32, 32, 0.1)
    sim = SMLMAnalysis.StaticSMLMConfig(density = 10.0, σ_psf = 0.13, nframes = 20, ndatasets = 1)
    _, sim_info = SMLMAnalysis.simulate(
        sim;
        pattern = SMLMAnalysis.Nmer2D(n = 6, d = 0.2),
        molecule = SMLMAnalysis.GenericFluor(photons = 5.0e4, k_off = 20.0, k_on = 1.0), camera = cam
    )
    psf = SMLMAnalysis.MicroscopePSFs.GaussianPSF(0.13)
    imgs, _ = SMLMAnalysis.gen_images(sim_info.smld_model, psf; bg = 20.0, poisson_noise = true)
    empty_model = SMLMAnalysis.BasicSMLD(similar(sim_info.smld_model.emitters, 0), cam, 20, 1, Dict{String, Any}())
    blank, _ = SMLMAnalysis.gen_images(empty_model, psf; bg = 20.0, poisson_noise = true)
    cfg = DetectFitConfig(
        camera = cam,
        boxer = BoxerConfig(boxsize = 7, psf_sigma = 0.13, backend = :cpu),
        fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu)
    )

    @testset "Blank dataset $nonblank" for (stacks, nonblank) in (([imgs, blank], 1), ([blank, imgs], 2))
        mktempdir() do outdir
            # Default verbosity runs info, stats, overlay, cache and per-frame writers.
            smld, info = analyze(stacks, cfg; outdir = outdir)
            @test smld.n_datasets == info.info.n_datasets == 2
            @test !isempty(smld.emitters)
            @test all(e -> e.dataset == nonblank, smld.emitters)
            @test length(info.info.fit_info) == 2
            @test eltype(info.info.fit_info) === GaussMLE.GaussMLEFitInfo
            @test info.info.fit_info[3 - nonblank].n_fits == 0
            @test info.info.fit_info[3 - nonblank].device_id == info.info.fit_info[nonblank].device_id
            @test info.info.fit_info[nonblank].n_fits ==
                length(smld.emitters) + info.info.n_out_of_field + info.info.n_nonfinite
            dir = joinpath(outdir, "01_detectfit")
            saved_info = TOML.parsefile(joinpath(dir, "info.toml"))
            @test saved_info["fit_info_$(3 - nonblank)"]["n_fits"] == 0
            @test isfile(joinpath(dir, "detection_overlay.png"))
            @test isfile(joinpath(dir, "localizations_per_frame.png"))
            cache = SMLMAnalysis.load_cache(outdir, "detectfit_samples.jld2")
            @test cache !== nothing
            @test any(f -> f <= 20, cache["sample_original_frames"])
            @test any(f -> f > 20, cache["sample_original_frames"])
            roi_abs_frames = cache["sample_original_frames"][cache["sample_roi_frames"]]
            @test all(f -> (nonblank - 1) * 20 < f <= nonblank * 20, roi_abs_frames)
        end
    end
    @test_throws r"detectfit: no localizations found across all 2 dataset\(s\)" analyze(
        [blank, blank], cfg; verbose = Verbosity.SILENT
    )

    @testset "Pipeline file logging respects the caller's logger" begin
        buffer = IOBuffer()
        for logger in (NullLogger(), ConsoleLogger(buffer, Logging.Error))
            mktempdir() do outdir
                pipeline = AnalysisConfig(
                    camera = cam, steps = [cfg, FilterConfig(photons = (1000.0, Inf))],
                    outdir = outdir, verbose = Verbosity.STANDARD, checkpoint = Checkpoint.NONE
                )
                result, info = with_logger(logger) do
                    analyze(imgs, pipeline)
                end
                @test !isempty(result.smld.emitters)
                @test length(info.step_infos) == 2
                log = read(joinpath(outdir, "log.txt"), String)
                @test occursin("Info", log)
                @test occursin("detectfit", log)
                @test occursin("filter", log)
            end
        end
        @test isempty(String(take!(buffer)))
    end

    @testset "Checkpoint before terminal Hopkins" begin
        mktempdir() do outdir
            Random.seed!(54)
            pipeline = AnalysisConfig(
                camera = cam, steps = [cfg, FilterConfig(photons = (1000.0, Inf)), HopkinsConfig()],
                outdir = outdir, checkpoint = Checkpoint.END, verbose = Verbosity.SILENT
            )
            result, _ = analyze(imgs, pipeline)
            checkpoints = [
                joinpath(dir, f) for (dir, _, files) in walkdir(outdir) for f in files
                    if startswith(f, "smld_") && endswith(f, ".h5")
            ]
            @test checkpoints == [joinpath(outdir, "02_filter", "smld_filtered.h5")]
            if length(checkpoints) == 1
                saved = SMLMAnalysis.load_smld(only(checkpoints))
                @test saved.n_frames == result.smld.n_frames
                @test saved.n_datasets == result.smld.n_datasets
                @test length(saved.emitters) == length(result.smld.emitters)
                for field in (:x, :y, :photons, :frame, :dataset)
                    @test getproperty.(saved.emitters, field) ≈ getproperty.(result.smld.emitters, field)
                end
            end
        end
    end
end

struct _UndeclaredCheckpointConfig <: SMLMAnalysis.AbstractSMLMConfig end

# A user-defined step following workflows/extending.md that saves its output when the
# orchestrator's checkpoint level is ALL, as the built-in filters do.
struct _CustomWriterConfig <: SMLMAnalysis.AbstractSMLMConfig
    keep::Int
end
function SMLMAnalysis.analyze(
        smld::SMLMAnalysis.BasicSMLD, cfg::_CustomWriterConfig;
        outdir = nothing, step_number::Int = 0, checkpoint::Int = Checkpoint.EXPENSIVE, kwargs...
    )
    out = SMLMAnalysis.BasicSMLD(smld.emitters[1:cfg.keep], smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)
    dir = SMLMAnalysis.step_outdir(outdir, step_number, cfg)
    if dir !== nothing && checkpoint >= Checkpoint.ALL
        mkpath(dir)
        save_smld(joinpath(dir, "smld_custom.h5"), out)
    end
    return (out, SMLMAnalysis.StepInfo(step_number, cfg, 0.0, Dict{Symbol, Any}()))
end

@testset "Checkpoint.END reaches a custom checkpoint-writing step" begin
    cam = IdealCamera(16, 16, 0.1)
    emitters = [
        SMLMAnalysis.Emitter2DFit{Float64}(0.1i, 0.1i, 1000.0 + i, 10.0, 0.01, 0.01, 30.0, 1.0; frame = i)
            for i in 1:10
    ]
    smld = SMLMAnalysis.BasicSMLD(emitters, cam, 10, 1, Dict{String, Any}())
    smld_files(outdir) = sort([joinpath(d, f) for (d, _, fs) in walkdir(outdir) for f in fs if startswith(f, "smld_")])
    for steps in (
            SMLMAnalysis.AbstractSMLMConfig[_CustomWriterConfig(4)],
            SMLMAnalysis.AbstractSMLMConfig[FilterConfig(photons = (0.0, Inf)), _CustomWriterConfig(4)],
            SMLMAnalysis.AbstractSMLMConfig[FilterConfig(photons = (0.0, Inf)), _CustomWriterConfig(4), HopkinsConfig()],
        )
        mktempdir() do outdir
            result, _ = SMLMAnalysis._run_pipeline(smld, steps, cam, outdir, Verbosity.SILENT, Checkpoint.END)
            files = smld_files(outdir)
            # Only the custom step's checkpoint, holding what the custom step kept.
            @test basename.(files) == ["smld_custom.h5"]
            if length(files) == 1
                @test length(SMLMAnalysis.load_smld(only(files)).emitters) == 4
            end
            @test length(result.smld.emitters) == 4
        end
    end
end

@testset "SMLD checkpoint writers are declared" begin
    # Covers every _save_step_smld call site, including each upstream config family.
    for cfg in (
            DetectFitConfig(), FilterConfig(), IntensityFilterConfig(), DensityFilterConfig(),
            OuterPolygonConfig(), KdeValleyConfig(), DBSCANConfig(eps_nm = 50.0), HDBSCANConfig(),
            VoronoiConfig(), HierarchicalConfig(n_clusters = 2), FrameConnectConfig(), DriftConfig(), BaGoLConfig(),
        )
        @test SMLMAnalysis._produces_smld(cfg)
    end
    @test !SMLMAnalysis._produces_smld(HopkinsConfig())
    @test !SMLMAnalysis._produces_smld(VoronoiDensityConfig())
    @test !SMLMAnalysis._produces_smld(RenderConfig())
    @test !SMLMAnalysis._produces_smld(CompositeRenderConfig())
    # A user-defined step may write a checkpoint, so it stays eligible for Checkpoint.END.
    @test SMLMAnalysis._produces_smld(_UndeclaredCheckpointConfig())
    steps = SMLMAnalysis.AbstractSMLMConfig[FilterConfig(), _UndeclaredCheckpointConfig(), HopkinsConfig(), RenderConfig()]
    @test SMLMAnalysis._find_final_smld_step(steps) == 2
end

@testset "crop axis conventions" begin
    # crop_images(imgs, roi_x, roi_y) == imgs[roi_y, roi_x, :]: roi_x indexes
    # columns (x, dim 2), roi_y indexes rows (y, dim 1). Encode (row, col) into
    # each pixel so an axis swap is caught. Locks a convention documented but
    # otherwise untested (SMART transposes on load, MIC does not, overlays
    # transpose before drawing — all easy to get backwards).
    nrow, ncol, nfr = 6, 8, 3
    img = [1000r + 10c + f for r in 1:nrow, c in 1:ncol, f in 1:nfr]
    roi_x, roi_y = 3:6, 2:4        # columns, rows
    cropped = SMLMAnalysis.crop_images(img, roi_x, roi_y)
    @test size(cropped) == (length(roi_y), length(roi_x), nfr)
    @test cropped == img[roi_y, roi_x, :]
    @test cropped[1, 1, 1] == 1000 * first(roi_y) + 10 * first(roi_x) + 1

    # crop_camera uses the same convention: roi_x → x-edges, roi_y → y-edges.
    cam = IdealCamera(ncol, nrow, 0.1)   # IdealCamera(nx=cols, ny=rows, px)
    cc = SMLMAnalysis.crop_camera(cam, roi_x, roi_y)
    @test cc.pixel_edges_x == cam.pixel_edges_x[first(roi_x):(last(roi_x) + 1)]
    @test cc.pixel_edges_y == cam.pixel_edges_y[first(roi_y):(last(roi_y) + 1)]
    @test length(cc.pixel_edges_x) - 1 == length(roi_x)   # x pixel count = #cols
    @test length(cc.pixel_edges_y) - 1 == length(roi_y)   # y pixel count = #rows
end

@testset "roi refuses a step-level DetectFit camera" begin
    # roi crops only the pipeline camera; a camera on DetectFitConfig would stay
    # full-frame and offset every localization by the crop origin.
    cam = IdealCamera(16, 16, 0.1)
    cfg = AnalysisConfig(
        camera = cam, roi = (x = 3:10, y = 3:10),
        steps = [DetectFitConfig(camera = cam)], outdir = nothing
    )
    # Match the message: all-zero frames also throw an (unrelated) ArgumentError downstream.
    @test_throws r"DetectFitConfig has its own camera" analyze(zeros(Float32, 16, 16, 2), cfg)
end
