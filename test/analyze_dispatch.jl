using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Random
using TOML
using Statistics

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
    @test SMLMAnalysis._inject_camera(DetectFitConfig(camera = cam, pixel_size = 0.1), cam).camera === cam  # camera already set: unchanged

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
