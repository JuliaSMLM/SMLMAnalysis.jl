using SMLMAnalysis
using GaussMLE
using Test
using TOML
using Logging

@testset "DetectFit physical field boundaries" begin
    parent = IdealCamera(64, 48, (0.125, 0.25))
    maps = reshape(collect(1.0:(64 * 48)), 48, 64)
    scmos = SMLMAnalysis.SCMOSCamera(
        parent.pixel_edges_x, parent.pixel_edges_y, maps, 2maps, 3maps, 0.9
    )
    for parent_camera in (parent, scmos)
        camera = SMLMAnalysis.crop_camera(parent_camera, 17:56, 9:40)
        xmin, xmax = extrema(camera.pixel_edges_x)
        ymin, ymax = extrema(camera.pixel_edges_y)
        xmid, ymid = (xmin + xmax) / 2, (ymin + ymax) / 2
        positions = [
            (xmin, ymin), (xmax, ymax), (xmin, ymax), (xmax, ymin), (xmid, ymid),
            (prevfloat(xmin), ymid), (nextfloat(xmax), ymid),
            (xmid, prevfloat(ymin)), (xmid, nextfloat(ymax)),
        ]
        emitters = [
            SMLMAnalysis.Emitter2DFit{Float64}(x, y, 1000.0, 20.0, 0.01, 0.02, 10.0, 1.0; id = i)
                for (i, (x, y)) in enumerate(positions)
        ]
        original = copy(emitters)
        @test SMLMAnalysis._drop_out_of_field_emitters!(emitters, camera) == 4
        @test length(emitters) == 5
        @test all(emitters[i] === original[i] for i in eachindex(emitters))
        if camera isa SMLMAnalysis.SCMOSCamera
            @test camera.offset == maps[9:40, 17:56]
            @test camera.gain == (2maps)[9:40, 17:56]
            @test camera.readnoise == (3maps)[9:40, 17:56]
        end
    end
end

# Deterministic spot images. The off-camera Gaussian has its brightest sampled
# pixel on the last row, so detection places its fitting box against that edge.
_field_stack(y0) = [
    round(UInt16, 20 + 20000 * exp(-((x - 20)^2 + (y - y0)^2) / (2 * 1.3^2)))
        for y in 1:32, x in 1:40, f in 1:2
]
_field_values(emitters) = [Tuple(getfield(e, f) for f in fieldnames(typeof(e))) for e in emitters]

@testset "DetectFit removes edge fits and preserves interior fits" begin
    boxer = BoxerConfig(boxsize = 9, psf_sigma = 0.13, backend = :cpu)
    fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), iterations = 50, backend = :cpu)
    edge, interior = _field_stack(35), _field_stack(12)
    full_camera = IdealCamera(70, 52, 0.1)
    cropped_camera = SMLMAnalysis.crop_camera(full_camera, 17:56, 9:40)

    for camera in (IdealCamera(40, 32, 0.1), cropped_camera)
        cfg = DetectFitConfig(camera = camera, boxer = boxer, fitter = fitter)
        in_field(e) = first(camera.pixel_edges_x) <= e.x <= last(camera.pixel_edges_x) &&
            first(camera.pixel_edges_y) <= e.y <= last(camera.pixel_edges_y)

        edge_boxes, _ = SMLMAnalysis.SMLMBoxer.getboxes(edge, camera, boxer)
        @test length(edge_boxes) == 2
        @test all(==(32 - boxer.boxsize + 1), edge_boxes.y_corners)
        raw_edge, _ = GaussMLE.fit(edge_boxes, fitter)
        @test all(e -> e.y > last(camera.pixel_edges_y), raw_edge.emitters)
        @test all(e -> isfinite(e.x) && isfinite(e.y), raw_edge.emitters)

        interior_boxes, _ = SMLMAnalysis.SMLMBoxer.getboxes(interior, camera, boxer)
        raw_interior, _ = GaussMLE.fit(interior_boxes, fitter)
        @test !isempty(raw_interior.emitters)
        @test all(in_field, raw_interior.emitters)
        kept, info = analyze(interior, cfg; verbose = Verbosity.SILENT)
        @test isequal(_field_values(kept.emitters), _field_values(raw_interior.emitters))
        @test info.info.n_out_of_field == info.summary[:n_out_of_field] == 0

        log = IOBuffer()
        mixed, mixed_info = with_logger(ConsoleLogger(log)) do
            analyze([edge, interior, edge], cfg; verbose = Verbosity.PROGRESS)
        end
        @test all(in_field, mixed.emitters)
        @test all(e -> e.dataset == 2, mixed.emitters)
        @test mixed_info.info.n_out_of_field == 2length(raw_edge.emitters)
        @test mixed_info.info.n_fits == 2length(raw_edge.emitters) + length(raw_interior.emitters)
        dataset_logs = filter(l -> occursin("ROIs ->", l), split(String(take!(log)), '\n'))
        @test length(dataset_logs) == 3
        @test occursin("2 out-of-field fits removed", dataset_logs[1])
        @test !occursin("removed", dataset_logs[2])
        @test occursin("2 out-of-field fits removed", dataset_logs[3])
        printed = sprint(show, MIME"text/plain"(), SMLMAnalysis.AnalysisInfo(0.0, [mixed_info]))
        @test occursin("n_out_of_field: 4", printed)

        mktempdir() do outdir
            # All fits can be removed; info.toml, stats.md and the checkpoint
            # must still write, including with the cropped camera's origin.
            empty_smld, empty_info = analyze(edge, cfg; outdir = outdir)
            @test isempty(empty_smld.emitters)
            @test empty_info.info.n_out_of_field == length(raw_edge.emitters)
            dir = joinpath(outdir, "01_detectfit")
            saved = TOML.parsefile(joinpath(dir, "info.toml"))
            @test saved["n_out_of_field"] == empty_info.info.n_out_of_field
            @test saved["fit_info"]["n_fits"] == length(raw_edge.emitters)
            stats = read(joinpath(dir, "stats.md"), String)
            @test occursin("**Out-of-field fits removed**: 2", stats)
            @test occursin("No localizations remain", stats)
            @test isempty(load_smld(joinpath(dir, "smld_raw.h5")).emitters)
        end

        mktempdir() do dir
            path = joinpath(dir, "field.h5")
            SMLMAnalysis.HDF5.h5open(path, "w") do f
                g = SMLMAnalysis.HDF5.create_group(f, "Main")
                g["data"] = cat(edge, interior, dims = 3)
                SMLMAnalysis.HDF5.attrs(g["data"])["dimension_order"] = "HWN"
            end
            file_cfg = DetectFitConfig(
                camera = camera, boxer = boxer, fitter = fitter, path = path,
                h5_format = :smart, dataset_frames = [1:2, 3:4]
            )
            from_file, file_info = analyze(file_cfg; verbose = Verbosity.SILENT)
            from_memory, memory_info = analyze([edge, interior], cfg; verbose = Verbosity.SILENT)
            @test all(in_field, from_file.emitters)
            @test isequal(_field_values(from_file.emitters), _field_values(from_memory.emitters))
            @test file_info.info.n_out_of_field == memory_info.info.n_out_of_field == 2
        end
    end
end

@testset "DetectFit warns when cleanup removes every fit" begin
    camera = IdealCamera(40, 32, 0.1)
    boxer = BoxerConfig(boxsize = 9, psf_sigma = 0.13, backend = :cpu)
    cases = (
        (images = _field_stack(35), model = GaussianXYNBS(), n_nonfinite = 0, n_outside = 2),
        # Zero PSF width leaves finite, in-field positions but NaN uncertainties.
        (images = _field_stack(12), model = GaussianXYNB(0.0f0), n_nonfinite = 2, n_outside = 0),
    )
    for case in cases
        cfg = DetectFitConfig(
            camera = camera, boxer = boxer,
            fitter = GaussMLEConfig(psf_model = case.model, iterations = 50, backend = :cpu)
        )
        warning = "detectfit: every fit was removed during cleanup ($(case.n_nonfinite) non-finite, " *
            "$(case.n_outside) outside the camera); returning an empty result."
        for write_outputs in (false, true)
            mktempdir() do tmp
                outdir = write_outputs ? tmp : nothing
                logger = Test.TestLogger(min_level = Logging.Warn)
                smld, info = with_logger(logger) do
                    analyze(case.images, cfg; outdir = outdir)
                end
                cleanup_warnings = filter(log -> occursin("every fit was removed", log.message), logger.logs)
                @test length(cleanup_warnings) == 1
                @test only(cleanup_warnings).level == Logging.Warn
                @test only(cleanup_warnings).message == warning
                @test isempty(smld.emitters)
                @test info.info.n_fits == 2
                @test info.info.n_nonfinite == case.n_nonfinite
                @test info.info.n_out_of_field == case.n_outside
                if write_outputs
                    dir = joinpath(outdir, "01_detectfit")
                    saved = TOML.parsefile(joinpath(dir, "info.toml"))
                    @test saved["n_nonfinite"] == case.n_nonfinite
                    @test saved["n_out_of_field"] == case.n_outside
                    @test occursin("No localizations remain after fit cleanup.", read(joinpath(dir, "stats.md"), String))
                    @test isempty(load_smld(joinpath(dir, "smld_raw.h5")).emitters)
                end
            end
        end
    end
end

@testset "DetectFitInfo constructor defaults" begin
    kwargs = (
        boxes_info = SMLMAnalysis.SMLMBoxer.BoxesInfo[], fit_info = GaussMLE.GaussMLEFitInfo[],
        n_datasets = 1, n_rois = 2, n_fits = 2, n_frames_per_dataset = 2, elapsed_s = 0.1,
    )
    @test SMLMAnalysis.DetectFitInfo(; kwargs...).n_out_of_field == 0
    @test SMLMAnalysis.DetectFitInfo(; kwargs..., n_out_of_field = 2).n_out_of_field == 2
    args = values(kwargs)
    for trailing in ((), (nothing,), (nothing, 1))
        @test SMLMAnalysis.DetectFitInfo(args..., trailing...).n_out_of_field == 0
    end
end
