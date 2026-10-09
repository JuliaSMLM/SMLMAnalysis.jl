using SMLMAnalysis
using GaussMLE
using MicroscopePSFs: MicroscopePSFs
using Test
using Random
using Statistics
using TOML

# MicroscopePSFs.save_psf tags the file with the module-qualified type name when the type
# is not visible from Main ("MicroscopePSFs.SplinePSF"), and load_psf cannot read that tag
# back. PSFLearning.save_psf rewrites the tag the same way; so does this helper.
function _save_psf_file(path, psf)
    MicroscopePSFs.save_psf(path, psf)
    SMLMAnalysis.HDF5.h5open(path, "r+") do f
        a = SMLMAnalysis.HDF5.attrs(f)
        a["psf_type"] = last(split(a["psf_type"], '.'))
    end
    return path
end

@testset "detectfit with a learned PSF file (psf_file)" begin
    # Integration plan S1: an astigmatic 3D SplinePSF saved as MicroscopePSFs HDF5,
    # data simulated with SMLMSim (pixel-integrated), then analyze() with psf_file.
    zc = MicroscopePSFs.ZernikeCoefficients(15)
    zc.phase[6] = 0.5   # vertical astigmatism, 0.5 rad RMS
    spline = MicroscopePSFs.SplinePSF(
        MicroscopePSFs.ScalarPSF(1.4, 0.6, 1.518; zernike_coeffs = zc);
        lateral_range = 1.0, axial_range = 0.8, lateral_step = 0.05, axial_step = 0.05
    )

    # 64x64 frames of 0.1 µm pixels; 9 emitters per frame on a jittered 3x3 grid, at
    # least 1.7 µm apart so 15-pixel boxes never overlap.
    # Frames 1-100 (the acceptance set): z uniform in ±0.5 µm.
    # Frames 101-140 (the outer band): |z| uniform in [0.4, 0.5] µm, either sign. Fits
    # there used to land on the wrong side of focus until GaussMLE started spline fits at
    # the best z of a likelihood scan (GaussMLE.jl#18); the band keeps that fixed.
    rng = Xoshiro(11)
    px, n_accept, n_frames = 0.1, 100, 140
    cam = IdealCamera(64, 64, px)
    centres = (1.1, 3.2, 5.3)
    truth = SMLMAnalysis.Emitter3DFit{Float64}[]
    for f in 1:n_frames, cx in centres, cy in centres
        x = cx + 0.4 * rand(rng) - 0.2
        y = cy + 0.4 * rand(rng) - 0.2
        z = f <= n_accept ? rand(rng) - 0.5 :
            (rand(rng, Bool) ? 1 : -1) * (0.4 + 0.1 * rand(rng))
        push!(
            truth,
            SMLMAnalysis.Emitter3DFit{Float64}(
                x, y, z, 3000.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0; frame = f
            )
        )
    end
    smld_true = SMLMAnalysis.BasicSMLD(truth, cam, n_frames, 1, Dict{String, Any}())
    Random.seed!(12)   # gen_images draws its Poisson noise from the global RNG
    (images, _) = SMLMAnalysis.gen_images(
        smld_true, spline; bg = 10.0, poisson_noise = true, support = 1.0
    )

    mktempdir() do dir
        psf_path = _save_psf_file(joinpath(dir, "psf.h5"), spline)
        cfg = DetectFitConfig(
            boxer = BoxerConfig(boxsize = 15, psf_sigma = 0.15, backend = :cpu),
            fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu),
            camera = cam,
            psf_file = psf_path,
        )
        outdir = joinpath(dir, "out")
        (smld, _) = analyze(images, cfg; outdir = outdir)
        fits = smld.emitters
        @test :z in fieldnames(eltype(fits))   # 3D fits (GaussMLE's Emitter3DFitGaussMLE)

        # Match each true emitter to the nearest fit in its frame within 200 nm laterally.
        by_frame = Dict{Int, Vector{Int}}()
        for (k, e) in enumerate(fits)
            push!(get!(by_frame, e.frame, Int[]), k)
        end
        # Errors (nm) of the matched fits, and the recall, over the true emitters `ts`.
        function errors(ts)
            dx, dy, dz = Float64[], Float64[], Float64[]
            for t in ts
                best, bestd = 0, Inf
                for k in get(by_frame, t.frame, Int[])
                    d = hypot(fits[k].x - t.x, fits[k].y - t.y)
                    d < bestd && ((best, bestd) = (k, d))
                end
                bestd <= 0.2 || continue
                push!(dx, (fits[best].x - t.x) * 1000)
                push!(dy, (fits[best].y - t.y) * 1000)
                push!(dz, (fits[best].z - t.z) * 1000)
            end
            return dx, dy, dz, length(dx) / length(ts)
        end
        rmse(v) = sqrt(mean(abs2, v))

        (dx, dy, dz, recall) = errors(filter(t -> t.frame <= n_accept, truth))
        @info "psf_file acceptance, z in ±0.5 µm" recall bias_nm = (mean(dx), mean(dy), mean(dz)) rmse_nm =
            (rmse(dx), rmse(dy), rmse(dz))
        @test recall >= 0.95
        @test abs(mean(dx)) <= 3
        @test abs(mean(dy)) <= 3
        @test abs(mean(dz)) <= 15
        @test rmse(dx) <= 15
        @test rmse(dy) <= 15
        @test rmse(dz) <= 50

        (bx, by, bz, brecall) = errors(filter(t -> t.frame > n_accept, truth))
        @info "psf_file outer band, |z| in [0.4, 0.5] µm" brecall bias_nm = (mean(bx), mean(by), mean(bz)) rmse_nm =
            (rmse(bx), rmse(by), rmse(bz))
        @test brecall >= 0.95
        @test rmse(bz) <= 50
        @test rmse(bx) <= 15 && rmse(by) <= 15

        # The step's config.toml records the file and the model's type, not the spline.
        toml_path = joinpath(SMLMAnalysis.step_outdir(outdir, 1, cfg), "config.toml")
        @test filesize(toml_path) < 10_000
        parsed = TOML.parsefile(toml_path)
        @test parsed["psf_file"] == psf_path
        @test parsed["fitter"]["psf_model"]["type"] == "SplinePSFModel"
        @test parsed["fitter"]["psf_model"]["pixel_size"] ≈ px
    end
end
