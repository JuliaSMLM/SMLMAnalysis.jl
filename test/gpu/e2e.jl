# SMLMAnalysis must load before PSFLearning: SMLMAnalysis loads CUDNN_jll (through SMLMBoxer),
# which fails to initialise once Reactant's bundled cuDNN is loaded.
using SMLMAnalysis
using SMLMDeepFit: SMLMDeepFit
using PSFLearning: PSFLearning
using MicroscopePSFs: MicroscopePSFs
using SMLMSim: SMLMSim
using Test
using Random
using Statistics

# Integration plan section 4, E2E: PSFLearning learns a PSF from a simulated bead stack and
# writes psf.h5; GaussMLE (DetectFitConfig.psf_file) and DeepFit (trained on psf.h5) then
# localize a movie simulated with the true PSF. Julia 1.12, CUDA, one process, sequential.

# MicroscopePSFs.save_psf tags the file with the module-qualified type name, which load_psf
# cannot read back (MicroscopePSFs 0.5.7 fixes it); PSFLearning.save_psf rewrites the tag
# the same way.
function _save_psf_file(path, psf)
    MicroscopePSFs.save_psf(path, psf)
    SMLMAnalysis.HDF5.h5open(path, "r+") do f
        a = SMLMAnalysis.HDF5.attrs(f)
        a["psf_type"] = last(split(a["psf_type"], '.'))
    end
    return path
end

# Second-moment widths (µm) of a PSF at depth z, over ±1 µm.
function _widths(psf, z)
    xs = -1.0:0.02:1.0
    img = [psf(x, y, z) for y in xs, x in xs]   # rows y, columns x
    s = sum(img)
    mx, my = sum(img .* xs') / s, sum(img .* xs) / s
    return sqrt(sum(img .* (xs' .- mx) .^ 2) / s), sqrt(sum(img .* (xs .- my) .^ 2) / s)
end

# One-to-one matching per frame: pairs within `lat` µm laterally and `ax` µm axially,
# closest first. Returns the matched (truth, fit) index pairs.
function _match(truth, fits; lat = 0.2, ax = 0.4)
    byframe = Dict{Int, Vector{Int}}()
    for (k, e) in enumerate(fits)
        push!(get!(byframe, e.frame, Int[]), k)
    end
    pairs = Tuple{Float64, Int, Int}[]
    for (i, t) in enumerate(truth), k in get(byframe, t.frame, Int[])
        d = hypot(fits[k].x - t.x, fits[k].y - t.y)
        d <= lat && abs(fits[k].z - t.z) <= ax && push!(pairs, (d, i, k))
    end
    sort!(pairs)
    usedt, usedf, out = falses(length(truth)), falses(length(fits)), Tuple{Int, Int}[]
    for (_, i, k) in pairs
        usedt[i] || usedf[k] || (usedt[i] = usedf[k] = true; push!(out, (i, k)))
    end
    return out
end

# Errors (nm) of the matched fits for the truth emitters selected by `keep`.
function _errors(truth, fits, pairs, keep)
    sel = [(i, k) for (i, k) in pairs if keep(truth[i])]
    dx = [(fits[k].x - truth[i].x) * 1000 for (i, k) in sel]
    dy = [(fits[k].y - truth[i].y) * 1000 for (i, k) in sel]
    dz = [(fits[k].z - truth[i].z) * 1000 for (i, k) in sel]
    return dx, dy, dz, length(sel) / count(keep, truth)
end
_rmse(v) = sqrt(mean(abs2, v))

@testset "E2E: PSFLearning -> psf.h5 -> GaussMLE and DeepFit" begin
    mktempdir() do dir
        px = 0.1
        # 1. Ground truth: vector model, astigmatism (Noll 6) and spherical (Noll 11).
        cfg = PSFLearning.PSFLearningConfig(
            model = :vector, pixelsize_x = Float32(px), pixelsize_y = Float32(px),
            roi_size = 21, backend = :reactant_gpu
        )
        fwd = PSFLearning.build_forward(cfg)
        gt = zeros(Float32, PSFLearning.n_pupil_coeffs(fwd))
        gt[6] = -0.5f0   # rad; this sign makes the PSF wider in x at z = +0.5 µm
        gt[11] = 0.1f0
        truth_psf = PSFLearning.to_psf(cfg, gt; law = :bead, z_range = (-0.8, 0.8))
        (sx, sy) = _widths(truth_psf, 0.5)
        @info "E2E truth PSF widths at z = +0.5 µm" sx sy
        @test sx / sy >= 1.3   # fails if the astigmatism is on the wrong mode (Noll 5: 1.0)

        # 2. Bead stack from PSFLearning's forward model: 20000 photons per plane, bg 20,
        # Poisson noise, background subtracted as PSFLearning's bead pipeline does.
        z_positions = collect(Float32, -0.8:0.1:0.8)
        raw = PSFLearning.forward_images(fwd, gt, z_positions)   # (N_z, M, M), rows y
        Random.seed!(21)   # poisson_noise draws from the global RNG
        counts = SMLMSim.poisson_noise(Float64.(raw) .* (20000 / sum(raw[9, :, :])) .+ 20.0)
        stack = Float32.(max.(counts .- 20.0, 0.0))
        recovered(c) = abs(c[6] + 0.5) <= 0.05 && abs(c[11] - 0.1) <= 0.05   # rad

        # 3. Learn through SMLMAnalysis's psflearning step, which writes psf.h5.
        # PSFLearning's bead-stack loss has no intensity parameter, so from this photon-scale
        # stack learn_psf converges to wrong coefficients (Noll 6 about -0.30 for -0.50).
        # PSFLearning owns the fix; this flips to an unexpected pass when it lands.
        (_, photon_info) = analyze(stack, cfg; z_positions = z_positions)
        pc = photon_info.info.coeffs
        @info "E2E learn from the photon-scale stack (known broken)" pc[6] pc[11]
        @test_broken recovered(pc)
        # Until then the arms below learn from the stack rescaled to the forward model's own
        # scale, which measures them against a correctly learned PSF.
        scale0 = sum(PSFLearning.forward_images(fwd, zero(gt), Float32[0])) /
            maximum(sum(stack; dims = (2, 3)))
        outdir = joinpath(dir, "out")
        (_, learn_info) = analyze(
            stack .* Float32(scale0), cfg; z_positions = z_positions, outdir = outdir, step_number = 1
        )
        psf_path = joinpath(SMLMAnalysis.step_outdir(outdir, 1, cfg), "psf.h5")
        @test isfile(psf_path)
        learned = MicroscopePSFs.load_psf(psf_path)
        @test learned isa MicroscopePSFs.SplinePSF
        @test learned.z_min <= -0.8 && learned.z_max >= 0.8
        coeffs = learn_info.info.coeffs
        @info "E2E learned coefficients (truth: Noll 6 = -0.5, Noll 11 = 0.1)" coeffs[6] coeffs[11] seconds =
            learn_info.elapsed_s
        @test recovered(coeffs)

        # 4. Movie with the true PSF: 64x64 pixels of 0.1 µm, 9 emitters per frame on a
        # jittered 3x3 grid (at least 1.7 µm apart), 3000 photons, bg 10, z uniform in ±0.5 µm.
        rng = Xoshiro(41)
        n_frames = 500
        cam = IdealCamera(64, 64, px)
        truth = SMLMAnalysis.Emitter3DFit{Float64}[]
        for f in 1:n_frames, cx in (1.1, 3.2, 5.3), cy in (1.1, 3.2, 5.3)
            x = cx + 0.4 * rand(rng) - 0.2
            y = cy + 0.4 * rand(rng) - 0.2
            z = rand(rng) - 0.5
            push!(truth, SMLMAnalysis.Emitter3DFit{Float64}(x, y, z, 3000.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0; frame = f))
        end
        smld_true = SMLMAnalysis.BasicSMLD(truth, cam, n_frames, 1, Dict{String, Any}())
        Random.seed!(42)
        (movie, _) = SMLMAnalysis.gen_images(smld_true, truth_psf; bg = 10.0, poisson_noise = true, support = 1.0)

        # 5. GaussMLE arm: DetectFit with psf_file, and a control with the true PSF saved the
        # same way, over the whole ±0.5 µm; the 0.4-0.5 µm band is also reported on its own.
        # Detection psf_sigma = 0.22 µm is the true PSF's equivalent Gaussian sigma (the sigma
        # with its peak-to-integral ratio) averaged over the tested ±0.5 µm: 0.16 µm in focus,
        # 0.34 µm at z = +0.5 µm. SMLMBoxer's photon estimate assumes that Gaussian; with the
        # in-focus width (0.15) about half of the z > +0.4 µm emitters gave no ROI, in both arms.
        control_path = _save_psf_file(joinpath(dir, "truth_psf.h5"), truth_psf)
        gmle = Dict{String, Any}()
        for (arm, path) in (("learned", psf_path), ("control", control_path))
            dcfg = DetectFitConfig(
                boxer = BoxerConfig(boxsize = 15, psf_sigma = 0.22, backend = :cpu),
                fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu),
                camera = cam, psf_file = path
            )
            (smld, _) = analyze(movie, dcfg; outdir = joinpath(dir, "gaussmle_$arm"))
            fits = smld.emitters
            pairs = _match(truth, fits)
            gmle[arm] = (
                accept = _errors(truth, fits, pairs, t -> true),
                band = _errors(truth, fits, pairs, t -> abs(t.z) >= 0.4),
            )
            (dx, dy, dz, recall) = gmle[arm].accept
            (bx, by, bz, brecall) = gmle[arm].band
            @info "E2E GaussMLE $arm arm" recall bias_nm = (mean(dx), mean(dy), mean(dz)) rmse_nm =
                (_rmse(dx), _rmse(dy), _rmse(dz)) band_recall = brecall band_rmse_nm = (_rmse(bx), _rmse(by), _rmse(bz))
        end
        (dx, dy, dz, recall) = gmle["learned"].accept
        (cx, cy, cz, _) = gmle["control"].accept
        @test recall >= 0.95
        @test abs(mean(dx)) <= 3
        @test abs(mean(dy)) <= 3
        @test abs(mean(dz)) <= 15
        @test _rmse(dx) <= min(15, 1.25 * _rmse(cx))
        @test _rmse(dy) <= min(15, 1.25 * _rmse(cy))
        @test _rmse(dz) <= min(50, 1.25 * _rmse(cz))

        # 6. DeepFit arm: train DECODE on psf.h5 (settings of SMLMDeepFit's
        # examples/train_decode.jl, photons and bg of the movie), then infer on the movie with
        # the 100-count pedestal the training data carry. The 100 epochs are that file's demo
        # length, so this arm shows integration and learning, not production quality.
        dfdir = joinpath(dir, "deepfit")
        decode = SMLMDeepFit.Decode(;
            sz = 40, ρ = 1.0, photons = 3000.0, bg = 10.0, minz = -0.5, maxz = 0.5, bgmaxz = 0.8,
            pixelsize = px, n_train = 2000, n_test = 200, psffile = psf_path, savepath = dfdir
        )
        tcfg = SMLMDeepFit.TrainConfig(;
            traintype = decode, epochs = 100, batchsize = 50,
            optimiser = SMLMDeepFit.Optimisers.Adam(1.0f-4), patience = 200, data_refresh_interval = 50,
            seed = 42, use_reactant = SMLMDeepFit.REACTANT_AVAILABLE[], use_cuda = true, savepath = dfdir, infotime = 5,
            checktime = 50, tblogger = false
        )
        (result, train_info) = analyze(tcfg; outdir = outdir, step_number = 2)
        @test isfile(result.model_path)
        @info "E2E DeepFit training" best_epoch = train_info.info.best_epoch seconds = train_info.elapsed_s

        icfg = SMLMDeepFit.DeepFitConfig(; model_path = result.model_path, ccdoffset = 0.0f0, camera = cam)
        (smld, _) = analyze(movie .+ 100, icfg; outdir = outdir, step_number = 3)
        # DeepFit's window n covers frames n..n+2 and predicts the centre frame, n+1, but
        # SMLMDeepFit stamps its fits with n. Fits sit in their own frame when shift 0 matches
        # best. Until then the arm is scored at shift +1; when this flips to an unexpected pass,
        # SMLMDeepFit's fix has landed: score at shift 0 (fits1 = shifted(0)) and drop this test.
        shifted(s) = [(x = e.x, y = e.y, z = e.z, frame = e.frame + s) for e in smld.emitters]
        jac(t, f) = (p = _match(t, f); length(p) / (length(t) + length(f) - length(p)))
        jshift = Dict(s => jac(truth, shifted(s)) for s in (-1, 0, 1))
        @test_broken jshift[0] >= max(jshift[-1], jshift[1])
        fits1 = shifted(1)
        # The first and last frames are no window's centre, so they never get fits.
        seen = [t for t in truth if 2 <= t.frame <= n_frames - 1]
        pairs = _match(seen, fits1)
        jaccard = length(pairs) / (length(seen) + length(fits1) - length(pairs))
        (dx, dy, dz, recall) = _errors(seen, fits1, pairs, t -> true)
        # Chance: the same fits against the truth 250 frames away, which shares only the grid.
        chance = jac(seen, [(; f..., frame = mod1(f.frame + 250, n_frames)) for f in fits1])
        # z matched laterally only (0.15 µm), so a wrong z is scored rather than dropped.
        zpairs = _match(seen, fits1; lat = 0.15, ax = Inf)
        zcorr = cor([seen[i].z for (i, _) in zpairs], [fits1[k].z for (_, k) in zpairs])
        @info "E2E DeepFit arm (fits at frame +1)" n_fits = length(fits1) jaccard_shift_m1_0_p1 = (jshift[-1], jshift[0], jshift[1]) jaccard chance recall bias_nm =
            (mean(dx), mean(dy), mean(dz)) rmse_nm = (_rmse(dx), _rmse(dy), _rmse(dz)) zcorr n_zpairs = length(zpairs)
        # Production bars (Jaccard 0.7, RMSE 40/40/80 nm), which demo training is not expected to meet.
        @info "E2E DeepFit arm against production bars" jaccard_ge_0_7 = jaccard >= 0.7 rmse_le_40_40_80 =
            _rmse(dx) <= 40 && _rmse(dy) <= 40 && _rmse(dz) <= 80
        # Floors that show learning:
        # Jaccard: nine emitters per frame on a fixed jittered grid put a fit inside some truth's
        # match window by chance; chance above is that level on this movie. The floor sits above it.
        @test chance < 0.25
        @test jaccard >= 0.25
        # Lateral: a fit placed uniformly at random in the 0.2 µm match disc has RMSE
        # 0.2 / 2 = 100 nm per axis; 60 nm shows positions learned, not just landing in the window.
        @test _rmse(dx) <= 60
        @test _rmse(dy) <= 60
        # Axial: z with no information correlates 0 with the truth (sd 1/sqrt(n), n above) and
        # a flipped axis correlates negatively; 0.4 rules out both.
        @test zcorr >= 0.4
    end
end
