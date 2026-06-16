#
# dev/learn_psf_sim.jl — Stage 1: IN-SITU PSF learning (GPU).
#
# Learns the PSF from single-molecule emitters (the SMLM data itself) via
# PSFLearning.learn_psf_insitu_relearn: a joint L-BFGS fit of all emitters + the shared
# pupil, with iterative rejection of outliers (out-of-z / dim / worst-residual). Writes
# psf.h5 + diagnostics (kept-vs-rejected "ROIs picked from data" overlay, recovered PSF
# z-stack, Zernike recovery, stats).
#
#   julia -t auto --project=. dev/learn_psf_sim.jl
#
# SIM standin: synthetic single molecules (PSFLearning.forward_images_insitu) at random FOV
# positions + injected outliers. REAL data: pre-filter the SMLM localizations to ~1e3-1e4
# bright/isolated/good-CRLB candidates, extract an ROI per localization, pass initial_* from
# the fit. Loads PSFLearning WITHOUT SMLMAnalysis (cuDNN: SMLMBoxer→CUDNN_jll vs Reactant 9.14).
#
using PSFLearning
using HDF5
using CairoMakie
using Random

include(joinpath(@__DIR__, "deepfit_sim_common.jl"))
mkpath(joinpath(OUTDIR, "01_psflearning"))

Random.seed!(1)

function main()
    cfg = PSFLearning.PSFLearningConfig(
        n_max = 4, pupilsize = 32, roi_size = ROI,
        NA = 1.4f0, λ = 0.68f0, n_imm = 1.516f0,
        pixelsize_x = Float32(PIXEL_SIZE), pixelsize_y = Float32(PIXEL_SIZE),
        model = :scalar, iterations = 40, mu_init = 0.0f0, backend = :reactant_gpu,
    )
    m  = PSFLearning.build_insitu_model(cfg)
    nc = PSFLearning.n_coeffs(m)
    true_coeffs = zeros(Float32, nc); true_coeffs[5] = 0.4f0; true_coeffs[7] = -0.2f0   # defocus j5, coma j7

    # --- synthetic single-molecule candidate set: ~1000 good in-focus emitters + outliers ---
    N_good, N_bad = 1100, 120
    good_z  = clamp.(randn(Float32, N_good) .* 0.25f0, -0.5f0, 0.5f0)
    good_bg = fill(15f0, N_good); good_I = fill(4000f0, N_good)
    good = PSFLearning.forward_images_insitu(true_coeffs, zeros(Float32, N_good), zeros(Float32, N_good),
                                             good_z, good_bg, good_I, m)
    n_oz = N_bad ÷ 2; n_dim = N_bad - n_oz                                  # half out-of-z, half too dim
    bad_z  = vcat((rand(Float32, n_oz) .* 2f0 .+ 1f0) .* rand(Float32[-1f0, 1f0], n_oz),   # |z|∈[1,3] → out of z_range
                  clamp.(randn(Float32, n_dim) .* 0.25f0, -0.5f0, 0.5f0))
    bad_bg = fill(15f0, N_bad)
    bad_I  = vcat(fill(4000f0, n_oz), rand(Float32, n_dim) .* 50f0 .+ 20f0)                 # dim I∈[20,70] < min_photons
    bad = PSFLearning.forward_images_insitu(true_coeffs, zeros(Float32, N_bad), zeros(Float32, N_bad),
                                            bad_z, bad_bg, bad_I, m)
    data = Float32.(cat(good, bad; dims = 1))                              # (N, M, M); outliers = rows N_good+1..N
    N = N_good + N_bad
    z0 = vcat(good_z, bad_z); bg0 = vcat(good_bg, bad_bg); I0 = vcat(good_I, bad_I)
    fov_x = rand(Float32, N) .* (FOV * PIXEL_SIZE); fov_y = rand(Float32, N) .* (FOV * PIXEL_SIZE)   # FOV positions for the overlay

    # --- in-situ learn + iterative outlier rejection ---
    (psf, rinfo) = PSFLearning.learn_psf_insitu_relearn(data, cfg;
        relearn_iterations = 3, z_range = (-0.6f0, 0.6f0), min_photons = 100f0,
        reject_fraction = 0.05f0, min_emitters = 10,
        initial_z = z0, initial_bg = bg0, initial_intensity = I0)
    PSFLearning.save_psf(PSF_PATH, rinfo.info)
    HDF5.h5open(PSF_PATH, "r+") do f                                       # strip module-qualified psf_type
        a = HDF5.attributes(f)
        if haskey(a, "psf_type")
            t = read(a["psf_type"])
            occursin(".", t) && (HDF5.delete_attribute(f, "psf_type"); a["psf_type"] = String(split(t, ".")[end]))
        end
    end

    kept = rinfo.kept_indices; rejected = setdiff(1:N, kept); info = rinfo.info

    try
        dir = joinpath(OUTDIR, "01_psflearning")

        # "ROIs picked from data" overlay: emitters in the FOV, kept (used for the PSF) vs rejected
        figo = Figure(size = (560, 560))
        axo = Axis(figo[1, 1], xlabel = "x (µm)", ylabel = "y (µm)", aspect = DataAspect(),
                   title = "in-situ PSF learning: kept ($(length(kept))) vs rejected ($(length(rejected)))")
        scatter!(axo, fov_x[rejected], fov_y[rejected], color = :red, markersize = 12, marker = :xcross, label = "rejected")
        scatter!(axo, fov_x[kept], fov_y[kept], color = :limegreen, markersize = 9, label = "kept (used for PSF)")
        axislegend(axo)
        CairoMakie.save(joinpath(dir, "insitu_selection.png"), figo)

        # recovered PSF z-stack (learned model at fixed z-planes)
        zp = Float32[-0.5, -0.25, 0.0, 0.25, 0.5]
        pp = PSFLearning.PupilFieldParams(n_max = cfg.n_max, grid_size = cfg.pupilsize, NA = cfg.NA, λ = cfg.λ, n_imm = cfg.n_imm)
        cp = PSFLearning.prechirpz(cfg.pupilsize, cfg.roi_size, 2f0 * cfg.NA / (cfg.λ * cfg.pupilsize), cfg.pixelsize_x, cfg.pixelsize_y)
        model = Array(PSFLearning.forward_images(info.coeffs, pp, cp, zp))
        fig = Figure(size = (210 * length(zp), 260))
        Label(fig[0, 1:length(zp)], "in-situ learned PSF model across z", fontsize = 13)
        for k in 1:length(zp)
            a = Axis(fig[1, k], title = "z=$(zp[k])µm", aspect = DataAspect(), yreversed = true)
            heatmap!(a, model[k, :, :]', colormap = :inferno); hidedecorations!(a)
        end
        CairoMakie.save(joinpath(dir, "psf_zstack.png"), fig)

        # Zernike recovery
        truec = zeros(Float32, pp.n_zernike); truec[5] = 0.4f0; truec[7] = -0.2f0
        fig2 = Figure(size = (760, 340))
        ax = Axis(fig2[1, 1], xlabel = "Zernike index j", ylabel = "phase coeff (rad)",
                  title = "Learned Zernike phase (bars) vs true (red)")
        barplot!(ax, 1:pp.n_zernike, info.coeffs[1:pp.n_zernike], color = (:steelblue, 0.85))
        scatter!(ax, 1:pp.n_zernike, truec, color = :red, markersize = 12, marker = :hline)
        CairoMakie.save(joinpath(dir, "zernike_coeffs.png"), fig2)

        open(joinpath(dir, "stats.md"), "w") do io
            println(io, "# PSF Learning (in-situ)\n")
            println(io, "- model=scalar  NA=$(cfg.NA)  λ=$(cfg.λ)µm  pixelsize=$(cfg.pixelsize_x)µm  roi_size=$(cfg.roi_size)")
            println(io, "- candidate emitters=$N  kept=$(length(kept))  rejected=$(length(rejected))")
            println(io, "- relearn iterations=$(info.iterations)  converged=$(info.converged)\n")
            println(io, "## Zernike phase recovery (learned vs true)\n\n| j | learned | true |\n|---|---------|------|")
            for j in 1:min(pp.n_zernike, 12)
                println(io, "| $j | $(round(info.coeffs[j], digits = 3)) | $(round(truec[j], digits = 3)) |")
            end
        end
        @info "psflearning (in-situ) figures written" dir kept = length(kept) rejected = length(rejected)
    catch e
        @warn "psflearning figure generation failed (psf.h5 still valid)" exception = e
    end

    @info "stage 1 (psflearning in-situ, GPU) done" psf = typeof(psf) kept = length(rinfo.kept_indices) isfile(PSF_PATH)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
