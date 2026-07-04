"""
    generate_data.jl

Shared data generation and caching for SMLMAnalysis examples.

Generates simulated SMLM datasets and saves them to `examples/data/` (gitignored)
so that examples can skip the slow image generation step on subsequent runs.

Usage:
  - As a standalone script: `julia --project=. generate_data.jl` to generate all datasets
  - As an include: `include("generate_data.jl")` then `load_or_generate("single_target")`
"""

using SMLMAnalysis
using SMLMSim: Nanoruler2D   # new pattern, not yet re-exported by SMLMAnalysis
using MicroscopePSFs
using JLD2

const DATA_DIR = joinpath(@__DIR__, "data")

# ============================================================================
# Helper: Generate images per dataset (workaround for SMLMSim bug)
# ============================================================================

"""
Generate images for a specific dataset by filtering emitters first.
Workaround for gen_images not respecting dataset parameter.
"""
function gen_images_for_dataset(smld, psf, dataset::Int; kwargs...)
    emitters_d = filter(e -> e.dataset == dataset, smld.emitters)
    smld_d = BasicSMLD(emitters_d, smld.camera, smld.n_frames, 1, smld.metadata)
    (images, _) = gen_images(smld_d, psf; dataset=1, kwargs...)
    images
end

# ============================================================================
# Single-target dataset generation
# ============================================================================

function generate_single_target()
    println("Generating single-target dataset...")

    camera = IdealCamera(256, 128, 0.1)
    n_frames = 2000
    n_datasets = 4
    psf_sigma = 0.13

    sim_params = StaticSMLMConfig(
        density = 2.0,
        σ_psf = psf_sigma,
        nframes = n_frames,
        ndatasets = n_datasets,
    )
    pattern = Nmer2D(n=8, d=0.05)
    fluor = GenericFluor(photons=50000.0, k_off=20.0, k_on=0.04)

    t = @elapsed begin
        (_, sim_info) = simulate(sim_params; pattern=pattern, molecule=fluor, camera=camera)
        smld_model = sim_info.smld_model

        psf = MicroscopePSFs.GaussianPSF(psf_sigma)
        image_stacks = [gen_images_for_dataset(smld_model, psf, d; bg=20.0, poisson_noise=true)
                        for d in 1:n_datasets]
        images = cat(image_stacks...; dims=3)
    end
    println("  $(n_datasets) datasets x $(n_frames) frames, $(size(images)) ($(round(t, digits=1))s)")

    data = Dict{String,Any}(
        "images" => images,
        "camera_nx" => 256,
        "camera_ny" => 128,
        "camera_pixelsize" => 0.1,
        "n_frames" => n_frames,
        "n_datasets" => n_datasets,
        "psf_sigma" => psf_sigma,
    )

    mkpath(DATA_DIR)
    path = joinpath(DATA_DIR, "single_target.jld2")
    jldopen(path, "w") do f
        for (k, v) in data
            f[k] = v
        end
    end
    println("  Saved: $path ($(round(filesize(path) / 1e6, digits=1)) MB)")

    data
end

# ============================================================================
# Line pattern dataset generation
# ============================================================================

function generate_lines()
    println("Generating line pattern dataset...")

    camera = IdealCamera(256, 128, 0.1)
    n_frames = 2000
    n_datasets = 4
    psf_sigma = 0.13

    sim_params = StaticSMLMConfig(
        density = 5.0,
        σ_psf = psf_sigma,
        nframes = n_frames,
        ndatasets = n_datasets,
    )
    pattern = Line2D(λ=50.0, endpoints=[(-0.4, 0.0), (0.4, 0.0)])
    fluor = GenericFluor(photons=50000.0, k_off=20.0, k_on=0.04)

    t = @elapsed begin
        (_, sim_info) = simulate(sim_params; pattern=pattern, molecule=fluor, camera=camera)
        smld_model = sim_info.smld_model

        psf = MicroscopePSFs.GaussianPSF(psf_sigma)
        image_stacks = [gen_images_for_dataset(smld_model, psf, d; bg=20.0, poisson_noise=true)
                        for d in 1:n_datasets]
        images = cat(image_stacks...; dims=3)
    end
    println("  $(n_datasets) datasets x $(n_frames) frames, $(size(images)) ($(round(t, digits=1))s)")

    data = Dict{String,Any}(
        "images" => images,
        "camera_nx" => 256,
        "camera_ny" => 128,
        "camera_pixelsize" => 0.1,
        "n_frames" => n_frames,
        "n_datasets" => n_datasets,
        "psf_sigma" => psf_sigma,
    )

    mkpath(DATA_DIR)
    path = joinpath(DATA_DIR, "lines.jld2")
    jldopen(path, "w") do f
        for (k, v) in data
            f[k] = v
        end
    end
    println("  Saved: $path ($(round(filesize(path) / 1e6, digits=1)) MB)")

    data
end

# ============================================================================
# Nanoruler pattern dataset generation
# ============================================================================

function generate_nanoruler()
    println("Generating nanoruler pattern dataset...")

    # Acquisition targets (realistic dSTORM-like blinking), translated to
    # SMLMSim's 2-state CTMC fluorophore model (q = [-k_off k_off; k_on -k_on],
    # state 1 = ON emitting at γ Hz) at framerate F = 50 Hz, T = nframes/F = 40 s:
    #   - mean blink (ON) lifetime 2 frames  -> k_off = F/2      = 25 Hz
    #   - 1000 photons / saturated frame     -> γ    = 1000·F    = 50000 Hz
    #   - ~5 blinks (ON-periods) per emitter -> k_on ≈ 0.126 Hz  (with state1=:equilibrium)
    # k_on sets blinks-per-emitter and is INDEPENDENT of pattern density. SMLMSim
    # verified k_on=0.126 + state1=:equilibrium gives a median of exactly 5 ON-periods
    # over the 2000-frame movie (the original 0.04 gave only ~1.6).
    # NOTE on blinks vs localizations: one ON-period (blink) spans ~3 frames above
    # minphotons, so ~5 blinks -> ~13 localizations per mark per movie.
    #
    # DENSITY: 2 rulers/μm² — sparse and well-separated (rulers ~700 nm apart vs
    # ~80 nm long), so the 3-mark structure stays cleanly isolated. Per-frame
    # activation is then very sparse (~0.03 active/μm²/frame, an OUTPUT not a
    # target), so within-frame PSF blends are negligible and single-emitter fitting
    # is clean. (A full 1 active/μm²/frame would instead need ~45 rulers/μm², which
    # packs the rulers into an unresolvable overlapping field — a different, denser
    # regime where blend rejection craters the yield; not what we want here.)
    camera = IdealCamera(256, 128, 0.1)
    n_frames = 2000
    n_datasets = 4
    framerate = 50.0
    psf_sigma = 0.13

    sim_params = StaticSMLMConfig(
        density = 2.0,            # rulers per μm² (sparse, well-separated)
        σ_psf = psf_sigma,
        nframes = n_frames,
        ndatasets = n_datasets,
        framerate = framerate,
    )
    # GATTAquant-style 3-mark ruler, 40 nm between adjacent marks, random
    # orientation/position per ruler (handled by SMLMSim's uniform2D placement).
    pattern = Nanoruler2D(n=3, spacing=0.04)
    fluor = GenericFluor(photons=50000.0, k_off=25.0, k_on=0.126)

    t = @elapsed begin
        # state1=:equilibrium (simulate's default) starts each emitter in its
        # steady-state on/off split so the ~5-blinks count is unbiased.
        (_, sim_info) = simulate(sim_params; pattern=pattern, molecule=fluor,
                                 camera=camera, state1=:equilibrium)
        smld_model = sim_info.smld_model

        psf = MicroscopePSFs.GaussianPSF(psf_sigma)
        image_stacks = [gen_images_for_dataset(smld_model, psf, d; bg=20.0, poisson_noise=true)
                        for d in 1:n_datasets]
        images = cat(image_stacks...; dims=3)
    end
    println("  $(n_datasets) datasets x $(n_frames) frames, $(size(images)) ($(round(t, digits=1))s)")

    data = Dict{String,Any}(
        "images" => images,
        "camera_nx" => 256,
        "camera_ny" => 128,
        "camera_pixelsize" => 0.1,
        "n_frames" => n_frames,
        "n_datasets" => n_datasets,
        "psf_sigma" => psf_sigma,
        "ruler_spacing" => 0.04,
        "ruler_marks" => 3,
    )

    mkpath(DATA_DIR)
    path = joinpath(DATA_DIR, "nanoruler.jld2")
    jldopen(path, "w") do f
        for (k, v) in data
            f[k] = v
        end
    end
    println("  Saved: $path ($(round(filesize(path) / 1e6, digits=1)) MB)")

    data
end

# ============================================================================
# Load or generate helper
# ============================================================================

"""
    load_or_generate(name::String; force=false) -> Dict{String,Any}

Load cached dataset from `examples/data/{name}.jld2`, or generate and cache it.

Returns a Dict with dataset-specific keys. Key additions:
- `"image_stacks"`: Vector of views into concatenated images (1 per dataset)
- `"images"`: Concatenated 3D array (kept for backward compat)

The `image_stacks` field encodes dataset boundaries in the data structure,
eliminating the need for `n_datasets` as a separate parameter.
"""
function load_or_generate(name::String; force=false)
    path = joinpath(DATA_DIR, "$name.jld2")

    if !force && isfile(path)
        println("Loading cached data: $path")
        t = @elapsed data = Dict{String,Any}(String(k) => v for (k, v) in pairs(load(path)))
        println("  Loaded ($(round(t, digits=1))s)")
    else
        if name == "single_target"
            data = generate_single_target()
        elseif name == "lines"
            data = generate_lines()
        elseif name == "nanoruler"
            data = generate_nanoruler()
        else
            error("Unknown dataset: $name. Available: single_target, lines, nanoruler")
        end
    end

    # Create image_stacks as Vector of views (dataset boundaries from data structure)
    n_ds = data["n_datasets"]
    n_fr = data["n_frames"]
    imgs = data["images"]
    data["image_stacks"] = [@view imgs[:, :, (d-1)*n_fr+1:d*n_fr] for d in 1:n_ds]

    return data
end

# ============================================================================
# Standalone: generate all datasets
# ============================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    println("="^60)
    println("Generating all cached datasets")
    println("="^60)
    println()

    t_total = @elapsed begin
        generate_single_target()
        println()
        generate_lines()
        println()
        generate_nanoruler()
    end

    println()
    println("="^60)
    println("Done ($(round(t_total, digits=1))s total)")
    println("="^60)
end
