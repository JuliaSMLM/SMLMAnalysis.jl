"""
Generate hero image for README.

Simulates a dense field of Nmer2D octamers on a small camera,
runs the full pipeline, and renders a cropped region.
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "examples"))

using SMLMAnalysis
using MicroscopePSFs
using Statistics

# Small camera, dense structures
camera = IdealCamera(64, 64, 0.1)  # 6.4x6.4 μm
psf_sigma = 0.13

sim = StaticSMLMConfig(
    density = 6.0,
    σ_psf = psf_sigma,
    nframes = 8000,
    ndatasets = 1,
)
pattern = Nmer2D(n=8, d=0.10)   # 100nm diameter octamers -- easier to resolve
fluor = GenericFluor(photons=50000.0, k_off=20.0, k_on=0.04)

println("Simulating...")
(_, sim_info) = simulate(sim; pattern=pattern, molecule=fluor, camera=camera)
smld_model = sim_info.smld_model

psf = MicroscopePSFs.GaussianPSF(psf_sigma)
println("Generating images...")
function gen_images_for_dataset(smld, psf, dataset::Int; kwargs...)
    emitters_d = filter(e -> e.dataset == dataset, smld.emitters)
    smld_d = BasicSMLD(emitters_d, smld.camera, smld.n_frames, 1, smld.metadata)
    (images, _) = gen_images(smld_d, psf; dataset=1, kwargs...)
    images
end
images = gen_images_for_dataset(smld_model, psf, 1; bg=20.0, poisson_noise=true)
println("Images: $(size(images))")

# Run pipeline
println("Running pipeline...")
(result, info) = analyze([images], AnalysisConfig(
    camera = camera,
    steps = [
        DetectFitConfig(
            boxer=BoxerConfig(boxsize=7, min_photons=500.0, psf_sigma=psf_sigma),
            fitter=GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20)),
        FilterConfig(photons=(500.0, Inf), precision=(0.0, 0.010), pvalue=(1e-3, 1.0), psf_sigma=:auto),
        FrameConnectConfig(max_frame_gap=5, calibration=CalibrationConfig(clamp_k_to_one=true)),
        DriftConfig(degree=1, dataset_mode=:continuous),
    ],
    verbose=Verbosity.PROGRESS,
))

smld = result.smld
println("Localizations: $(length(smld.emitters))")

# Crop to center ~3x3 μm for visible structures
cx, cy = 3.2, 3.2  # center of 6.4 μm FOV
hw = 1.5  # half-width
roi_emitters = filter(e -> abs(e.x - cx) < hw && abs(e.y - cy) < hw, smld.emitters)
roi_smld = BasicSMLD(roi_emitters, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)
println("ROI emitters: $(length(roi_emitters))")

outdir = joinpath(@__DIR__, "src", "assets")
mkpath(outdir)

(img, rinfo) = render(roi_smld, RenderConfig(
    pixel_size=2.0,
    colormap=:inferno,
    clip_percentile=0.95,
    scalebar=true,
    scalebar_color=:white,
    filename=joinpath(outdir, "render_gaussian.png"),
))
println("Rendered: $(rinfo.output_size)")
println("Saved: $(joinpath(outdir, "render_gaussian.png"))")
