"""
    workflows.jl

High-level workflow functions that chain multiple SMLM analysis steps together.
Each workflow function handles data transformations and optionally tracks
processing steps in an SMLMWorkflow object.
"""

using SMLMData
using SMLMSim
using SMLMBoxer
using GaussMLE
using MicroscopePSFs

"""
    simulate_detect_fit_workflow(sim_params, camera; detect_params=NamedTuple(),
                                  fit_params=NamedTuple(), workflow=nothing)
                                  → (smld_fitted, smld_ground_truth, workflow)

Complete SMLM localization workflow: simulation → detection → fitting.

This workflow:
1. Simulates SMLM data with ground truth localizations
2. Generates realistic camera images
3. Detects particles using difference of Gaussians
4. Fits Gaussian PSFs to get precise localizations

# Arguments
- `sim_params`: SMLMSim parameters (StaticSMLMParams or DiffusionSMLMParams)
- `camera`: Camera object (IdealCamera or SCMOSCamera) for image generation and detection

## Keyword Arguments
- `detect_params`: NamedTuple of parameters for `getboxes()` (boxsize, overlap, sigma_small, sigma_large, minval, use_gpu)
- `fit_params`: NamedTuple of parameters for `GaussMLEFitter()` (psf_model, iterations, device, batch_size)
- `workflow`: Optional SMLMWorkflow object to track processing steps

# Returns
Named tuple with:
- `smld_fitted`: BasicSMLD with fitted localizations
- `smld_ground_truth`: BasicSMLD with ground truth from simulation
- `images`: Generated camera images (for inspection)
- `detections`: Raw detection results from boxer
- `fit_results`: Raw GaussMLE fit results
- `workflow`: Updated workflow object (or new one if none provided)

# Example
```julia
using SMLMAnalysis, SMLMData, SMLMSim

# Setup simulation
camera = IdealCamera(1:256, 1:256, 0.1f0)  # 256×256, 100nm pixels
sim_params = StaticSMLMParams(
    density = 1.0,      # molecules per μm²
    σ_psf = 0.13,       # PSF width in μm
    photons = 1000,     # photons per emitter
    background = 10     # background photons
)

# Run complete workflow
result = simulate_detect_fit_workflow(sim_params, camera)

# Access results
println("Ground truth: ", length(result.smld_ground_truth.emitters), " emitters")
println("Detected: ", size(result.detections.coords_pixels, 1), " particles")
println("Fitted: ", length(result.smld_fitted.emitters), " localizations")

# Inspect workflow
println(result.workflow)
```
"""
function simulate_detect_fit_workflow(
    sim_params::SMLMSim.SMLMSimParams,
    camera::SMLMData.AbstractCamera;
    pattern = nothing,
    molecule = SMLMSim.Core.GenericFluor(photons=1e4, k_off=50.0, k_on=1e-2),
    detect_params::NamedTuple = NamedTuple(),
    fit_params::NamedTuple = NamedTuple(),
    workflow::Union{SMLMWorkflow,Nothing} = nothing
)
    # Create workflow tracker if not provided
    if workflow === nothing
        workflow = SMLMWorkflow("Simulate → Detect → Fit")
    end

    # Step 1: Simulate SMLM data
    simulated = simulate(sim_params; pattern=pattern, molecule=molecule, camera=camera)

    # Handle different simulation return types
    if simulated isa Tuple && length(simulated) == 3
        # StaticSMLM returns (pattern, smld_true, smld_noisy)
        smld_ground_truth = simulated[3]  # Use noisy version as "ground truth" with applied localization noise
    else
        # DiffusionSMLM returns just SMLD
        smld_ground_truth = simulated
    end

    add_step!(workflow, "Simulate SMLM Data", :simulate,
             Dict{Symbol,Any}(:sim_params => string(typeof(sim_params))),
             :SMLMSim, summarize_smld(smld_ground_truth))

    # Step 2: Generate camera images with noise
    psf = MicroscopePSFs.GaussianPSF(sim_params.σ_psf)
    images = gen_images(smld_ground_truth, psf;
                       bg=10.0,              # Background photons/pixel
                       poisson_noise=true,   # Add shot noise
                       camera_noise=false)   # IdealCamera (no read noise)

    add_step!(workflow, "Generate Camera Images", :gen_images,
             Dict{Symbol,Any}(:psf_sigma => sim_params.σ_psf, :n_frames => size(images, 3)),
             :SMLMSim, "Image stack size $(size(images))")

    # Step 3: Detect particles
    # Merge default params with user params
    default_detect = (
        boxsize = 11,
        overlap = 2.0,
        sigma_small = 1.0,
        sigma_large = 2.0,
        minval = 10.0,
        use_gpu = true
    )
    detect_params_full = merge(default_detect, detect_params)

    detections = getboxes(images, camera;
                         boxsize = detect_params_full.boxsize,
                         overlap = detect_params_full.overlap,
                         sigma_small = detect_params_full.sigma_small,
                         sigma_large = detect_params_full.sigma_large,
                         minval = detect_params_full.minval,
                         use_gpu = detect_params_full.use_gpu)

    # getboxes now returns ROIBatch directly
    roi_batch = detections

    add_step!(workflow, "Detect Particles", :getboxes,
             Dict{Symbol,Any}(pairs(detect_params_full)...),
             :SMLMBoxer, summarize_boxer_result(roi_batch))

    # Check if any detections were made
    if length(roi_batch) == 0
        @warn "No particles detected! Returning empty results. Check detection threshold (minval) or simulation parameters."
        empty_smld = SMLMData.BasicSMLD(
            SMLMData.Emitter2DFit{Float32}[],
            camera,
            sim_params.nframes,
            1,
            Dict{String,Any}("workflow" => "simulate_detect_fit", "warning" => "no detections")
        )
        return (
            smld_fitted = empty_smld,
            smld_ground_truth = smld_ground_truth,
            images = images,
            detections = roi_batch,
            fit_results = nothing,
            workflow = workflow
        )
    end

    # Step 5: Fit Gaussian PSFs
    # Merge default params with user params
    default_fit = (
        psf_model = GaussianXYNB(Float32(sim_params.σ_psf)),
        iterations = 20,
        device = nothing,  # auto-detect
        batch_size = 10_000
    )
    fit_params_full = merge(default_fit, fit_params)

    fitter = GaussMLEFitter(;
        psf_model = fit_params_full.psf_model,
        iterations = fit_params_full.iterations,
        device = fit_params_full.device,
        batch_size = fit_params_full.batch_size
    )

    # fit() now returns BasicSMLD directly
    smld_fitted = fit(fitter, roi_batch)

    add_step!(workflow, "Fit Gaussian PSFs", :fit,
             Dict{Symbol,Any}(:psf_model => string(typeof(fit_params_full.psf_model)),
                  :iterations => fit_params_full.iterations,
                  :device => string(fitter.device)),
             :GaussMLE, summarize_smld(smld_fitted))

    # Return complete results
    return (
        smld_fitted = smld_fitted,
        smld_ground_truth = smld_ground_truth,
        images = images,
        detections = roi_batch,
        workflow = workflow
    )
end

"""
    standard_localization_workflow(images, camera; kwargs...) → (smld, workflow)

Standard localization workflow for experimental data: detection → fitting.

Similar to simulate_detect_fit_workflow but starts with existing images
rather than simulation.

# Arguments
- `images`: Camera image stack (2D or 3D array)
- `camera`: Camera object (IdealCamera or SCMOSCamera)

## Keyword Arguments
- `detect_params`: NamedTuple of detection parameters
- `fit_params`: NamedTuple of fitting parameters
- `workflow`: Optional SMLMWorkflow object

# Returns
Named tuple with:
- `smld`: BasicSMLD with fitted localizations
- `detections`: Raw detection results
- `fit_results`: Raw GaussMLE results
- `workflow`: Updated workflow object

# Example
```julia
# Load experimental data
images = load("experiment.tif")
camera = IdealCamera(1:512, 1:512, 0.065f0)

# Run workflow
result = standard_localization_workflow(images, camera;
    detect_params = (minval=50.0, use_gpu=true),
    fit_params = (iterations=30,)
)

# Save results
save_smld(result.smld, "localizations.h5")
```
"""
function standard_localization_workflow(
    images::AbstractArray,
    camera::SMLMData.AbstractCamera;
    detect_params::NamedTuple = NamedTuple(),
    fit_params::NamedTuple = NamedTuple(),
    workflow::Union{SMLMWorkflow,Nothing} = nothing
)
    # Create workflow tracker if not provided
    if workflow === nothing
        workflow = SMLMWorkflow("Detect → Fit")
    end

    # Step 1: Detect particles
    default_detect = (
        boxsize = 11,
        overlap = 2.0,
        sigma_small = 1.0,
        sigma_large = 2.0,
        minval = 10.0,
        use_gpu = true
    )
    detect_params_full = merge(default_detect, detect_params)

    detections = getboxes(images, camera;
                         boxsize = detect_params_full.boxsize,
                         overlap = detect_params_full.overlap,
                         sigma_small = detect_params_full.sigma_small,
                         sigma_large = detect_params_full.sigma_large,
                         minval = detect_params_full.minval,
                         use_gpu = detect_params_full.use_gpu)

    # getboxes now returns ROIBatch directly
    roi_batch = detections

    add_step!(workflow, "Detect Particles", :getboxes,
             Dict{Symbol,Any}(pairs(detect_params_full)...),
             :SMLMBoxer, summarize_boxer_result(roi_batch))

    # Step 3: Fit Gaussian PSFs
    default_fit = (
        psf_model = GaussianXYNB(1.3f0),  # Default PSF sigma
        iterations = 20,
        device = nothing,
        batch_size = 10_000
    )
    fit_params_full = merge(default_fit, fit_params)

    fitter = GaussMLEFitter(;
        psf_model = fit_params_full.psf_model,
        iterations = fit_params_full.iterations,
        device = fit_params_full.device,
        batch_size = fit_params_full.batch_size
    )

    # fit() now returns BasicSMLD directly
    smld = fit(fitter, roi_batch)

    add_step!(workflow, "Fit Gaussian PSFs", :fit,
             Dict{Symbol,Any}(:psf_model => string(typeof(fit_params_full.psf_model)),
                  :iterations => fit_params_full.iterations,
                  :device => string(fitter.device)),
             :GaussMLE, summarize_smld(smld))

    return (
        smld = smld,
        detections = roi_batch,
        workflow = workflow
    )
end
