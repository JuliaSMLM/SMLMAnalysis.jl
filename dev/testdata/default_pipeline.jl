# Default SMLMAnalysis pipeline for dev/testdata full-cell runs.
#
# CAMERA follows genmab (same scope, Keith-confirmed: 647nm/1.35NA silicon oil, 97.8nm px,
# qe=0.82) -- set in each dataset's run.jl. DETECTION is NOT inherited from genmab: genmab's
# min_photons=200 + tight quality filter FLOODS this data's noise floor (~494k candidate fits
# -> 29 locs on PFA/Cell_03) -- the thresholds are sample/SNR-specific. These are the prior
# working TNFaR1 values; still provisional, tune per dataset. Drift is :continuous (the 4 MIC
# blocks are one acquisition), vs genmab's :registered.
#
# Run at Verbosity.STANDARD so every stage emits its figures (incl. BaGoL's report).
# Render set (Keith's spec): HistogramRender 10x turbo/time, GaussianRender 20x, CircleRender 50x,
# then BaGoL, then GaussianRender 50x of the MAP-N emitters.
#
# Returns the AnalysisConfig.steps vector for one cell.

using SMLMAnalysis

function default_steps(h5;
        # detection + fitting (provisional TNFaR1 values; NOT genmab's -- see header)
        psf_sigma::Float64 = 0.130, min_photons::Float64 = 500.0, boxsize::Int = 9,
        fit_iterations::Int = 20,
        # quality filter
        max_precision::Float64 = 0.015, min_pvalue::Float64 = 1e-6,
        # frame connection
        max_frame_gap::Int = 5, max_sigma_dist::Float64 = 5.0,
        # drift -- :continuous for this acquisition (NOT genmab's :registered)
        dataset_mode::Symbol = :continuous, drift_degree::Int = 2)
    return [
        DetectFitConfig(path = h5, h5_format = :mic,
            boxer  = BoxerConfig(boxsize = boxsize, min_photons = min_photons, psf_sigma = psf_sigma),
            fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), iterations = fit_iterations)),
        FilterConfig(photons = (min_photons, Inf), precision = (0.0, max_precision),
            pvalue = (min_pvalue, 1.0)),
        FrameConnectConfig(max_frame_gap = max_frame_gap, max_sigma_dist = max_sigma_dist,
            calibration = CalibrationConfig(clamp_k_to_one = true)),
        DriftConfig(degree = drift_degree, dataset_mode = dataset_mode,
            quality = :iterative, auto_roi = false),

        # --- three pre-BaGoL renders ---
        RenderConfig(strategy = HistogramRender(), zoom = 10, colormap = :turbo,
            color_by = :absolute_frame, clip_percentile = nothing),       # hist / time / turbo 10x
        RenderConfig(strategy = GaussianRender(),  zoom = 20, colormap = :inferno),  # gaussian 20x
        RenderConfig(strategy = CircleRender(),    zoom = 50, colormap = :turbo,
            color_by = :absolute_frame),                                   # circles 50x

        # --- BaGoL grouping (genmab BAGOL_CFG: learn μ/shape, defaults otherwise) ---
        BaGoLConfig(max_partition_size = 500, learn_distribution = true,
            posterior_pixel_size = 0.002, overlap = 0.025),

        # --- post-BaGoL render: gaussian MAP-N at 50x ---
        RenderConfig(strategy = GaussianRender(), zoom = 50, colormap = :inferno),   # gaussian MAP-N 50x
    ]
end
