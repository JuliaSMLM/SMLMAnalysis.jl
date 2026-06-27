# Sebastian Kg1a TNFaR1 (1:250, ImStrF1) dSTORM -- acquired 2026-06-25.
# FULL-DATASET run using papers-geometry's EXACT SR/dSTORM config.
#
# Config provenance: ~/julia_shared_dev/papers/papers-geometry/src/config.jl
#   (SR_* consts + SR_STEPS_BAGOL = SR_STEPS_PREBAGOL + SR_BAGOL_CONFIG + SR_STEPS_POSTBAGOL),
#   verified against papers-geometry @ main 60c5f6a (via @papers-geometry-kitt, 2026-06-25).
#   Replicated VERBATIM below so this harness does not depend on the geometry repo being
#   loadable from the SMLMAnalysis project. Keith: "same config geometry is using (exactly)".
#
# Camera: per-pixel sCMOS model identical to geometry (pixel=0.0978um, qe=0.82, 0.130um PSF),
#   but the per-pixel offset/gain/readnoise come from EACH cell's OWN embedded Calibration/
#   group (CCDOffset/CCDVar/Gain) via build_camera_from_mic_h5 -- same camera as geometry,
#   freshly calibrated with this acquisition (geometry only uses a fixed cal file because its
#   older Rachel data lacks an embedded Calibration group; this data carries its own).
#
# DEVIATION NOTE: drift dataset_mode=:registered is geometry's value (used here for an exact
#   match). The prior 20260610 TNFaR1 harness used :continuous because the MIC blocks are one
#   movie. Flagged for Keith; switch the const below to :continuous if drift looks wrong.
#
# Output -> dev/output/ (under the workspace root, so the 5nm render is nav-previewable;
#   shared filesystem => descent writes, kitt previews). Per cell: <condition>/Cell_NN/<step>/.
#
# Run (descent, all cells):
#   CUDA_VISIBLE_DEVICES=0 julia -t auto \
#     --project=/home/kalidke/julia_shared_dev/SMLMAnalysis \
#     dev/testdata/sebastian_kg1a_tnfar1_20260625/run.jl
#   (optional cell filter:  ... run.jl 1nM_F1_30ms/Cell_01  500pM)

using SMLMAnalysis
using Statistics

# --- geometry SR camera constants (config.jl:74-76) ---
const SR_PIXEL_SIZE = 0.0978f0   # um
const SR_QE         = 0.82f0
const SR_PSF_SIGMA  = 0.130      # um -- 647nm, 1.35NA silicone (target-independent: AF647)

# --- geometry SR detection/filter/connect/drift consts (config.jl:80-93) ---
const SR_BOXSIZE        = 9
const SR_MIN_PHOTONS    = 200.0
const SR_FIT_ITERATIONS = 20
const SR_MAX_PRECISION  = 0.010   # um (10 nm)
const SR_MIN_PVALUE     = 1e-6
const SR_MAX_FRAME_GAP  = 5
const SR_NSIGMADEV      = 5.0
const SR_DRIFT_DEGREE      = 3
const SR_DRIFT_QUALITY     = :iterative
const SR_DRIFT_AUTO_ROI    = false
const SR_DRIFT_SHIFT_SCALE  = 0.050
const SR_DRIFT_DATASET_MODE = :registered   # geometry's value (see DEVIATION NOTE)

# Geometry's full single-machine SR pipeline (config.jl SR_STEPS_BAGOL = 01..09), with the
# per-cell H5 path injected into the DetectFit step. Steps:
#   01 detectfit  02 filter  03 frameconnect  04 driftcorrect
#   05 render (gaussian 20x inferno = 4.89 nm/px ~= the "5nm render")
#   06 render (histogram 10x turbo, time-colored)  07 render (circle 50x turbo, time-colored)
#   08 bagol      09 render (gaussian 50x inferno, post-BaGoL MAP-N)
function sr_steps(h5)
    return [
        DetectFitConfig(
            path = h5, h5_format = :mic,
            boxer  = BoxerConfig(boxsize = SR_BOXSIZE, psf_sigma = SR_PSF_SIGMA,
                                 min_photons = SR_MIN_PHOTONS),
            fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), iterations = SR_FIT_ITERATIONS)),
        FilterConfig(
            photons   = (SR_MIN_PHOTONS, Inf),
            precision = (0.0, SR_MAX_PRECISION),
            pvalue    = (SR_MIN_PVALUE, 1.0),
            psf_sigma = (0.100, 0.150)),
        FrameConnectConfig(
            max_frame_gap  = SR_MAX_FRAME_GAP,
            max_sigma_dist = SR_NSIGMADEV,
            calibration    = CalibrationConfig(clamp_k_to_one = true)),
        DriftConfig(
            degree       = SR_DRIFT_DEGREE,
            dataset_mode = SR_DRIFT_DATASET_MODE,
            quality      = SR_DRIFT_QUALITY,
            auto_roi     = SR_DRIFT_AUTO_ROI,
            shift_scale  = SR_DRIFT_SHIFT_SCALE),
        # --- three pre-BaGoL renders ---
        RenderConfig(zoom = 20, colormap = :inferno, scalebar = true),                 # 05 = 5nm gaussian
        RenderConfig(strategy = HistogramRender(), zoom = 10, colormap = :turbo,
                     color_by = :absolute_frame, clip_percentile = nothing, scalebar = true),  # 06
        RenderConfig(strategy = CircleRender(), zoom = 50, colormap = :turbo,
                     color_by = :absolute_frame, scalebar = true),                      # 07
        # --- BaGoL (config.jl SR_BAGOL_CONFIG) ---
        BaGoLConfig(
            se_adjust            = :auto,
            n_iterations         = 6000,
            burn_in              = 2000,
            sync_interval        = 100,
            μ                    = 10.0,
            shape                = 1.0,
            learn_distribution   = true,
            max_partition_size   = 500,
            overlap              = 0.025,
            posterior_pixel_size = 0.0),                                                # 08
        # --- post-BaGoL MAP-N render ---
        RenderConfig(zoom = 50, colormap = :inferno, scalebar = true),                 # 09
    ]
end

const BASE = @__DIR__
const REPO = abspath(joinpath(BASE, "..", "..", ".."))      # SMLMAnalysis repo root
const RAW  = "/mnt/nas/gillette/Sebastian/20260625_Kg1a_TNFaR1-1-250_ImStrF1"
const OUT  = joinpath(REPO, "dev", "output", "sebastian_kg1a_tnfar1_20260625")
const MARK = joinpath(OUT, "_markers")
mkpath(MARK)

isdir(RAW) || error("raw data dir missing: $RAW")

# --- discover every cell (Cell_NN/Label_01/Data_*.h5), skipping the Results/ tree ---
h5files = String[]
for (root, _, files) in walkdir(RAW)
    (occursin("/Results", root) || occursin("ipynb_checkpoints", root)) && continue
    endswith(root, "Label_01") || continue
    for f in files
        (startswith(f, "Data_") && endswith(f, ".h5")) && push!(h5files, joinpath(root, f))
    end
end
sort!(h5files)

filters = ARGS  # optional: only cells whose "<condition>/Cell_NN" tag contains an arg

println("="^70)
println("SEBASTIAN Kg1a TNFaR1 (1:250, ImStrF1) -- geometry SR config  ($(length(h5files)) cells)",
        isempty(filters) ? "" : "  filter=$(filters)")
println("raw     -> ", RAW)
println("results -> ", OUT)
println("config  -> papers-geometry SR_STEPS_BAGOL @ main 60c5f6a (verbatim)")
println("="^70); flush(stdout)

safe(tag) = replace(tag, "/" => "__")

rows = NamedTuple[]
for h5 in h5files
    # tag = "<condition>/Cell_NN"  from .../<condition>/Cell_NN/Label_01/Data_*.h5
    parts = splitpath(relpath(h5, RAW))
    cond, cell = parts[1], parts[2]
    tag = "$cond/$cell"
    (!isempty(filters) && !any(f -> occursin(f, tag), filters)) && continue

    donemark = joinpath(MARK, safe(tag) * ".done")
    if isfile(donemark)
        println("  SKIP $tag (already done)"); flush(stdout)
        push!(rows, (cond=cond, cell=cell, n=missing, status="skip"))
        continue
    end

    outdir = joinpath(OUT, cond, cell)
    rm(outdir; force=true, recursive=true)
    info = load_mic_h5_info(h5)
    print("  $tag  ($(info.n_frames)f / $(info.n_blocks) blocks, ",
          "$(round(info.file_size_gb, digits=2))GB) ... "); flush(stdout)
    try
        camera = build_camera_from_mic_h5(h5; pixel_size = SR_PIXEL_SIZE, qe = SR_QE)
        cfg = AnalysisConfig(camera = camera, steps = sr_steps(h5),
                             outdir = outdir, verbose = Verbosity.STANDARD)
        (result, _) = analyze(cfg)
        n = length(result.smld.emitters)

        # locate the 5nm render (gaussian zoom=20) for the nav preview; store repo-relative path
        png5 = ""
        for d in filter(d -> occursin("_render", d), readdir(outdir; join=true))
            for p in filter(p -> endswith(p, ".png") && occursin("gaussianrender", p) &&
                                 occursin("20x", p), readdir(d; join=true))
                png5 = relpath(p, REPO)
            end
        end
        write(donemark, png5)
        println("$n emitters | 5nm: $png5"); flush(stdout)
        push!(rows, (cond=cond, cell=cell, n=n, status="ok"))
    catch err
        write(joinpath(MARK, safe(tag) * ".fail"), sprint(showerror, err))
        println("FAILED"); flush(stdout)
        @error "cell $tag failed" exception=(err, catch_backtrace())
        push!(rows, (cond=cond, cell=cell, n=missing, status="error"))
    end
end

println("\n" * "="^70); println("SUMMARY"); println("="^70)
for r in rows
    println(rpad(r.cond, 18), rpad(r.cell, 10),
            lpad(r.n === missing ? "-" : string(r.n), 9), "  ", r.status)
end
ok = count(r -> r.status == "ok", rows)
println("\n$ok/$(length(rows)) cells ok -> ", OUT)
write(joinpath(MARK, "_ALL_DONE"), string(ok, "/", length(rows)))
flush(stdout)
