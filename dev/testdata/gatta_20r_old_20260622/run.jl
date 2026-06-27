# GATTAquant PAINT HiRes 20R (OLD ruler, 2D) on HS-TIRF Mic -- newest acquisition 2026-06-22.
# Same pipeline as the gatta_3dy run: papers-geometry SR config (2D) + Orca-Flash4.0 camera.
#
# Sample: GATTA-PAINT HiRes 20R (2D 20nm nanoruler). MIC H5, 500x500, 20 blocks
#   (last block Data0020 corrupt in every HS-TIRF file -> dropped, 1:19). No embedded
#   Calibration group -> scalar Orca-Flash4.0 camera (Keith: gain 1/2.2, RN 1.4).
#   Two exposures: 0.025s (25ms) and 0.1s (100ms), 4 files each.
#
# Config: papers-geometry SR_STEPS_BAGOL @ main 60c5f6a (verbatim, 2D GaussianXYNBS).
# Output -> dev/output/ (nav-previewable). Run (kitt, local GPU):
#   CUDA_VISIBLE_DEVICES=0 julia -t auto \
#     --project=/home/kalidke/julia_shared_dev/SMLMAnalysis \
#     dev/testdata/gatta_20r_old_20260622/run.jl     (optional: append a path substring)

using SMLMAnalysis

# --- camera: Orca-Flash4.0 (Keith 2026-06-26: gain 1/2.2, RN 1.4) ---
const PIXEL_SIZE = 0.078f0
const READNOISE  = 1.4f0
const OFFSET     = 100.0f0
const GAIN       = 1.0f0/2.2f0
const QE         = 1.0f0

# --- geometry SR consts (config.jl), verbatim ---
const SR_PSF_SIGMA = 0.130
const SR_BOXSIZE, SR_MIN_PHOTONS, SR_FIT_ITERATIONS = 9, 200.0, 20
const SR_MAX_PRECISION, SR_MIN_PVALUE = 0.010, 1e-6
const SR_MAX_FRAME_GAP, SR_NSIGMADEV = 5, 5.0
const SR_DRIFT_DEGREE, SR_DRIFT_QUALITY, SR_DRIFT_AUTO_ROI, SR_DRIFT_SHIFT_SCALE = 3, :iterative, false, 0.050
const SR_DRIFT_DATASET_MODE = :registered

function sr_steps(h5, datasets)
    return [
        DetectFitConfig(path = h5, h5_format = :mic, datasets = datasets,
            boxer  = BoxerConfig(boxsize = SR_BOXSIZE, psf_sigma = SR_PSF_SIGMA, min_photons = SR_MIN_PHOTONS),
            fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), iterations = SR_FIT_ITERATIONS)),
        FilterConfig(photons = (SR_MIN_PHOTONS, Inf), precision = (0.0, SR_MAX_PRECISION),
            pvalue = (SR_MIN_PVALUE, 1.0), psf_sigma = (0.100, 0.150)),
        FrameConnectConfig(max_frame_gap = SR_MAX_FRAME_GAP, max_sigma_dist = SR_NSIGMADEV,
            calibration = CalibrationConfig(clamp_k_to_one = true)),
        DriftConfig(degree = SR_DRIFT_DEGREE, dataset_mode = SR_DRIFT_DATASET_MODE,
            quality = SR_DRIFT_QUALITY, auto_roi = SR_DRIFT_AUTO_ROI, shift_scale = SR_DRIFT_SHIFT_SCALE),
        RenderConfig(zoom = 20, colormap = :inferno, scalebar = true),                  # 05 = 5nm gaussian
        RenderConfig(strategy = HistogramRender(), zoom = 10, colormap = :turbo,
            color_by = :absolute_frame, clip_percentile = nothing, scalebar = true),     # 06
        RenderConfig(strategy = CircleRender(), zoom = 50, colormap = :turbo,
            color_by = :absolute_frame, scalebar = true),                                # 07
        BaGoLConfig(se_adjust = :auto, n_iterations = 6000, burn_in = 2000, sync_interval = 100,
            μ = 10.0, shape = 1.0, learn_distribution = true,
            max_partition_size = 500, overlap = 0.025, posterior_pixel_size = 0.0),      # 08
        RenderConfig(zoom = 50, colormap = :inferno, scalebar = true),                   # 09 post-BaGoL
    ]
end

const BASE = @__DIR__
const REPO = abspath(joinpath(BASE, "..", "..", ".."))
const RAW  = "/mnt/nas/lidkelab/Projects/HS-TIRF/Data/Gattaquant_ruler_PAINT_HiRes_20R/6-22-2026"
const OUT  = joinpath(REPO, "dev", "output", "gatta_20r_old_20260622")
const MARK = joinpath(OUT, "_markers"); mkpath(MARK)

isdir(RAW) || error("raw data dir missing: $RAW")

# discover H5 under exposure subdirs (0.025s/, 0.1s/)
h5files = String[]
for (root,_,files) in walkdir(RAW), f in files
    endswith(f, ".h5") && push!(h5files, joinpath(root, f))
end
sort!(h5files)
filters = ARGS

println("="^70)
println("GATTA PAINT HiRes 20R (OLD, HS-TIRF, 6-22-2026) -- geometry SR config 2D  ($(length(h5files)) files)",
        isempty(filters) ? "" : "  filter=$(filters)")
println("raw     -> ", RAW)
println("results -> ", OUT)
println("camera  -> Orca-Flash4.0: 78nm, gain $(round(GAIN,digits=4)) (1/2.2), offset $OFFSET, RN $READNOISE, qe $QE")
println("="^70); flush(stdout)

rows = NamedTuple[]
for h5 in h5files
    tag = first(splitext(relpath(h5, RAW)))          # e.g. "0.1s/1-2026-..."
    (!isempty(filters) && !any(f -> occursin(f, tag), filters)) && continue
    safe = replace(tag, r"[^A-Za-z0-9_.-]" => "_")
    donemark = joinpath(MARK, safe * ".done")
    if isfile(donemark)
        println("  SKIP $tag (done)"); flush(stdout); push!(rows,(tag=tag,n=missing,status="skip")); continue
    end
    outdir = joinpath(OUT, safe); rm(outdir; force=true, recursive=true)
    info = load_mic_h5_info(h5)
    ds = 1:max(1, info.n_blocks - 1)                  # drop corrupt final block
    print("  $tag  ($(info.n_frames)f / $(info.n_blocks) blocks -> $(length(ds)), $(info.width)x$(info.height), ",
          "$(round(info.file_size_gb, digits=2))GB) ... "); flush(stdout)
    try
        camera = SCMOSCamera(info.width, info.height, PIXEL_SIZE, READNOISE; offset=OFFSET, gain=GAIN, qe=QE)
        cfg = AnalysisConfig(camera=camera, steps=sr_steps(h5, ds), outdir=outdir, verbose=Verbosity.STANDARD)
        (result,_) = analyze(cfg)
        n = length(result.smld.emitters)
        png5 = ""
        for d in filter(d->occursin("_render",d), readdir(outdir;join=true)),
            p in filter(p->endswith(p,".png") && occursin("gaussianrender",p) && occursin("20x",p), readdir(d;join=true))
            png5 = relpath(p, REPO)
        end
        write(donemark, png5)
        println("$n emitters | 5nm: $png5"); flush(stdout)
        push!(rows,(tag=tag,n=n,status="ok"))
    catch err
        write(joinpath(MARK, safe*".fail"), sprint(showerror,err))
        println("FAILED"); flush(stdout)
        @error "file $tag failed" exception=(err, catch_backtrace())
        push!(rows,(tag=tag,n=missing,status="error"))
    end
end

println("\n" * "="^70); println("SUMMARY"); println("="^70)
for r in rows
    println(rpad(r.tag,30), lpad(r.n===missing ? "-" : string(r.n),9), "  ", r.status)
end
ok = count(r->r.status=="ok", rows)
println("\n$ok/$(length(rows)) ok -> ", OUT)
write(joinpath(MARK,"_ALL_DONE"), string(ok,"/",length(rows))); flush(stdout)
