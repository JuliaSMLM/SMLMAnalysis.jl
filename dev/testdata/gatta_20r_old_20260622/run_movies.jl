# DETECTFIT-ONLY at DEBUG -> per-gallery-frame realtime MP4s for the old-20R data.
# Produces <file>/01_detectfit/frame_movies/detection_frame_<N>.mp4 matching the real
# (full 19-block) gallery frames. fps = realtime = 1/exposure (40fps @25ms, 10fps @100ms).
#
# Run alone (single 4090; concurrent GPU jobs OOM):
#   CUDA_VISIBLE_DEVICES=0 julia -t auto --project=<repo> \
#     dev/testdata/gatta_20r_old_20260622/run_movies.jl   (optional: append a path substring)

using SMLMAnalysis

const PIXEL_SIZE = 0.078f0
const READNOISE  = 1.4f0
const OFFSET     = 100.0f0
const GAIN       = 1.0f0/2.2f0
const QE         = 1.0f0

const BASE = @__DIR__
const REPO = abspath(joinpath(BASE, "..", "..", ".."))
const RAW  = "/mnt/nas/lidkelab/Projects/HS-TIRF/Data/Gattaquant_ruler_PAINT_HiRes_20R/6-22-2026"
const OUT  = joinpath(REPO, "dev", "output", "gatta_20r_old_20260622_movies")
const MARK = joinpath(OUT, "_markers"); mkpath(MARK)

isdir(RAW) || error("raw data dir missing: $RAW")
h5files = String[]
for (root,_,files) in walkdir(RAW), f in files
    endswith(f, ".h5") && push!(h5files, joinpath(root, f))
end
sort!(h5files)
filters = ARGS

println("="^70)
println("OLD-20R detectfit DEBUG frame-movies  ($(length(h5files)) files)")
println("raw -> $RAW"); println("out -> $OUT"); println("="^70); flush(stdout)

rows = NamedTuple[]
for h5 in h5files
    tag = first(splitext(relpath(h5, RAW)))
    (!isempty(filters) && !any(f -> occursin(f, tag), filters)) && continue
    safe = replace(tag, r"[^A-Za-z0-9_.-]" => "_")
    isfile(joinpath(MARK, safe*".done")) && (println("  SKIP $tag (done)"); flush(stdout); continue)
    outdir = joinpath(OUT, safe); rm(outdir; force=true, recursive=true)
    info = load_mic_h5_info(h5)
    ds = 1:max(1, info.n_blocks - 1)                      # drop corrupt final block
    m = match(r"([0-9.]+)s", tag)                         # exposure from condition dir (e.g. "0.025s")
    fps = m === nothing ? 20.0 : 1.0/parse(Float64, m.captures[1])
    print("  $tag  ($(info.n_frames)f/$(info.n_blocks)blk -> $(length(ds)), $(round(info.file_size_gb,digits=2))GB, $(round(fps,digits=1))fps) ... "); flush(stdout)
    try
        camera = SCMOSCamera(info.width, info.height, PIXEL_SIZE, READNOISE; offset=OFFSET, gain=GAIN, qe=QE)
        cfg = AnalysisConfig(camera=camera, outdir=outdir, verbose=Verbosity.DEBUG, steps=[
            DetectFitConfig(path=h5, h5_format=:mic, datasets=ds, movie_fps=fps,
                boxer=BoxerConfig(boxsize=9, psf_sigma=0.130, min_photons=200.0),
                fitter=GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20)),
        ])
        analyze(cfg)
        mdir = joinpath(outdir, "01_detectfit", "frame_movies")
        nmov = isdir(mdir) ? count(f->endswith(f,".mp4"), readdir(mdir)) : 0
        write(joinpath(MARK, safe*".done"), relpath(mdir, REPO))
        # surface this file's first movie in the FE as it comes in (badge; explicit slug — detached run has no tmux)
        try
            mp4s = sort(filter(f -> endswith(f, ".mp4"), readdir(mdir; join=true)))
            isempty(mp4s) || run(`$(homedir())/.sot-comm/bin/devenv-fe preview smlmanalysis $(relpath(mp4s[1], REPO))`)
        catch e
            @warn "FE preview failed for $tag" exception=e
        end
        println("$nmov movies @ $(round(fps,digits=1))fps"); flush(stdout)
        push!(rows, (tag=tag, n=nmov, status="ok"))
    catch err
        write(joinpath(MARK, safe*".fail"), sprint(showerror, err))
        println("FAILED"); flush(stdout); @error "movie pass $tag failed" exception=(err, catch_backtrace())
        push!(rows, (tag=tag, n=missing, status="error"))
    end
end

println("\n" * "="^70); println("SUMMARY"); println("="^70)
for r in rows; println(rpad(r.tag,30), lpad(r.n===missing ? "-" : string(r.n),5), " movies  ", r.status); end
ok = count(r->r.status=="ok", rows)
println("\n$ok/$(length(rows)) ok -> $OUT")
write(joinpath(MARK,"_ALL_DONE"), string(ok,"/",length(rows))); flush(stdout)
