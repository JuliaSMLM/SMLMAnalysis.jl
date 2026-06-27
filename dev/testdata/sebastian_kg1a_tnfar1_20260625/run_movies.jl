# DETECTFIT-ONLY at DEBUG -> per-gallery-frame realtime MP4s for the Sebastian Kg1a TNFaR1 data.
# Mirrors run.jl's discovery + camera (per-cell embedded sCMOS calibration), but runs ONLY
# detectfit at Verbosity.DEBUG so the frame-movie feature fires (no SR pipeline / BaGoL).
# Realtime fps = 1000/exposure_ms, parsed from the condition dir name (e.g. 1nM_F1_30ms -> 33.3fps).
#
# Run alone (single GPU): CUDA_VISIBLE_DEVICES=0 julia -t auto --project=<repo> \
#   dev/testdata/sebastian_kg1a_tnfar1_20260625/run_movies.jl   (optional: append a path substring)

using SMLMAnalysis

const PIXEL_SIZE = 0.0978f0      # genmab/TNFaR1 scope
const QE         = 0.82f0
const SR_PSF_SIGMA = 0.130
const BASE = @__DIR__
const REPO = abspath(joinpath(BASE, "..", "..", ".."))
const RAW  = "/mnt/nas/gillette/Sebastian/20260625_Kg1a_TNFaR1-1-250_ImStrF1"
const OUT  = joinpath(REPO, "dev", "output", "sebastian_kg1a_tnfar1_20260625_movies")
const MARK = joinpath(OUT, "_markers"); mkpath(MARK)

isdir(RAW) || error("raw data dir missing: $RAW")
h5files = String[]
for (root,_,files) in walkdir(RAW)
    (occursin("/Results", root) || occursin("ipynb_checkpoints", root)) && continue
    endswith(root, "Label_01") || continue
    for f in files
        (startswith(f, "Data_") && endswith(f, ".h5")) && push!(h5files, joinpath(root, f))
    end
end
sort!(h5files)
filters = ARGS

println("="^70); println("SEBASTIAN Kg1a TNFaR1 detectfit-DEBUG frame-movies  ($(length(h5files)) cells)")
println("raw -> $RAW"); println("out -> $OUT"); println("="^70); flush(stdout)

rows = NamedTuple[]
for h5 in h5files
    parts = splitpath(relpath(h5, RAW)); cond, cell = parts[1], parts[2]
    tag = "$cond/$cell"
    (!isempty(filters) && !any(f -> occursin(f, tag), filters)) && continue
    safe = replace(tag, r"[^A-Za-z0-9_.-]" => "_")
    isfile(joinpath(MARK, safe*".done")) && (println("  SKIP $tag (done)"); flush(stdout); continue)
    outdir = joinpath(OUT, safe); rm(outdir; force=true, recursive=true)
    info = load_mic_h5_info(h5)
    m = match(r"_(\d+)ms", cond)                          # exposure ms from condition name
    fps = m === nothing ? 20.0 : 1000.0/parse(Int, m.captures[1])
    print("  $tag  ($(info.n_frames)f/$(info.n_blocks)blk, $(round(info.file_size_gb,digits=2))GB, $(round(fps,digits=1))fps) ... "); flush(stdout)
    try
        camera = build_camera_from_mic_h5(h5; pixel_size = PIXEL_SIZE, qe = QE)
        cfg = AnalysisConfig(camera=camera, outdir=outdir, verbose=Verbosity.DEBUG, steps=[
            DetectFitConfig(path=h5, h5_format=:mic, movie_fps=fps,
                boxer=BoxerConfig(boxsize=9, psf_sigma=SR_PSF_SIGMA, min_photons=200.0),
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

println("\n"*"="^70); println("SUMMARY"); println("="^70)
for r in rows; println(rpad(r.tag,26), lpad(r.n===missing ? "-" : string(r.n),5), " movies  ", r.status); end
ok = count(r->r.status=="ok", rows)
println("\n$ok/$(length(rows)) ok -> $OUT")
write(joinpath(MARK,"_ALL_DONE"), string(ok,"/",length(rows))); flush(stdout)
