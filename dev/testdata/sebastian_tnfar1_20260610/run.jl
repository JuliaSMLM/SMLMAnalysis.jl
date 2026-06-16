# Sebastian TNFaR1 dSTORM -- FULL-DATASET run, papers-style data/results layout.
#
#   data/raw/tnfar1        -> symlink to the NAS acquisition dir
#   data/results/juliasmlm -> symlink to a NAS results dir; per-cell output goes BACK
#                             to the NAS as <condition>/Cell_NN/<step>/  (see README.md)
#
# Per cell runs the shared DEFAULT pipeline (../default_pipeline.jl): all per-stage
# figures + 3 renders (hist/time 10x, gaussian 20x, circle 50x) + BaGoL + gaussian MAP-N 50x.
#
# Run (all cells):      julia -t auto --project=<repo> .../run.jl
#     (subset by tag):  julia -t auto --project=<repo> .../run.jl PFA/Cell_03 MeOH

using SMLMAnalysis
using Statistics

include(joinpath(@__DIR__, "..", "default_pipeline.jl"))   # default_steps(h5; ...)

const BASE    = @__DIR__
const RAW     = joinpath(BASE, "data", "raw", "tnfar1")          # -> NAS acquisition dir
const RESULTS = joinpath(BASE, "data", "results", "juliasmlm")   # -> NAS results (per-cell)
const MAX_GB  = 3.0          # skip > this (PFA_250nM/Cell_02=3.5GB, Cell_03=16.6GB); raise to include
const PIXEL_SIZE = 0.0978f0  # 97.8 nm -- genmab camera (Keith: TNFaR1 uses the same scope), qe=0.82

isdir(RAW)     || error("data/raw/tnfar1 symlink missing -- see README.md")
isdir(RESULTS) || error("data/results/juliasmlm symlink missing -- see README.md")

function build_config(h5, info, outdir)
    cal = load_mic_h5_calibration_for_scmos(h5)
    camera = SCMOSCamera(info.width, info.height, PIXEL_SIZE, cal.readnoise;
        offset = cal.offset, gain = cal.gain, qe = 0.82f0)
    # STANDARD so every stage emits its diagnostic figures -- incl. BaGoL's report/
    # circles/partitions, which bagol_step gates behind verbose >= STANDARD.
    AnalysisConfig(camera = camera, steps = default_steps(h5),
                   outdir = outdir, verbose = Verbosity.STANDARD)
end

# --- discover every cell (Data_*.h5) under data/raw ---
h5files = String[]
for (root, _, files) in walkdir(RAW)
    occursin("juliasmlm", root) && continue          # skip the co-located results tree
    for f in files
        (startswith(f, "Data_") && endswith(f, ".h5")) && push!(h5files, joinpath(root, f))
    end
end
sort!(h5files)

filters = ARGS   # optional: only run cells whose "<condition>/Cell_NN" tag contains an arg

println("="^64)
println("SEBASTIAN TNFaR1 -- FULL DATASET  ($(length(h5files)) cells)",
        isempty(filters) ? "" : "  filter=$(filters)")
println("results -> ", readlink(RESULTS))
println("="^64)

rows = Vector{NamedTuple}()
for h5 in h5files
    parts = splitpath(relpath(h5, RAW))
    cond, cell = parts[1], parts[2]
    tag = "$cond/$cell"
    (!isempty(filters) && !any(f -> occursin(f, tag), filters)) && continue
    info = load_mic_h5_info(h5)
    if info.file_size_gb > MAX_GB
        println("  SKIP $tag  ($(round(info.file_size_gb, digits=1)) GB > MAX_GB=$MAX_GB)")
        push!(rows, (cond=cond, cell=cell, nloc=missing, status="skipped(size)"))
        continue
    end
    outdir = joinpath(RESULTS, cond, cell)
    rm(outdir; force=true, recursive=true)
    print("  $tag  ($(info.n_frames)f, $(round(info.file_size_gb, digits=2))GB) ... ")
    try
        (result, _) = analyze(build_config(h5, info, outdir))
        n = length(result.smld.emitters)
        println("$n emitters")
        push!(rows, (cond=cond, cell=cell, nloc=n, status="ok"))
    catch err
        println("FAILED")
        @error "cell $tag failed" exception=(err, catch_backtrace())
        push!(rows, (cond=cond, cell=cell, nloc=missing, status="error"))
    end
end

println("\n" * "="^64); println("SUMMARY"); println("="^64)
println(rpad("condition", 22), rpad("cell", 10), lpad("emitters", 9), "  status")
for r in rows
    println(rpad(r.cond, 22), rpad(r.cell, 10),
            lpad(r.nloc === missing ? "-" : string(r.nloc), 9), "  ", r.status)
end
ok = count(r -> r.status == "ok", rows)
println("\n$ok/$(count(r -> true, rows)) cells -> ", readlink(RESULTS))
println("NOTE: pixel size + detection params still placeholders -- counts low until tuned.")
