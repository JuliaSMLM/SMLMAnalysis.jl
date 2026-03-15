# Quick test: detectfit -> filter -> intensityfilter only
# Stops after intensity filter to iterate on diagnostics fast

using SMLMAnalysis

h5file = joinpath(@__DIR__, "..", "data", "gatta_ruler", "2025-10-23",
    "20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5")

info = load_smart_h5_info(h5file)
println("Loading $(info.nframes) frames...")
t_load = @elapsed images, _ = smart_h5_to_array(h5file)
println("  Loaded in $(round(t_load, digits=1))s")

camera = SCMOSCamera(info.width, info.height, 0.078f0, 0.7f0;
    offset = 100.0f0, gain = 0.24f0, qe = 1.0f0)

roi = (x = info.width÷2+1:info.width, y = 1:info.height)
outdir = joinpath(@__DIR__, "output", "test_intensityfilter")
rm(outdir; force=true, recursive=true)

config = AnalysisConfig(
    camera = camera,
    steps = [
        DetectFitConfig(
            boxer = BoxerConfig(boxsize=11, min_photons=1000.0, psf_sigma=0.135),
            fitter = GaussMLEConfig(psf_model=GaussianXYNBS(), iterations=20),
        ),
        FilterConfig(
            photons = (500.0, Inf),
            precision = (0.0, 0.007),
            pvalue = (1e-6, 1.0),
            psf_sigma = :auto
        ),
        IntensityFilterConfig(),
    ],
    roi = roi,
    outdir = outdir,
    verbose = Verbosity.DETAILED,
)

(result, analysis_info) = analyze([images], config)
println("\nDone. Output: $outdir")
println("Localizations: $(length(result.smld.emitters))")
