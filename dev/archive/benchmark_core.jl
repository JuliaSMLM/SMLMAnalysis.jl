"""
Minimal benchmark script for SMLM workflow performance analysis.

No plotting - just times core operations with allocation tracking.
Run twice: first for precompilation, second for actual timing.
"""

using SMLMAnalysis
using SMLMDriftCorrection
using Statistics

# Config
h5file = "../data/gatta_ruler/2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5"
max_frames = 1000

println("="^60)
println("SMLM Workflow Benchmark")
println("="^60)

# Load data (outside timing loop - just once)
println("\nLoading data...")
stats = @timed begin
    data, info = smart_h5_to_array(h5file, max_frames=max_frames)
end
data, info = stats.value
println("  Load: $(round(stats.time, digits=2))s, $(round(stats.bytes/1e6, digits=1))MB allocated")
println("  Data size: $(size(data))")

# Setup camera
camera = SCMOSCamera(size(data,2), size(data,1), 0.1f0, 0.7f0;
    offset=100.0f0, gain=0.24f0, qe=0.80f0)

# Fitter setup - use GPU
fitter = GaussMLEFitter(psf_model=GaussianXYNBS(), iterations=20)

println("\n" * "-"^60)
println("Benchmarking core operations (2 runs)")
println("-"^60)

for run in 1:2
    println("\n--- Run $run $(run==1 ? "(includes precompile)" : "(warmed up)") ---")

    # Detection (GPU)
    stats = @timed roi_batch = getboxes(data, camera;
        boxsize=9, overlap=2.0, sigma_small=1.0, sigma_large=2.0,
        minval=500.0, use_gpu=true)
    n_rois = length(roi_batch)
    t_detect = stats.time
    println("  Detection: $(round(stats.time, digits=2))s, $(round(stats.bytes/1e6, digits=1))MB, $n_rois ROIs")

    # Fitting
    stats = @timed smld = fit(fitter, roi_batch)
    n_fits = length(smld.emitters)
    t_fit = stats.time
    println("  Fitting:   $(round(stats.time, digits=2))s, $(round(stats.bytes/1e6, digits=1))MB, $n_fits fits")

    # Filter (simple p-value filter)
    good_emitters = filter(e -> e.pvalue > 1e-7, smld.emitters)
    smld_filt = BasicSMLD(good_emitters, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)
    n_filt = length(smld_filt.emitters)
    println("  Filtered:  $n_filt localizations ($(round(100*n_filt/n_fits, digits=1))%)")

    # Drift correction
    stats = @timed smld_dc = driftcorrect(smld_filt)
    println("  Drift:     $(round(stats.time, digits=2))s, $(round(stats.bytes/1e6, digits=1))MB")

    # Total
    println("  ---")
    println("  Throughput: $(round(n_rois/t_detect, digits=0)) ROIs/s (detection)")
    println("  Throughput: $(round(n_fits/t_fit, digits=0)) fits/s (fitting)")
end

println("\n" * "="^60)
println("Benchmark complete")
println("="^60)
