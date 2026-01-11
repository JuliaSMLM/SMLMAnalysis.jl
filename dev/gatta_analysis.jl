# Gatta Ruler Analysis - Isotropic Sigma
# Uses GaussianXYNBS (fits single σ parameter)
# See gatta_analysis_sxsy.jl for anisotropic version

using SMLMAnalysis

# =============================================================================
# Data
# =============================================================================
h5file = "../data/gatta_ruler/2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5"

println("Loading data...")
data, info = smart_h5_to_array(h5file)  # All 20k frames
# data, info = smart_h5_to_array(h5file; max_frames=1000)  # Quick test
println("  Loaded $(size(data, 3)) frames ($(size(data, 1))×$(size(data, 2)))")

# =============================================================================
# Camera Setup - ORCA-Fusion (C14440-20UP)
# =============================================================================
# Calibrated parameters for this microscope:
#   - Pixel size: 78nm effective (6.5μm sensor / ~83x magnification)
#   - Gain: 0.24 e-/ADU (SMITE convention: 1/4.16 ADU/e-)
#   - Offset: 100 ADU (dark level)
#   - Readnoise: 0.7 e- RMS (ultra quiet mode)
#   - QE: 80% at 642nm

camera = SCMOSCamera(size(data, 2), size(data, 1), 0.078f0, 0.7f0;
    offset = 100.0f0, gain = 0.24f0, qe = 0.8f0)

# =============================================================================
# Analysis
# =============================================================================
result = analyze(data, camera;
    # Detection
    minval = 500.0,

    # Fitting - uses GaussianXYNBS (variable sigma) by default
    # fit_model = :variable,  # default

    # Filtering
    min_photons = 500.0,
    min_pvalue = 1e-3,
    max_sigma = 0.015,  # 15nm precision threshold

    # Frame connection
    frameconnect = false,  # Skip for now

    # Rendering
    render = true,
    render_zoom = 20,

    # Output
    outdir = "output/gatta/"
)

# =============================================================================
# Summary
# =============================================================================
println("\n" * "="^60)
println("Analysis Complete")
println("="^60)
println(result)

# Quick stats
emitters = result.smld.emitters
if length(emitters) > 0
    using Statistics
    σ = [e.σ for e in emitters]
    photons = [e.photons for e in emitters]
    println("\nFiltered localizations:")
    println("  Count: $(length(emitters))")
    println("  Photons: median=$(round(median(photons), digits=0))")
    println("  PSF σ: median=$(round(median(σ)*1000, digits=1)) nm")
end
