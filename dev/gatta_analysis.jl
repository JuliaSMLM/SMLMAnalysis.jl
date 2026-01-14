# Gatta Ruler Analysis - Anisotropic Sigma
# Uses GaussianXYNBSXSY (fits separate σx, σy parameters)
# This handles elliptical PSFs and gives better pvalues

using SMLMAnalysis

# =============================================================================
# Data
# =============================================================================
h5file = joinpath(@__DIR__, "..", "data", "gatta_ruler", "2025-10-23", "20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5")

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
    offset = 100.0f0, gain = 0.24f0, qe = 1.0f0)

# =============================================================================
# Analysis
# =============================================================================
result = analyze(data, camera;
    # Detection - PSF sigma in microns (DoG uses σ and 2σ)
    psf_sigma = 0.135f0,  # ~135nm typical for TIRF
    detect_min_photons = 1000.0,  # Higher threshold to reduce false detections

    # Fitting - anisotropic PSF (separate σx, σy)
    fit_model = :anisotropic,

    # Filtering
    min_photons = 500.0,
    max_precision = 0.015,        # 15nm precision threshold
    psf_sigma_mode_tolerance = 0.10,  # Keep PSF σx,σy within ±10% of their modes
    min_pvalue = 1e-6,

    # Drift correction - sequential Legendre order 3
    drift_model = "LegendrePoly",
    drift_degree = 3,
    drift_cost_fun = "Kdtree",

    # Isolated emitter filter - remove noise
    filter_isolated = true,
    isolated_n_sigma = 2.0,    # Neighbor if dist < 2σ_combined
    isolated_min_neighbors = :auto, # Triangle method for automatic threshold

    # Frame connection - combine repeated localizations
    frameconnect = true,
    fc_maxframegap = 5,
    fc_nsigmadev = 5.0,

    # Rendering (uses defaults: gaussian@20x, histogram@10x, circles@50x)
    render = true,

    # Output
    outdir = joinpath(@__DIR__, "output", "gatta")
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
    σx = [e.σ_x for e in emitters]
    σy = [e.σ_y for e in emitters]
    photons = [e.photons for e in emitters]
    println("\nFiltered localizations:")
    println("  Count: $(length(emitters))")
    println("  Photons: median=$(round(median(photons), digits=0))")
    println("  PSF σx: median=$(round(median(σx)*1000, digits=1)) nm")
    println("  PSF σy: median=$(round(median(σy)*1000, digits=1)) nm")
end
