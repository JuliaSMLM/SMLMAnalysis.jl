# Gatta Ruler Analysis - Anisotropic Sigma
# Uses GaussianXYNBSXSY (fits σx, σy separately)
# Useful for detecting PSF asymmetry, astigmatism, or z-dependent widths

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
camera = SCMOSCamera(size(data, 2), size(data, 1), 0.078f0, 0.7f0;
    offset = 100.0f0, gain = 0.24f0, qe = 0.8f0)

# =============================================================================
# Analysis - Anisotropic PSF model
# =============================================================================
result = analyze(data, camera;
    # Detection - PSF sigma in microns (DoG uses σ and 2σ)
    psf_sigma = 0.135f0,  # ~135nm typical for TIRF
    detect_min_photons = 500.0,  # Detection threshold in photons

    # Fitting - anisotropic sigma (σx, σy)
    fit_model = :anisotropic,

    # Filtering
    min_photons = 500.0,
    max_precision = 0.015,        # 15nm precision threshold
    psf_sigma_mode_tolerance = 0.10,  # Keep PSF sigma within ±10% of mode
    min_pvalue = 1e-3,

    # Frame connection
    frameconnect = false,

    # Rendering
    render = true,
    render_zoom = 20,

    # Output
    outdir = joinpath(@__DIR__, "output", "gatta_sxsy")
)

# =============================================================================
# Summary
# =============================================================================
println("\n" * "="^60)
println("Analysis Complete (Anisotropic σx/σy)")
println("="^60)
println(result)

# Anisotropic stats - check for PSF asymmetry
emitters = result.smld.emitters
if length(emitters) > 0
    using Statistics

    # GaussianXYNBSXSY emitters have σx, σy fields
    σx = [e.σx for e in emitters]
    σy = [e.σy for e in emitters]
    photons = [e.photons for e in emitters]

    println("\nFiltered localizations:")
    println("  Count: $(length(emitters))")
    println("  Photons: median=$(round(median(photons), digits=0))")
    println("  PSF σx: median=$(round(median(σx)*1000, digits=1)) nm")
    println("  PSF σy: median=$(round(median(σy)*1000, digits=1)) nm")
    println("  Asymmetry (σx/σy): $(round(median(σx)/median(σy), digits=2))")
end
