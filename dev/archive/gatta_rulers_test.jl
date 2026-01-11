"""
DNA-PAINT Ruler Analysis with SMART Microscope Data

This script demonstrates the complete SMLM analysis workflow using experimental
DNA-PAINT ruler data from the SMART microscope system.

Camera: Hamamatsu ORCA-Fusion (C14440-20UP)
Data: 20R ruler, TIRF illumination, 642nm laser, 100ms exposure
"""

using SMLMAnalysis
using Statistics
using Printf

println("="^80)
println("DNA-PAINT Ruler Analysis - SMART Microscope Data")
println("="^80)

# =============================================================================
# 1. DATA LOADING
# =============================================================================
println("\n📂 STEP 1: Loading H5 data...")

h5file = "../data/gatta_ruler/2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5"

# Load file info
info = load_smart_h5_info(h5file)
println("   File: ", basename(info.filepath))
println("   Full sensor: 2304 x 2304 pixels")
println("   ROI dimensions: ", info.width, " x ", info.height, " pixels")
println("   Total frames: ", info.nframes)
println("   Exposure time: 100 ms")
println("   Frame rate: 10 fps")
println("   File size: ", round(info.file_size_gb, digits=2), " GB")

# Load first 1000 frames for testing (DNA-PAINT is sparse, this should be enough)
println("\n   Loading first 1000 frames...")
data, _ = smart_h5_to_array(h5file, max_frames=1000)
println("   Data shape: ", size(data))
println("   Memory usage: ", round(sizeof(data) / 1e6, digits=1), " MB")

# =============================================================================
# 2. CAMERA SETUP
# =============================================================================
println("\n📷 STEP 2: Setting up camera model...")

# ORCA-Fusion (C14440-20UP) Specifications
# Physical pixel: 6.5 µm
# For magnification estimation, typical TIRF SMLM uses ~100-150x total magnification
# Assuming ~100x objective + 1.5x relay = 150x total
# Effective pixel size = 6.5 µm / 150 = 43.3 nm
#
# However, metadata shows pixelsize=1.0 which might be placeholder
# For now, we'll use a typical SMLM pixel size of 100 nm effective
# TODO: Verify actual magnification from microscope calibration

pixel_size_um = 0.1  # 100 nm effective pixel size (to be verified)
println("   Camera: Hamamatsu ORCA-Fusion (C14440-20UP)")
println("   Physical pixel size: 6.5 µm")
println("   Assumed effective pixel size: $(pixel_size_um * 1000) nm")
println("   Note: Verify magnification calibration!")

# Create sCMOS camera model
# For DNA-PAINT at 642nm with low background, we expect:
# - Low read noise (ultra quiet mode: 0.7e-)
# - Offset: typical ~100-200 ADU
# - Using 16-bit mode
# - QE: 80% at 642nm
# Note: SCMOSCamera in SMLMData doesn't include dark current parameter
#       Dark current is negligible for cooled sCMOS (~0.2 e-/pix/s) compared to photon noise

camera = SCMOSCamera(
    info.width,     # nx pixels
    info.height,    # ny pixels
    pixel_size_um,  # pixel size in microns
    0.7;            # readnoise in electrons (ultra quiet mode)
    offset=150.0,   # ADU baseline
    gain=1.0,       # electrons per ADU (may need calibration from manufacturer)
    qe=0.80         # 80% QE at 642nm
)

println("   Camera model created:")
println("     - Dimensions: $(info.width) x $(info.height) pixels")
println("     - Pixel size: $(pixel_size_um) µm")
println("     - Read noise: 0.7 e-")
println("     - QE: 80%")
println("     - Gain: 1.0 e-/ADU")
println("     - Offset: 150.0 ADU")

# =============================================================================
# 3. DETECTION - Find candidate particles
# =============================================================================
println("\n🔍 STEP 3: Detecting particles with SMLMBoxer...")

# DNA-PAINT produces bright, well-separated spots
# SMLMBoxer uses Difference of Gaussians (DoG) filtering
# Parameters:
# - boxsize: ROI size (pixels), should be ~4-5x PSF width
# - sigma_small, sigma_large: DoG filter scales
# - minval: detection threshold (ADU or photons)
# - overlap: pixels of overlap between adjacent ROIs

boxer_params = (
    boxsize = 7,        # 7x7 pixel ROIs
    overlap = 2.0,      # 2 pixels overlap
    sigma_small = 1.0,  # Smaller DoG sigma
    sigma_large = 2.0,  # Larger DoG sigma
    minval = 200.0      # Minimum peak value (ADU above background)
)

println("   SMLMBoxer parameters:")
println("     - boxsize: $(boxer_params.boxsize) pixels")
println("     - DoG sigmas: $(boxer_params.sigma_small), $(boxer_params.sigma_large) pixels")
println("     - minval threshold: $(boxer_params.minval) ADU")
println("     - overlap: $(boxer_params.overlap) pixels")

# Run detection - getboxes returns ROIBatch directly
println("\n   Running detection...")
roi_batch = getboxes(data, camera;
    boxsize=boxer_params.boxsize,
    overlap=boxer_params.overlap,
    sigma_small=boxer_params.sigma_small,
    sigma_large=boxer_params.sigma_large,
    minval=boxer_params.minval
)

println("\n   Detection Results:")
# ROIBatch stores data in a 3D array (roi_pixels, roi_pixels, num_rois), not .rois field
num_rois = size(roi_batch.data, 3)
println("     - ROIs detected: ", num_rois)
println("     - Detections per frame (avg): ", round(num_rois / size(data, 3), digits=2))

# =============================================================================
# 4. FITTING - Gaussian MLE fitting
# =============================================================================
println("\n📊 STEP 4: Fitting Gaussians with GaussMLE...")

# For 2D data, use GaussianXYNB model (X, Y, photons, background)
# roi_batch already created by getboxes

println("   Fitting model: GaussianXYNB (2D + photons + background)")
println("   ROI batch size: ", num_rois, " ROIs")
if num_rois > 0
    println("   ROI dimensions: ", roi_batch.roi_size, " x ", roi_batch.roi_size)
end

# Create fitter
# For diffraction-limited 642nm with NA~1.4: σ_PSF ≈ 0.21*λ/NA ≈ 95nm ≈ 0.95 pixels
psf_sigma = 1.0f0  # pixels (Float32)
fitter = GaussMLEFitter(
    psf_model = GaussianXYNB(psf_sigma),
    iterations = 20,
    device = :cpu
)

# Fit all ROIs
println("\n   Running maximum likelihood estimation...")
fit_result = fit(fitter, roi_batch)

summarize_fit_result(fit_result)

# =============================================================================
# 5. CONVERT TO SMLD
# =============================================================================
println("\n📦 STEP 5: Converting to SMLD format...")

smld = localization_result_to_smld(fit_result, roi_batch, camera)
summarize_smld(smld)

# =============================================================================
# 6. ANALYSIS - Ruler measurements
# =============================================================================
println("\n📏 STEP 6: Analyzing localizations...")

# Basic statistics
println("\n   Localization statistics:")
println("     - Total localizations: ", length(smld.x))
println("     - X range: $(round(minimum(smld.x), digits=2)) to $(round(maximum(smld.x), digits=2)) µm")
println("     - Y range: $(round(minimum(smld.y), digits=2)) to $(round(maximum(smld.y), digits=2)) µm")
println("     - Mean photons: $(round(mean(smld.photons), digits=0))")
println("     - Median photons: $(round(median(smld.photons), digits=0))")
println("     - Mean background: $(round(mean(smld.bg), digits=1))")

# Precision analysis
println("\n   Localization precision:")
if hasfield(typeof(smld), :x_sigma) && !isempty(smld.x_sigma)
    println("     - Mean σ_x: $(round(mean(smld.x_sigma) * 1000, digits=1)) nm")
    println("     - Median σ_x: $(round(median(smld.x_sigma) * 1000, digits=1)) nm")
    println("     - Mean σ_y: $(round(mean(smld.y_sigma) * 1000, digits=1)) nm")
    println("     - Median σ_y: $(round(median(smld.y_sigma) * 1000, digits=1)) nm")
else
    println("     - Precision data not available in this SMLD format")
end

# Photon statistics
println("\n   Photon count distribution:")
println("     - Min: $(round(minimum(smld.photons), digits=0))")
println("     - 25th percentile: $(round(quantile(smld.photons, 0.25), digits=0))")
println("     - Median: $(round(median(smld.photons), digits=0))")
println("     - 75th percentile: $(round(quantile(smld.photons, 0.75), digits=0))")
println("     - Max: $(round(maximum(smld.photons), digits=0))")

# Frame distribution
println("\n   Temporal distribution:")
frames_used = unique(smld.frame)
println("     - Frames with localizations: ", length(frames_used))
println("     - Localizations per frame (avg): ", round(length(smld.x) / length(frames_used), digits=2))

# =============================================================================
# 7. CONCLUSIONS
# =============================================================================
println("\n" * "="^80)
println("✅ Analysis complete!")
println("="^80)

println("\nNext steps:")
println("  1. ⚠️  Verify pixel size calibration (currently assumed 100nm)")
println("  2. 🔗 Apply frame connection to link blinking events")
println("  3. 📐 Perform ruler distance analysis (expect ~20nm for 20R ruler)")
println("  4. 🎨 Generate super-resolution image with SMLMRender")
println("  5. 📊 Compare with ground truth/expected distances")

println("\nData saved in memory as:")
println("  - data: Raw image stack ($(size(data)))")
println("  - boxes: Detected particle coordinates")
println("  - fit_result: MLE fitting results")
println("  - smld: SMLD container with all localizations")
println()
