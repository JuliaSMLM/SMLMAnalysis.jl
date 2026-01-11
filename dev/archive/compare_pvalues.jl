"""
Quick comparison of p-values: simulated vs real data

The key diagnostic is that p-values should be uniform [0,1] for good fits.
Real data shows all p-values = 0, indicating complete model failure.
"""

using SMLMAnalysis
using Statistics

println("="^80)
println("P-VALUE DIAGNOSTIC")
println("="^80)

# Load and analyze real data (quick - just 100 frames, 1000 ROIs)
println("\n📷 Real Data Analysis...")
h5file = "../data/gatta_ruler/2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5"
data, _ = smart_h5_to_array(h5file, max_frames=100)
camera = IdealCamera(size(data,2), size(data,1), 0.1)

# Get just a subset of ROIs for quick testing
roi_batch = getboxes(data, camera;
    boxsize=7, overlap=2.0, sigma_small=1.0, sigma_large=2.0, minval=500.0
)

nrois_total = size(roi_batch.data, 3)
println("  Total ROIs detected: $nrois_total")

# Fit with different PSF sigma values to test model sensitivity
println("\n🔬 Testing different PSF sigma values...")
println("(Testing first 100 ROIs for speed)")
sigma_values = [0.8f0, 1.0f0, 1.2f0, 1.5f0, 2.0f0]

for psf_sigma in sigma_values
    # Create subset of first 100 ROIs by slicing the data
    roi_subset = ROIBatch(
        roi_batch.data[:, :, 1:100],
        roi_batch.x_corners[1:100],
        roi_batch.y_corners[1:100],
        roi_batch.frame_indices[1:100],
        roi_batch.camera
    )

    fitter = GaussMLEFitter(
        psf_model = GaussianXYNB(psf_sigma),
        iterations = 20,
        device = :cpu
    )

    smld = fit(fitter, roi_subset)
    pvals = [e.pvalue for e in smld.emitters]

    println("\n  PSF σ = $(psf_sigma) pixels:")
    println("    p-value: min=$(round(minimum(pvals), digits=4)), max=$(round(maximum(pvals), digits=4))")
    println("    p-value: mean=$(round(mean(pvals), digits=4)), median=$(round(median(pvals), digits=4))")
    println("    Non-zero p-values: $(sum(pvals .> 0)) / $(length(pvals))")

    # Check some example fits
    println("    Example results (first 3):")
    for i in 1:min(3, length(smld.emitters))
        e = smld.emitters[i]
        roi = roi_subset.data[:, :, i]
        println("      ROI $i: photons=$(round(e.photons, digits=0)), bg=$(round(e.bg, digits=1)), pval=$(e.pvalue)")
        println("             ROI range: $(minimum(roi))-$(maximum(roi)) ADU")
    end
end

println("\n" * "="^80)
println("DIAGNOSIS:")
println("="^80)
println("If p-values are all 0 across different PSF sigma:")
println("  → Model fundamentally doesn't fit the data")
println("  → Possible causes:")
println("    1. Background not uniform (varying across ROI)")
println("    2. PSF shape not Gaussian")
println("    3. Camera calibration completely wrong")
println("    4. Numerical issues in MLE computation")
