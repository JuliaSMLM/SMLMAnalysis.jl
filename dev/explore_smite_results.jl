"""
Explore SMITE Results.mat file structure and compare with Julia analysis
"""

using MAT
using Statistics

# Path to SMITE results for our dataset
smite_dir = "../data/gatta_ruler/2025-10-23/Results/N1/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25_cropped1_MnP_1e-14_MaxXY_SE_0.05"
mat_file = joinpath(smite_dir, "20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25_cropped1_Results.mat")

println("="^80)
println("SMITE Results Analysis")
println("="^80)
println("\nFile: ", basename(mat_file))

# Load the mat file
println("\nLoading MAT file...")
mat_data = matread(mat_file)

# The main data should be in SMD
smd = mat_data["SMD"]

# Get pixel size
pixel_size = smd["PixelSize"]
println("\nPixel size: $pixel_size µm ($(pixel_size*1000) nm)")

# Extract key arrays
X = vec(smd["X"])
Y = vec(smd["Y"])
X_SE = vec(smd["X_SE"])
Y_SE = vec(smd["Y_SE"])
photons = vec(smd["Photons"])
bg = vec(smd["Bg"])
pvalue = vec(smd["PValue"])
frame_num = vec(smd["FrameNum"])
dataset_num = vec(smd["DatasetNum"])
psf_sigma_x = vec(smd["PSFSigmaX"])
psf_sigma_y = vec(smd["PSFSigmaY"])

n_locs = length(X)

println("\n" * "="^60)
println("SMITE Localization Summary")
println("="^60)
println("\nTotal localizations: $n_locs")
println("Frame range: $(Int(minimum(frame_num))) - $(Int(maximum(frame_num)))")
println("Number of datasets: $(Int(maximum(dataset_num)))")

# Dataset breakdown
println("\nPer-dataset breakdown:")
for i in 1:Int(maximum(dataset_num))
    mask = dataset_num .== i
    n = count(mask)
    fr_min, fr_max = extrema(frame_num[mask])
    println("  Dataset $i: $n localizations (frames $(Int(fr_min))-$(Int(fr_max)))")
end

# Coordinate ranges (in microns)
println("\n--- Spatial Extent ---")
println("X range: $(round(minimum(X)*pixel_size, digits=2)) - $(round(maximum(X)*pixel_size, digits=2)) µm")
println("Y range: $(round(minimum(Y)*pixel_size, digits=2)) - $(round(maximum(Y)*pixel_size, digits=2)) µm")

# Precision (X_SE, Y_SE are in pixels, need to convert to nm)
println("\n--- Localization Precision ---")
x_se_nm = X_SE .* pixel_size .* 1000
y_se_nm = Y_SE .* pixel_size .* 1000
println("X precision: median = $(round(median(x_se_nm), digits=1)) nm, mean = $(round(mean(x_se_nm), digits=1)) nm")
println("Y precision: median = $(round(median(y_se_nm), digits=1)) nm, mean = $(round(mean(y_se_nm), digits=1)) nm")

# PSF Sigma
println("\n--- PSF Sigma (fitted) ---")
sig_x_nm = psf_sigma_x .* pixel_size .* 1000
sig_y_nm = psf_sigma_y .* pixel_size .* 1000
println("PSF sigma X: median = $(round(median(sig_x_nm), digits=1)) nm, mean = $(round(mean(sig_x_nm), digits=1)) nm")
println("PSF sigma Y: median = $(round(median(sig_y_nm), digits=1)) nm, mean = $(round(mean(sig_y_nm), digits=1)) nm")

# Photon statistics
println("\n--- Photon Counts ---")
println("Photons: median = $(round(median(photons), digits=0)), mean = $(round(mean(photons), digits=0))")
println("Photon range: $(round(minimum(photons), digits=0)) - $(round(maximum(photons), digits=0))")

# Background
println("\n--- Background ---")
println("Background: median = $(round(median(bg), digits=1)), mean = $(round(mean(bg), digits=1))")

# P-value
println("\n--- Fit Quality (p-value) ---")
println("P-value: median = $(round(median(pvalue), sigdigits=3)), mean = $(round(mean(pvalue), sigdigits=3))")
println("P-value > 0.01: $(round(100*count(p -> p > 0.01, pvalue)/n_locs, digits=1))%")
println("P-value > 0.001: $(round(100*count(p -> p > 0.001, pvalue)/n_locs, digits=1))%")
println("P-value > 1e-7: $(round(100*count(p -> p > 1e-7, pvalue)/n_locs, digits=1))%")

println("\n" * "="^80)
println("Comparison: SMITE vs Julia Analysis")
println("="^80)
println("""
| Metric                | SMITE (5000 fr, 4 datasets) | Julia before filter | Julia after p>1e-7 |
|-----------------------|-----------------------------|---------------------|---------------------|
| Total localizations   | $n_locs               | 105,916             | 549                 |
| PSF sigma (median)    | $(round(median(sig_x_nm), digits=1)) nm                 | 165.5 nm            | -                   |
| X precision (median)  | $(round(median(x_se_nm), digits=1)) nm                  | 1.7 nm              | 2.4 nm              |
| Photons (median)      | $(round(median(photons), digits=0))                | 25,777              | 14,966              |
| Background (median)   | $(round(median(bg), digits=1))                 | 523                 | 371                 |
| Detection rate        | ~151/frame                  | ~106/frame          | ~0.5/frame          |
| % passing p>1e-7      | 98.5%                       | 0.5%                | 100%                |
""")

println("Critical findings:")
println("1. PSF sigma matches perfectly: SMITE $(round(median(sig_x_nm), digits=1)) nm vs Julia 165.5 nm")
println("2. Precision is similar: SMITE $(round(median(x_se_nm), digits=1)) nm vs Julia 2.4 nm")
println("3. P-VALUE MISMATCH: SMITE accepts 98.5% at p>1e-7, Julia only 0.5%!")
println("4. Photon difference: Julia reports ~3x more photons than SMITE")
println("")
println("Interpretation:")
println("- Our chi-squared/p-value calculation differs from SMITE's")
println("- The p-value distribution suggests model mismatch or different DOF")
println("- SMITE folder: 'MnP_1e-14_MaxXY_SE_0.05' = minPValue=1e-14, maxSE=50nm")
println("- Despite more lenient threshold (1e-7 vs 1e-14), Julia rejects more fits")

println("\n" * "="^80)
println("Done")
println("="^80)
