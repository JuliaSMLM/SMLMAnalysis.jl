"""
DNA-PAINT Ruler Analysis - Simplified Version with IdealCamera

Uses IdealCamera to avoid SCMOSCamera compatibility issues.
This version works end-to-end with the experimental data.
"""

using SMLMAnalysis
using Statistics
using Printf
using CairoMakie

println("="^80)
println("DNA-PAINT Ruler Analysis - Simplified (IdealCamera)")
println("="^80)

# Output directory
outdir = "output"
mkpath(outdir)
println("Output directory: dev/$outdir/")

# Load data
println("\n📂 Loading data...")
h5file = "../data/gatta_ruler/2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5"
info = load_smart_h5_info(h5file)
println("   ROI: $(info.width) x $(info.height) pixels, $(info.nframes) frames")

# Load first 1000 frames
data, _ = smart_h5_to_array(h5file, max_frames=1000)
println("   Loaded: $(size(data,3)) frames ($(round(sizeof(data)/1e6, digits=1)) MB)")

# Camera setup - SCMOSCamera with Float32 parameters to avoid type mismatch bug
# data is (rows, cols, frames) = (256, 860, 1000)
# Hamamatsu ORCA-Fusion (C14440-20UP) specs:
#   - Physical pixel: 6.5 µm
#   - Effective pixel size: ~0.1 µm (with magnification - needs calibration)
#   - Readnoise: 0.7 e⁻ RMS (ultra quiet mode)
#   - QE: 80% at 642nm
#   - Gain: ~0.24 e⁻/ADU (16-bit mode: 15000 e⁻ full well / 65536 ADU)
#   - Offset: ~100 ADU (typical baseline)
# Note: SMITE uses inverted gain convention (ADU/e⁻ = 4.16), our convention is e⁻/ADU
println("\n📷 Camera: SCMOSCamera (ORCA-Fusion specs, Float32)")
pixel_size_um = 0.1f0  # Effective pixel size (needs magnification calibration)
readnoise = 0.7f0      # e⁻ RMS (ultra quiet mode)
offset_adu = 100.0f0   # ADU
gain = 0.24f0          # e⁻/ADU (1/4.16 to match SMITE's convention)
qe = 0.80f0            # QE at 642nm

camera = SCMOSCamera(size(data,2), size(data,1), pixel_size_um, readnoise;
    offset=offset_adu, gain=gain, qe=qe)
println("   Dimensions: $(size(data,2)) x $(size(data,1)) pixels")
println("   Pixel size: $pixel_size_um µm")
println("   Readnoise: $readnoise e⁻ RMS")
println("   Offset: $offset_adu ADU, Gain: $gain e⁻/ADU, QE: $qe")

# Detection with adjusted threshold to reduce false positives
println("\n🔍 Detection...")
println("   Parameters: boxsize=9, sigma_small=1.0, sigma_large=2.0, minval=500")
roi_batch = getboxes(data, camera;
    boxsize=9,  # Must be odd
    overlap=2.0,
    sigma_small=1.0,
    sigma_large=2.0,
    minval=500.0  # Increased from 200 to reduce false positives
)
num_rois = size(roi_batch.data, 3)
println("   Detected: $num_rois ROIs ($(round(num_rois/size(data,3), digits=1)) per frame)")

# Fitting - use GaussianXYNBS to estimate sigma from data
println("\n📊 Fitting...")
# Use XYNBS model: fits X, Y, N (photons), Background, and Sigma
# This lets us estimate the actual PSF width from the data
fitter = GaussMLEFitter(
    psf_model = GaussianXYNBS(),  # Variable sigma, estimated from data
    iterations = 20,
    device = :cpu
)

t_fit = @elapsed fit_result = fit(fitter, roi_batch)
println("   Fitted $num_rois ROIs in $(round(t_fit, digits=1))s")

# fit_result is already a BasicSMLD with emitters
println("\n📦 Fit results (BasicSMLD format)...")
smld = fit_result
println("   Localizations: $(length(smld.emitters))")

# Extract data from emitters (before filtering)
x_all = [e.x for e in smld.emitters]
y_all = [e.y for e in smld.emitters]
photons_all = [e.photons for e in smld.emitters]
bg_all = [e.bg for e in smld.emitters]
σ_x_all = [e.σ_x for e in smld.emitters]
σ_y_all = [e.σ_y for e in smld.emitters]
pvalue_all = [e.pvalue for e in smld.emitters]

# Extract fitted PSF sigma (from XYNBS model)
sigma_all = [e.σ for e in smld.emitters]

println("\n📏 Analysis (before filtering):")
println("   Total localizations: $(length(smld.emitters))")
println("   X range: $(round(minimum(x_all), digits=2)) - $(round(maximum(x_all), digits=2)) µm")
println("   Y range: $(round(minimum(y_all), digits=2)) - $(round(maximum(y_all), digits=2)) µm")
println("   Mean photons: $(round(mean(photons_all), digits=0))")
println("   Median photons: $(round(median(photons_all), digits=0))")
println("   Mean background: $(round(mean(bg_all), digits=1))")
println("\n   Fitted PSF sigma:")
println("     Mean: $(round(mean(sigma_all)*1000, digits=1)) nm")
println("     Median: $(round(median(sigma_all)*1000, digits=1)) nm")
println("     Min: $(round(minimum(sigma_all)*1000, digits=1)) nm")
println("     Max: $(round(maximum(sigma_all)*1000, digits=1)) nm")
println("\n   X uncertainty (σ_x):")
println("     Mean: $(round(mean(σ_x_all)*1000, digits=1)) nm")
println("     Median: $(round(median(σ_x_all)*1000, digits=1)) nm")
println("   Y uncertainty (σ_y):")
println("     Mean: $(round(mean(σ_y_all)*1000, digits=1)) nm")
println("     Median: $(round(median(σ_y_all)*1000, digits=1)) nm")

# P-value statistics
println("\n📊 P-value statistics:")
println("   Min: $(minimum(pvalue_all))")
println("   Max: $(maximum(pvalue_all))")
println("   Mean: $(round(mean(pvalue_all), digits=6))")
println("   Median: $(round(median(pvalue_all), digits=6))")
println("   % > 0.01: $(round(100*count(p -> p > 0.01, pvalue_all)/length(pvalue_all), digits=2))%")
println("   % > 0.001: $(round(100*count(p -> p > 0.001, pvalue_all)/length(pvalue_all), digits=2))%")
println("   % > 1e-5: $(round(100*count(p -> p > 1e-5, pvalue_all)/length(pvalue_all), digits=3))%")
println("   % > 1e-7: $(round(100*count(p -> p > 1e-7, pvalue_all)/length(pvalue_all), digits=3))%")

# Generate diagnostic plots
println("\n📈 Generating diagnostic plots...")

# Helper: compute empirical CDF
function ecdf_values(x)
    sorted = sort(x)
    n = length(sorted)
    return sorted, (1:n) ./ n
end

# Create fit quality diagnostic panel (2x3) with log scales and CDFs
fig = Figure(size=(1800, 1200))

# 1. P-value: log histogram + CDF overlay (CRITICAL diagnostic)
ax1 = Axis(fig[1, 1], xlabel="log₁₀(p-value)", ylabel="Density", title="P-value Distribution")
ax1_cdf = Axis(fig[1, 1], ylabel="CDF", yaxisposition=:right, yticklabelcolor=:red)
hidespines!(ax1_cdf)
hidexdecorations!(ax1_cdf)
pval_nonzero = pvalue_all[pvalue_all .> 0]
if !isempty(pval_nonzero)
    log_pval = log10.(pval_nonzero)
    hist!(ax1, log_pval, bins=50, color=(:steelblue, 0.7), normalization=:pdf)
    # CDF overlay
    cdf_x, cdf_y = ecdf_values(log_pval)
    lines!(ax1_cdf, cdf_x, cdf_y, color=:red, linewidth=2, label="CDF")
    # Threshold lines
    vlines!(ax1, [log10(1e-7)], color=:orange, linestyle=:dash, linewidth=2)
    text!(ax1, log10(1e-7)+0.3, 0.1, text="p=1e-7", fontsize=12)
end

# 2. Photons: log scale histogram + CDF
ax2 = Axis(fig[1, 2], xlabel="Photons", ylabel="Density", title="Photon Distribution", xscale=log10)
ax2_cdf = Axis(fig[1, 2], ylabel="CDF", yaxisposition=:right, yticklabelcolor=:red, xscale=log10)
hidespines!(ax2_cdf)
hidexdecorations!(ax2_cdf)
photons_pos = photons_all[photons_all .> 0]
hist!(ax2, photons_pos, bins=50, color=(:green, 0.7), normalization=:pdf)
cdf_x, cdf_y = ecdf_values(photons_pos)
lines!(ax2_cdf, cdf_x, cdf_y, color=:red, linewidth=2)
vlines!(ax2, [median(photons_pos)], color=:black, linestyle=:dash)

# 3. Fitted PSF sigma histogram + CDF
ax3 = Axis(fig[1, 3], xlabel="PSF σ (nm)", ylabel="Density", title="Fitted PSF Sigma")
ax3_cdf = Axis(fig[1, 3], ylabel="CDF", yaxisposition=:right, yticklabelcolor=:red)
hidespines!(ax3_cdf)
hidexdecorations!(ax3_cdf)
sigma_nm = sigma_all .* 1000
hist!(ax3, sigma_nm, bins=50, color=(:purple, 0.7), normalization=:pdf)
cdf_x, cdf_y = ecdf_values(sigma_nm)
lines!(ax3_cdf, cdf_x, cdf_y, color=:red, linewidth=2)
vlines!(ax3, [median(sigma_nm)], color=:black, linestyle=:dash)

# 4. Background: log scale + CDF
ax4 = Axis(fig[2, 1], xlabel="Background (photons)", ylabel="Density", title="Background Distribution", xscale=log10)
ax4_cdf = Axis(fig[2, 1], ylabel="CDF", yaxisposition=:right, yticklabelcolor=:red, xscale=log10)
hidespines!(ax4_cdf)
hidexdecorations!(ax4_cdf)
bg_pos = bg_all[bg_all .> 0]
hist!(ax4, bg_pos, bins=50, color=(:orange, 0.7), normalization=:pdf)
cdf_x, cdf_y = ecdf_values(bg_pos)
lines!(ax4_cdf, cdf_x, cdf_y, color=:red, linewidth=2)
vlines!(ax4, [median(bg_pos)], color=:black, linestyle=:dash)

# 5. Precision (σ_x): log scale + CDF
ax5 = Axis(fig[2, 2], xlabel="σ_x (nm)", ylabel="Density", title="X Precision (CRLB)", xscale=log10)
ax5_cdf = Axis(fig[2, 2], ylabel="CDF", yaxisposition=:right, yticklabelcolor=:red, xscale=log10)
hidespines!(ax5_cdf)
hidexdecorations!(ax5_cdf)
prec_nm = σ_x_all .* 1000
prec_pos = prec_nm[prec_nm .> 0]
hist!(ax5, prec_pos, bins=50, color=(:teal, 0.7), normalization=:pdf)
cdf_x, cdf_y = ecdf_values(prec_pos)
lines!(ax5_cdf, cdf_x, cdf_y, color=:red, linewidth=2)
vlines!(ax5, [median(prec_pos)], color=:black, linestyle=:dash)

# 6. Photons vs Precision scatter (log-log)
ax6 = Axis(fig[2, 3], xlabel="Photons", ylabel="σ_x (nm)", title="Photons vs Precision",
           xscale=log10, yscale=log10)
scatter!(ax6, photons_all[1:min(5000, end)], σ_x_all[1:min(5000, end)] .* 1000,
         markersize=2, alpha=0.3, color=:blue)

save("$outdir/fit_quality_panel.png", fig)
println("   ✓ Saved fit_quality_panel.png")

# Filter by pvalue (reject bad fits with low p-value)
println("\n🔬 Filtering with pvalue > 1e-7 (reject bad fits)...")
pval_threshold = 1e-7
good_idx = pvalue_all .> pval_threshold
filtered_emitters = smld.emitters[good_idx]
println("   Kept: $(length(filtered_emitters)) / $(length(smld.emitters)) localizations ($(round(100*length(filtered_emitters)/length(smld.emitters), digits=1))%)")

# Handle case where no localizations pass filter
if length(filtered_emitters) == 0
    println("\n⚠️  WARNING: No localizations passed the p-value filter!")
    println("   This likely indicates a model mismatch:")
    println("   - PSF sigma (currently 1.0 pixel) may not match real PSF")
    println("   - Noise model (IdealCamera = Poisson) doesn't match sCMOS")
    println("   - Background estimation may be off")
    println("\n   Proceeding with unfiltered data for visualization...")
    filtered_emitters = smld.emitters
end

# Create filtered SMLD
smld_filtered = BasicSMLD(filtered_emitters, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)

# Extract filtered data
x = [e.x for e in filtered_emitters]
y = [e.y for e in filtered_emitters]
photons = [e.photons for e in filtered_emitters]
bg = [e.bg for e in filtered_emitters]
σ_x = [e.σ_x for e in filtered_emitters]
σ_y = [e.σ_y for e in filtered_emitters]

println("\n📏 After filtering:")
println("   X range: $(round(minimum(x), digits=2)) - $(round(maximum(x), digits=2)) µm")
println("   Y range: $(round(minimum(y), digits=2)) - $(round(maximum(y), digits=2)) µm")
println("   Mean photons: $(round(mean(photons), digits=0))")
println("   Median photons: $(round(median(photons), digits=0))")
println("   Mean background: $(round(mean(bg), digits=1))")
println("\n   X uncertainty (σ_x):")
println("     Mean: $(round(mean(σ_x)*1000, digits=1)) nm")
println("     Median: $(round(median(σ_x)*1000, digits=1)) nm")
println("   Y uncertainty (σ_y):")
println("     Mean: $(round(mean(σ_y)*1000, digits=1)) nm")
println("     Median: $(round(median(σ_y)*1000, digits=1)) nm")

# Render super-resolution images (using filtered data)
println("\n🎨 Rendering images (filtered data)...")

# Gaussian render with inferno colormap
println("   Rendering Gaussian (inferno colormap)...")
t_render = @elapsed result = render(smld_filtered;
    strategy=GaussianRender(),
    zoom=20,  # 20x zoom for high resolution
    colormap=:inferno,
    filename="$outdir/ruler_gaussian_inferno.png"
)
println("   ✓ Saved ruler_gaussian_inferno.png ($(round(t_render, digits=1))s, $(size(result.image)))")

# Histogram render
println("   Rendering histogram (viridis colormap)...")
t_render = @elapsed result = render(smld_filtered;
    strategy=HistogramRender(),
    zoom=20,
    colormap=:viridis,
    filename="$outdir/ruler_histogram_viridis.png"
)
println("   ✓ Saved ruler_histogram_viridis.png ($(round(t_render, digits=1))s, $(size(result.image)))")

# Circle render colored by frame (turbo colormap) - higher zoom for detail
println("   Rendering circles colored by frame (turbo colormap, zoom=40)...")
t_render = @elapsed result = render(smld_filtered;
    strategy=CircleRender(),
    zoom=40,  # Higher zoom for circle visibility (100 causes OOM)
    color_by=:frame,
    colormap=:turbo,
    filename="$outdir/ruler_circles_frame.png"
)
println("   ✓ Saved ruler_circles_frame.png ($(round(t_render, digits=1))s, $(size(result.image)))")

println("\n" * "="^80)
println("✅ Analysis complete!")
println("="^80)
println("\nResults in $outdir/:")
println("  • fit_quality_panel.png - Diagnostic panel (p-value, photons, sigma, etc.)")
println("  • ruler_gaussian_inferno.png - Super-resolution with Gaussian blur")
println("  • ruler_histogram_viridis.png - Histogram visualization")
println("  • ruler_circles_frame.png - Circles colored by frame (zoom=40)")
println("\nAnalysis summary:")
println("  • Localizations: $(length(filtered_emitters)) (after pvalue < $pval_threshold filter)")
println("  • Mean precision: $(round(mean(σ_x)*1000, digits=1)) nm")
println("  • Mean photons: $(round(mean(photons), digits=0))")
println("\nNext steps:")
println("  • Apply frame connection for blinking events")
println("  • Measure ruler distances (expect ~20nm for 20R)")
println("  • Verify/calibrate pixel size")
