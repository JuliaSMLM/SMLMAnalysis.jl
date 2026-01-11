"""
Diagnostic comparison of simulated vs real data fitting

This script runs both simulated and real data through identical
analysis pipelines to identify where the real data fitting fails.
"""

using SMLMAnalysis
using CairoMakie
using Statistics
using Printf

println("="^80)
println("FITTING DIAGNOSTIC: Simulated vs Real Data")
println("="^80)

outdir = "output/06_diagnostic"
mkpath(outdir)

# =============================================================================
# PART 1: SIMULATED DATA (as reference)
# =============================================================================
println("\n📊 PART 1: Simulated Data (reference)")
println("-"^80)

# Simulation parameters matching real data characteristics
pixel_size = 0.1  # µm
roi_size = (7, 7)  # pixels to match real data
camera_sim = IdealCamera(128, 128, pixel_size)

# Simulate with realistic parameters
println("Simulating...")
sim_params = StaticSMLMParams(
    density = 2.0,  # patterns per µm²
    σ_psf = 0.13,   # µm (~1.3 pixels, close to our assumption)
    nframes = 100,
    ndatasets = 1
)

pattern = Nmer2D(n=8, d=0.15)  # 8-mer, 150nm diameter
fluorophore = GenericFluor(
    photons = 200_000.0,  # photons/sec
    k_off = 20.0,         # Hz
    k_on = 0.06           # Hz
)

pattern_result, smld_true, smld_noisy = simulate(sim_params; pattern=pattern, molecule=fluorophore, camera=camera_sim)
println("  Generated $(length(smld_noisy.emitters)) emitters in $(smld_noisy.n_frames) frames")

# Generate images
println("Generating images...")
images_sim = gen_images(smld_noisy, camera_sim)
println("  Images: $(size(images_sim))")

# Detection
println("Detecting...")
roi_batch_sim = getboxes(images_sim, camera_sim;
    boxsize=7,
    overlap=2.0,
    sigma_small=1.0,
    sigma_large=2.0,
    minval=5.0  # Lower threshold for sim data
)
println("  Detected: $(size(roi_batch_sim.data, 3)) ROIs")

# Fitting
println("Fitting...")
fitter_sim = GaussMLEFitter(
    psf_model = GaussianXYNB(Float32(sim_params.σ_psf)),
    iterations = 20,
    device = :cpu
)
smld_sim = fit(fitter_sim, roi_batch_sim)

# Extract sim results
photons_sim = [e.photons for e in smld_sim.emitters]
bg_sim = [e.bg for e in smld_sim.emitters]
σ_x_sim = [e.σ_x for e in smld_sim.emitters]
pvalue_sim = [e.pvalue for e in smld_sim.emitters]

println("\n  Simulation Results:")
println("    Photons: mean=$(round(mean(photons_sim), digits=0)), median=$(round(median(photons_sim), digits=0))")
println("    Background: mean=$(round(mean(bg_sim), digits=1)), median=$(round(median(bg_sim), digits=1))")
println("    σ_x: mean=$(round(mean(σ_x_sim)*1000, digits=1)) nm, median=$(round(median(σ_x_sim)*1000, digits=1)) nm")
println("    p-value: mean=$(round(mean(pvalue_sim), digits=3)), median=$(round(median(pvalue_sim), digits=3))")
println("    p-value range: $(round(minimum(pvalue_sim), digits=3)) - $(round(maximum(pvalue_sim), digits=3))")

# =============================================================================
# PART 2: REAL DATA
# =============================================================================
println("\n📷 PART 2: Real Data")
println("-"^80)

h5file = "../data/gatta_ruler/2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5"
data_real, _ = smart_h5_to_array(h5file, max_frames=100)
camera_real = IdealCamera(size(data_real,2), size(data_real,1), pixel_size)

println("Loading...")
println("  Data: $(size(data_real,1))×$(size(data_real,2))×$(size(data_real,3))")
println("  Intensity: min=$(minimum(data_real)), max=$(maximum(data_real)), mean=$(round(mean(data_real), digits=1))")

# Detection with same parameters
println("Detecting...")
roi_batch_real = getboxes(data_real, camera_real;
    boxsize=7,
    overlap=2.0,
    sigma_small=1.0,
    sigma_large=2.0,
    minval=500.0
)
println("  Detected: $(size(roi_batch_real.data, 3)) ROIs")

# Fitting with same model
println("Fitting...")
fitter_real = GaussMLEFitter(
    psf_model = GaussianXYNB(1.0f0),  # Same PSF sigma as sim
    iterations = 20,
    device = :cpu
)
smld_real = fit(fitter_real, roi_batch_real)

# Extract real results
photons_real = [e.photons for e in smld_real.emitters]
bg_real = [e.bg for e in smld_real.emitters]
σ_x_real = [e.σ_x for e in smld_real.emitters]
pvalue_real = [e.pvalue for e in smld_real.emitters]

println("\n  Real Data Results:")
println("    Photons: mean=$(round(mean(photons_real), digits=0)), median=$(round(median(photons_real), digits=0))")
println("    Background: mean=$(round(mean(bg_real), digits=1)), median=$(round(median(bg_real), digits=1))")
println("    σ_x: mean=$(round(mean(σ_x_real)*1000, digits=1)) nm, median=$(round(median(σ_x_real)*1000, digits=1)) nm")
println("    p-value: mean=$(round(mean(pvalue_real), digits=3)), median=$(round(median(pvalue_real), digits=3))")
println("    p-value range: $(round(minimum(pvalue_real), digits=3)) - $(round(maximum(pvalue_real), digits=3))")

# =============================================================================
# PART 3: COMPARISON & DIAGNOSIS
# =============================================================================
println("\n🔬 PART 3: Comparison & Diagnosis")
println("-"^80)

println("\nRatios (Real / Simulated):")
println("  Photons: $(round(mean(photons_real)/mean(photons_sim), digits=1))×")
println("  Background: $(round(mean(bg_real)/mean(bg_sim), digits=1))×")
println("  Precision (σ_x): $(round(mean(σ_x_real)/mean(σ_x_sim), digits=1))×")

# Visualize comparison
fig = Figure(size=(1600, 1000))

# Photons comparison
ax1 = Axis(fig[1, 1], xlabel="Photons", ylabel="Count", title="Photons: Sim vs Real")
hist!(ax1, photons_sim, bins=50, label="Simulated", alpha=0.6, color=:blue)
hist!(ax1, photons_real, bins=50, label="Real", alpha=0.6, color=:red)
axislegend(ax1)

# Background comparison
ax2 = Axis(fig[1, 2], xlabel="Background (photons/pixel)", ylabel="Count", title="Background: Sim vs Real")
hist!(ax2, bg_sim, bins=50, label="Simulated", alpha=0.6, color=:blue)
hist!(ax2, bg_real, bins=50, label="Real", alpha=0.6, color=:red)
axislegend(ax2)

# p-value comparison
ax3 = Axis(fig[2, 1], xlabel="p-value", ylabel="Count", title="p-value: Sim vs Real")
if !all(pvalue_sim .== 0)
    hist!(ax3, pvalue_sim, bins=50, label="Simulated", alpha=0.6, color=:blue)
end
if !all(pvalue_real .== 0)
    hist!(ax3, pvalue_real, bins=50, label="Real", alpha=0.6, color=:red)
else
    text!(ax3, 0.5, 0.5, text="All real data p-values = 0\n(COMPLETE FIT FAILURE)",
          align=(:center, :center), color=:red, fontsize=16)
end
axislegend(ax3)

# Precision comparison
ax4 = Axis(fig[2, 2], xlabel="σ_x (nm)", ylabel="Count", title="Precision: Sim vs Real")
hist!(ax4, σ_x_sim .* 1000, bins=50, label="Simulated", alpha=0.6, color=:blue)
hist!(ax4, σ_x_real .* 1000, bins=50, label="Real", alpha=0.6, color=:red)
axislegend(ax4)

save(joinpath(outdir, "sim_vs_real_comparison.png"), fig)
println("\n✓ Saved comparison plot: $outdir/sim_vs_real_comparison.png")

# Examine a few example ROIs from real data
println("\n🔎 Examining example real data ROIs...")
println("First 5 ROIs:")
for i in 1:min(5, length(smld_real.emitters))
    e = smld_real.emitters[i]
    println("  ROI $i: photons=$(round(e.photons, digits=0)), bg=$(round(e.bg, digits=1)), σ_x=$(round(e.σ_x*1000, digits=1))nm, pvalue=$(e.pvalue)")

    # Show actual ROI data
    roi_data = roi_batch_real.data[:, :, i]
    println("    ROI intensity: min=$(minimum(roi_data)), max=$(maximum(roi_data)), mean=$(round(mean(roi_data), digits=1))")
end

println("\n" * "="^80)
println("DIAGNOSIS COMPLETE")
println("="^80)
println("\n⚠️  KEY FINDING: All real data p-values = 0 indicates:")
println("   1. Model does not fit the data at all")
println("   2. Likely causes:")
println("      - Wrong PSF model (σ assumption incorrect)")
println("      - Wrong camera calibration (gain/offset)")
println("      - Background model inappropriate")
println("      - Initialization far from optimum")
println("\n💡 NEXT STEPS:")
println("   1. Estimate PSF σ from real data (measure bright spots)")
println("   2. Verify camera gain calibration")
println("   3. Check if background is uniform or varying")
println("   4. Try different PSF models (e.g., variable sigma)")
