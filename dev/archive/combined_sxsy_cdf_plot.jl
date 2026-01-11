"""
Combined σx/σy CDF Plot

Load data and re-run sxsy fitting for 1ch and 2ch files.
Creates combined CDF plot with all PSF width distributions.

Visual encoding:
- Color: Dataset (1ch=black, pol1=blue, pol2=orange)
- Linestyle: Dimension (σx=solid, σy=dashed)
"""

using SMLMAnalysis
using SMLMBoxer
using GaussMLE
using SMLMData
using CairoMakie
using Statistics
using HDF5

# Configuration - match workflow parameters
const PARAMS = (
    psf_sigma = 0.150,           # PSF sigma in microns
    min_photons = 500.0,
    boxsize = 9,
    pixel_size_um = 0.078f0,
    readnoise = 1.38f0,
    offset_adu = 100.0f0,
    gain = 0.73f0,
    qe = 0.95f0,
    pvalue_threshold = 0.001,
    max_sigma_xy_nm = 10.0,
    psf_tolerance = 0.20,
)

# File paths
base_path = "/mnt/nas/adapt/projects/smart-microscope/data/DNA paint ruler/2025-10-23"
file_1ch = joinpath(base_path, "20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5")
file_2ch = joinpath(base_path, "20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-2ch--2025-10-23_12-04-53.h5")

# Output
outdir = joinpath(base_path, "Results_JuliaSMLM")
outfile = joinpath(outdir, "combined_sxsy_cdf.png")

# Helper: ecdf values
function ecdf_values(data)
    sorted = sort(data)
    n = length(sorted)
    cdf = (1:n) ./ n
    return sorted, cdf
end

# Process a single file and extract PSF widths
function process_file(filepath::String, roi_x_pixels, is_2ch::Bool)
    println("\n" * "="^60)
    println("Processing: $(basename(filepath))")
    println("="^60)

    # Load data using workflow helper
    print("  Loading... ")
    data, _ = smart_h5_to_array(filepath; verbose=false)
    println("$(size(data)) loaded")

    # Apply ROI if specified
    # Data is (height, width, frames) from smart_h5_to_array
    if roi_x_pixels !== nothing
        x_start, x_end = roi_x_pixels
        data = data[:, x_start:x_end, :]  # ROI on second dimension (width/x)
        println("  ROI: x=$roi_x_pixels -> $(size(data))")
    end

    # Create camera
    readnoise_array = fill(PARAMS.readnoise, size(data, 2), size(data, 1))
    camera = SCMOSCamera(size(data, 2), size(data, 1), PARAMS.pixel_size_um, readnoise_array;
        offset=PARAMS.offset_adu, gain=PARAMS.gain, qe=PARAMS.qe)

    # Detection
    print("  Detecting... ")
    t_det = @elapsed roi_batch = getboxes(data, camera;
        psf_sigma=PARAMS.psf_sigma,
        min_photons=PARAMS.min_photons,
        boxsize=PARAMS.boxsize,
        use_gpu=true)
    println("$(length(roi_batch)) ROIs in $(round(t_det, digits=1))s")

    # Fitting with sxsy model (same API as workflow_standard.jl)
    print("  Fitting sxsy... ")
    psf_model_type = GaussianXYNBSXSY()
    fitter = GaussMLEFitter(psf_model=psf_model_type, iterations=50, device=:gpu)
    t_fit = @elapsed smld_fit = GaussMLE.fit(fitter, roi_batch)
    println("done in $(round(t_fit, digits=1))s")

    # Extract data for filtering
    emitters = smld_fit.emitters
    psf_sx = [e.σx for e in emitters]
    psf_sy = [e.σy for e in emitters]
    sigma = sqrt.(psf_sx .* psf_sy)
    x = [e.x for e in emitters]
    pvalues = [e.pvalue for e in emitters]

    # Apply loose filter: just p-value to remove complete failures
    # (workflow applies strict filtering AFTER plotting the sxsy CDF)
    mask = pvalues .> PARAMS.pvalue_threshold

    n_pass = sum(mask)
    println("  P-value filter: $n_pass / $(length(mask)) ($(round(100*n_pass/length(mask), digits=1))%)")

    # Apply filter
    psf_sx_filt = psf_sx[mask]
    psf_sy_filt = psf_sy[mask]
    x_filt = x[mask]

    if is_2ch
        # Split by polarization at x=43μm
        split_x = 43.0
        left_mask = x_filt .< split_x
        right_mask = x_filt .>= split_x

        sx_left = psf_sx_filt[left_mask] .* 1000  # nm
        sy_left = psf_sy_filt[left_mask] .* 1000
        sx_right = psf_sx_filt[right_mask] .* 1000
        sy_right = psf_sy_filt[right_mask] .* 1000

        println("  Pol1 (left):  $(sum(left_mask)) locs, σx=$(round(median(sx_left), digits=1))nm, σy=$(round(median(sy_left), digits=1))nm")
        println("  Pol2 (right): $(sum(right_mask)) locs, σx=$(round(median(sx_right), digits=1))nm, σy=$(round(median(sy_right), digits=1))nm")

        return (pol1_sx=sx_left, pol1_sy=sy_left, pol2_sx=sx_right, pol2_sy=sy_right)
    else
        sx_nm = psf_sx_filt .* 1000
        sy_nm = psf_sy_filt .* 1000
        println("  1ch: $(length(sx_nm)) locs, σx=$(round(median(sx_nm), digits=1))nm, σy=$(round(median(sy_nm), digits=1))nm")
        return (sx=sx_nm, sy=sy_nm)
    end
end

# Main execution
println("\n" * "="^80)
println("COMBINED σx/σy CDF ANALYSIS")
println("="^80)

# 1-channel: right side ROI (pixels 430-860)
data_1ch = process_file(file_1ch, (430, 860), false)

# 2-channel: full width for polarization split
data_2ch = process_file(file_2ch, nothing, true)

# Create combined plot
println("\n" * "="^60)
println("Creating combined plot...")
println("="^60)

fig = Figure(size=(1000, 600))

ax = Axis(fig[1, 1],
    xlabel="PSF width (nm)",
    ylabel="CDF",
    title="PSF Width Distributions: σx (solid) vs σy (dashed)")

# Colors
c_1ch = :black
c_pol1 = :steelblue
c_pol2 = :darkorange

# Line styles
ls_sx = :solid
ls_sy = :dash

# Plot 1ch
cx, cy = ecdf_values(data_1ch.sx)
lines!(ax, cx, cy, color=c_1ch, linestyle=ls_sx, linewidth=2.5, label="1ch σx ($(round(median(data_1ch.sx), digits=0))nm)")
cx, cy = ecdf_values(data_1ch.sy)
lines!(ax, cx, cy, color=c_1ch, linestyle=ls_sy, linewidth=2.5, label="1ch σy ($(round(median(data_1ch.sy), digits=0))nm)")

# Plot 2ch pol1 (left - blue)
cx, cy = ecdf_values(data_2ch.pol1_sx)
lines!(ax, cx, cy, color=c_pol1, linestyle=ls_sx, linewidth=2.5, label="Pol1 σx ($(round(median(data_2ch.pol1_sx), digits=0))nm)")
cx, cy = ecdf_values(data_2ch.pol1_sy)
lines!(ax, cx, cy, color=c_pol1, linestyle=ls_sy, linewidth=2.5, label="Pol1 σy ($(round(median(data_2ch.pol1_sy), digits=0))nm)")

# Plot 2ch pol2 (right - orange)
cx, cy = ecdf_values(data_2ch.pol2_sx)
lines!(ax, cx, cy, color=c_pol2, linestyle=ls_sx, linewidth=2.5, label="Pol2 σx ($(round(median(data_2ch.pol2_sx), digits=0))nm)")
cx, cy = ecdf_values(data_2ch.pol2_sy)
lines!(ax, cx, cy, color=c_pol2, linestyle=ls_sy, linewidth=2.5, label="Pol2 σy ($(round(median(data_2ch.pol2_sy), digits=0))nm)")

# Legend
axislegend(ax, position=:rb, nbanks=2)

# Set reasonable limits
xlims!(ax, 50, 300)

save(outfile, fig)
println("\nSaved: $outfile")

# Print summary table
println("\n" * "="^60)
println("SUMMARY TABLE")
println("="^60)
println("Dataset      σx (nm)    σy (nm)    Ratio σy/σx")
println("-"^60)
med_1ch_sx = round(median(data_1ch.sx), digits=1)
med_1ch_sy = round(median(data_1ch.sy), digits=1)
med_pol1_sx = round(median(data_2ch.pol1_sx), digits=1)
med_pol1_sy = round(median(data_2ch.pol1_sy), digits=1)
med_pol2_sx = round(median(data_2ch.pol2_sx), digits=1)
med_pol2_sy = round(median(data_2ch.pol2_sy), digits=1)

println("1ch          $med_1ch_sx       $med_1ch_sy       $(round(med_1ch_sy/med_1ch_sx, digits=2))")
println("Pol1 (left)  $med_pol1_sx      $med_pol1_sy      $(round(med_pol1_sy/med_pol1_sx, digits=2))")
println("Pol2 (right) $med_pol2_sx      $med_pol2_sy      $(round(med_pol2_sy/med_pol2_sx, digits=2))")
println("="^60)

println("\nExpected pattern:")
println("- Pol1: narrow σx, wide σy  (ratio > 1)")
println("- Pol2: wide σx, narrow σy  (ratio < 1)")
println("- 1ch:  intermediate/mixed  (ratio ~ 1)")
