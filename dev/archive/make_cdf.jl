using SMLMAnalysis, CairoMakie, Statistics

h5file = "../data/gatta_ruler/2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5"
data, _ = smart_h5_to_array(h5file, max_frames=100)
pixels = vec(data)
sorted_pixels = sort(pixels)
cdf_values = (1:length(sorted_pixels)) ./ length(sorted_pixels)

fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1],
    xlabel="Intensity (ADU)",
    ylabel="Cumulative Probability",
    title="Intensity CDF (first 100 frames)",
    xscale=log10)
lines!(ax, sorted_pixels, cdf_values, linewidth=2)

percentiles = [0.01, 0.25, 0.50, 0.75, 0.99]
for p in percentiles
    val = quantile(pixels, p)
    vlines!(ax, [val], color=:red, linestyle=:dash, alpha=0.5)
    text!(ax, val, p, text="$(Int(round(p*100)))%", align=(:left, :center))
end

save("output/01_raw/intensity_cdf.png", fig)
println("Saved output/01_raw/intensity_cdf.png")
