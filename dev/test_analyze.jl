# Test the new analyze() function
using SMLMAnalysis

# Load test data
h5file = "../data/gatta_ruler/2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5"
println("Loading data...")
data, info = smart_h5_to_array(h5file; max_frames=500)
println("  Loaded $(size(data, 3)) frames ($(size(data, 1))×$(size(data, 2)))")

# Create camera (ORCA-Fusion specs)
# SCMOSCamera(width, height, pixelsize, readnoise; offset, gain, qe)
camera = SCMOSCamera(size(data, 2), size(data, 1), 0.1f0, 0.7f0;
    offset = 100.0f0, gain = 0.46f0, qe = 0.8f0)

# Test 1: All defaults
println("\n" * "="^60)
println("Test 1: All defaults")
println("="^60)
result = analyze(data, camera)
println(result)

# Test 2: With custom params and output
println("\n" * "="^60)
println("Test 2: Custom params + output")
println("="^60)
result = analyze(data, camera;
    minval = 500.0,
    min_photons = 500.0,
    render = true,
    outdir = "output/test_analyze/"
)
println(result)
println("\nTimings:")
for (k, v) in result.timings
    println("  $k: $(round(v, digits=2))s")
end

# Test 3: Config object
println("\n" * "="^60)
println("Test 3: Config object")
println("="^60)
config = AnalysisConfig(
    minval = 500.0,
    min_photons = 1000.0,
    frameconnect = false,
    render = true,
    render_zoom = 10,
    outdir = "output/test_analyze_config/"
)
result = analyze(data, camera, config)
println(result)

# Test 4: Load config from saved TOML
println("\n" * "="^60)
println("Test 4: Load config from TOML")
println("="^60)
loaded_config = load_config("output/test_analyze/config.toml")
println("Loaded config:")
println("  minval = $(loaded_config.minval)")
println("  min_photons = $(loaded_config.min_photons)")
println("  render = $(loaded_config.render)")

println("\n✅ All tests passed!")
