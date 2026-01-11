# Full Gatta Ruler Analysis
using SMLMAnalysis

# Load ALL data
h5file = "../data/gatta_ruler/2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5"
println("Loading data...")
data, info = smart_h5_to_array(h5file)  # All frames
println("  Loaded $(size(data, 3)) frames ($(size(data, 1))×$(size(data, 2)))")

# Camera setup
camera = SCMOSCamera(size(data, 2), size(data, 1), 0.1f0, 0.7f0;
    offset = 100.0f0, gain = 0.46f0, qe = 0.8f0)

# Run analysis - skip frame connection, enable render
result = analyze(data, camera;
    minval = 500.0,
    min_photons = 500.0,
    frameconnect = false,
    render = true,
    outdir = "output/gatta_full/"
)

println("\n" * "="^60)
println("DONE")
println("="^60)
println(result)
