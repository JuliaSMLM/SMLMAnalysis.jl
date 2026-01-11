using SMLMAnalysis

# Test file path
h5file = "data/gatta_ruler/2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5"

println("="^80)
println("Testing SMART microscope H5 data import")
println("="^80)

# Test 1: Load file info
println("\n1. Loading file info...")
info = load_smart_h5_info(h5file)
println("   File: ", basename(info.filepath))
println("   Dimensions: ", info.width, " x ", info.height, " x ", info.nframes, " frames")
println("   Data type: ", info.dtype)
println("   File size: ", round(info.file_size_gb, digits=2), " GB")

# Test 2: Load a single frame
println("\n2. Loading single frame (frame 1)...")
frame1 = load_smart_h5_frame(h5file, 1)
println("   Frame shape: ", size(frame1))
println("   Min/Max values: ", minimum(frame1), " / ", maximum(frame1))

# Test 3: Load a range of frames
println("\n3. Loading frame range (1:100)...")
data_subset = load_smart_h5(h5file, frame_range=1:100)
println("   Data shape: ", size(data_subset))
println("   Memory usage: ", round(sizeof(data_subset) / 1e6, digits=2), " MB")

# Test 4: Load and transpose data for SMLM processing
println("\n4. Loading data in standard format (first 500 frames)...")
data, info = smart_h5_to_array(h5file, max_frames=500)
println("   Data shape (rows, cols, frames): ", size(data))
println("   Memory usage: ", round(sizeof(data) / 1e6, digits=2), " MB")
println("   Min/Max values: ", minimum(data), " / ", maximum(data))

println("\n" * "="^80)
println("All tests passed!")
println("="^80)
