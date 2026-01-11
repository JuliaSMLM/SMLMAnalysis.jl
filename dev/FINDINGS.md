# DNA-PAINT Ruler Data Analysis - Findings

## Camera Identification

**Camera**: Hamamatsu ORCA-Fusion (C14440-20UP)
**Interface**: DCAM4

### Specifications (from metadata and Hamamatsu docs)
- **Sensor**: 2304 x 2304 pixels
- **Physical pixel size**: 6.5 µm
- **Readnoise**: 0.7 - 1.4 electrons RMS (scan mode dependent)
  - Ultra quiet: 0.7 e-
  - Standard: 1.0 e-
  - Fast: 1.4 e-
- **Quantum Efficiency**: 80% (peak)
- **Full Well**: 15,000 electrons
- **Dark current**: 0.2 - 0.5 e-/pixel/s (cooling dependent)
- **ADC**: 16-bit / 12-bit / 8-bit (selectable)

### Actual Acquisition Settings (from H5 metadata)
- **ROI**: 860 x 256 pixels
- **ROI offset**: x=516, y=1480
- **Exposure time**: 100.014 ms
- **Frame rate**: 10 fps
- **Sequence length**: 20,000 frames
- **Trigger mode**: 1 (external trigger)
- **Gain**: 1.0
- **Pixel format gain**: 1.0
- **Laser**: 642nm at 100% power (DNA-PAINT)
- **Mode**: TIRF illumination

## Data Import

Successfully created H5 import functions in `src/import_smart_h5.jl`:

- `load_smart_h5_info(filepath)` - Get metadata without loading data
- `load_smart_h5_frame(filepath, frame)` - Load single frame
- `load_smart_h5(filepath; frame_range)` - Load specific frames
- `smart_h5_to_array(filepath; max_frames)` - Load and transpose for processing

### H5 File Structure
```
/Main/data                    - (860, 256, 20000) UInt16 array
/Main/camera                  - Empty group, but has attributes:
  - camera_format_x_pixels: 2304
  - camera_format_y_pixels: 2304
  - roi_width: 860
  - roi_height: 256
  - exposure_time: 0.100014
  - frame_rate: 10.0
  - unique_id: "DCAM4Camera"
/Main/laser_XXX               - Laser control groups
/Main/stage_XXX               - Stage position groups
```

## Known Issues

### 1. Pixel Size Calibration Required

The metadata shows `camera_format_pixelsize: 1.0` which is likely a placeholder or already accounts for magnification.

For typical TIRF-SMLM setups:
- Objective: 100x oil immersion
- Additional magnification: ~1.5x relay
- **Total magnification**: ~150x
- **Calculated effective pixel size**: 6.5 µm / 150 = **43.3 nm**

However, without calibration data, this is an estimate. **Action required**: Verify actual magnification from microscope calibration files.

### 2. SCMOSCamera Type Mismatch Bug

**Bug Location**: `SMLMBoxer/src/filter.jl:213` in `convolve_variance_weighted()`

**Issue**: Function signature expects:
```julia
convolve_variance_weighted(::AbstractArray{T}, ::AbstractMatrix{T}, ::Float32, ::Bool) where T<:Real
```

But gets called with:
```julia
convolve_variance_weighted(::Array{Float32, 4}, ::Matrix{Float64}, ::Float32, ::Bool)
```

**Root Cause**: SCMOSCamera stores readnoise/variance as Float64, but imagestack is converted to Float32, causing type mismatch.

**Workarounds**:
1. Use `IdealCamera` instead (works correctly)
2. Fix SMLMBoxer to handle mixed Float32/Float64 types
3. Ensure SCMOSCamera parameters are Float32

**Status**: Reported in dev/gatta_rulers_test.jl comments. Needs fix in SMLMBoxer.

### 3. Background Processing Status

Several background Julia processes are still running from previous tests:
- 6dbd60, 7a8c80, 27234e, 52a442, 3ff9fe, 0fa147, cbea5c, 20b91c, 40ba4e

**Action**: Clean up these processes or let them complete.

## Success: Example Workflows

The `examples/step_by_step_workflow.jl` completed successfully with:
- **Detection**: 10,623 ROIs from 12,014 ground truth emitters
- **Fitting**: 10,564 emitters after precision filtering
- **Jaccard Index**: 0.609 (100nm cutoff)
- **Mean precision**: 3.8 nm (both x and y)

This demonstrates the pipeline works with IdealCamera.

## Next Steps

### Immediate
1. Fix SCMOSCamera type mismatch in SMLMBoxer
2. Verify pixel size calibration
3. Re-run `gatta_rulers_test.jl` with corrected camera

### Analysis
4. Apply frame connection (SMLMFrameConnection.jl)
5. Measure ruler distances (expect ~20nm spacing for "20R" ruler)
6. Generate super-resolution image with SMLMRender
7. Compare measured vs. expected distances

### Code Integration
8. Move import functions from dev/ to src/ once tested
9. Add unit tests for H5 import
10. Document camera setup patterns in CLAUDE.md

## Files Created

- `dev/gatta_rulers_test.jl` - Main analysis script (needs SCMOSCamera fix)
- `src/import_smart_h5.jl` - H5 import functions (working)
- `dev/explore_h5.jl` - H5 structure exploration
- `dev/explore_camera_data.jl` - Camera metadata extraction
- `dev/test_smart_import.jl` - Import function tests (passing)
- `dev/FINDINGS.md` - This document

## Data Location

- **Symlink**: `data/gatta_ruler` → `/mnt/nas/adapt/projects/smart-microscope/data/DNA paint ruler`
- **Test dataset**: `2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5`
- **Size**: 8.81 GB (20,000 frames)
