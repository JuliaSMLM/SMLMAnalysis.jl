"""
    helpers.jl

Helper functions to connect outputs between ecosystem packages.
These functions handle data transformations and format conversions.
"""

using SMLMData
using SMLMBoxer
using GaussMLE

"""
    boxer_to_roi_batch(boxer_result, camera) → SMLMData.ROIBatch

Convert SMLMBoxer getboxes() output to GaussMLE ROIBatch format.

# Arguments
- `boxer_result`: NamedTuple from `getboxes()` containing `boxes`, `boxcoords`, etc.
- `camera`: AbstractCamera used in detection

# Returns
SMLMData.ROIBatch ready for GaussMLE fitting

# Example
```julia
boxer_result = getboxes(images, camera)
roi_batch = boxer_to_roi_batch(boxer_result, camera)
fitter = GaussMLEFitter()
results = fit(fitter, roi_batch)
```
"""
function boxer_to_roi_batch(boxer_result, camera::SMLMData.AbstractCamera)
    # Extract data from boxer result
    boxes = boxer_result.boxes  # (boxsize, boxsize, nboxes)
    boxcoords = boxer_result.boxcoords  # N×3 matrix (row, col, frame)

    # GaussMLE expects x_corners, y_corners (col, row) not (row, col)
    # boxcoords from boxer is (row, col, frame)
    x_corners = Int32.(boxcoords[:, 2])  # col → x
    y_corners = Int32.(boxcoords[:, 1])  # row → y
    frame_indices = Int32.(boxcoords[:, 3])  # frame

    # Create ROIBatch
    return SMLMData.ROIBatch(boxes, x_corners, y_corners, frame_indices, camera)
end

"""
    localization_result_to_smld(result::GaussMLE.LocalizationResult, roi_batch::SMLMData.ROIBatch;
                                dataset::Int=1, metadata::Dict=Dict()) → SMLMData.BasicSMLD

Convert GaussMLE LocalizationResult to SMLMData BasicSMLD format.

This is a convenience wrapper around GaussMLE.to_smld().

# Arguments
- `result`: LocalizationResult from GaussMLE fit
- `roi_batch`: Original ROIBatch used for fitting (contains camera info)
- `dataset`: Dataset number (default: 1)
- `metadata`: Additional metadata to include

# Returns
SMLMData.BasicSMLD with fitted localizations in micron coordinates
"""
function localization_result_to_smld(result::GaussMLE.LocalizationResult,
                                     roi_batch::SMLMData.ROIBatch;
                                     dataset::Int=1,
                                     metadata::Dict{String,Any}=Dict{String,Any}())
    return GaussMLE.to_smld(result, roi_batch; dataset=dataset, metadata=metadata)
end

"""
    summarize_boxer_result(result::SMLMData.ROIBatch) → String

Create a summary string of SMLMBoxer detection result (ROIBatch).
"""
function summarize_boxer_result(result::SMLMData.ROIBatch)
    n_detections = length(result)
    boxsize = result.roi_size
    return "Detected $n_detections particles (boxsize=$(boxsize)×$(boxsize))"
end

"""
    summarize_fit_result(result::LocalizationResult) → String

Create a summary string of GaussMLE fit result.
"""
function summarize_fit_result(result::GaussMLE.LocalizationResult)
    n_fits = result.n_fits
    psf_model = typeof(result.psf_model).name.name
    return "Fitted $n_fits localizations using $psf_model"
end

"""
    summarize_smld(smld::BasicSMLD) → String

Create a summary string of an SMLD object.
"""
function summarize_smld(smld::SMLMData.BasicSMLD)
    n_emitters = length(smld.emitters)
    n_frames = smld.n_frames
    return "BasicSMLD with $n_emitters emitters across $n_frames frames"
end
