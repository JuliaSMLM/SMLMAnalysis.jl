# HDF5 I/O for BasicSMLD
# Native Julia format for SMLD serialization
#
# Supports emitter types from SMLMData and GaussMLE:
# - Emitter2DFit (standard 2D)
# - Emitter3DFit (standard 3D)
# - Emitter2DFitSigma (GaussMLE: isotropic fitted σ)
# - Emitter2DFitSigmaXY (GaussMLE: anisotropic fitted σx, σy)

using HDF5
using Dates
import Pkg

const SMLD_FORMAT_VERSION = "1.1"  # v1.1: Added PSF width fields for GaussMLE emitter types

# Get package version safely
function _get_package_version()
    try
        # Try to get version from the SMLMAnalysis package directly
        for (uuid, pkg) in Pkg.dependencies()
            if pkg.name == "SMLMAnalysis"
                return string(pkg.version)
            end
        end
        # Fallback: try project version
        proj = Pkg.project()
        if proj.version !== nothing
            return string(proj.version)
        end
    catch
    end
    return "unknown"
end

"""
    save_smld(filepath::String, smld::BasicSMLD;
              source_file::Union{String,Nothing}=nothing,
              drift_model=nothing,
              compression::Int=3)

Save BasicSMLD to HDF5 file with full metadata for reproducibility.

Supports all emitter types including GaussMLE types with fitted PSF widths:
- Emitter2DFitSigma: saves fitted σ (isotropic PSF width)
- Emitter2DFitSigmaXY: saves fitted σx, σy (anisotropic PSF widths)

# Arguments
- `filepath`: Output .h5 file path
- `smld`: BasicSMLD object to save
- `source_file`: Original data file path (for provenance tracking)
- `drift_model`: Optional drift correction model (LegendrePolynomial, etc.)
- `compression`: HDF5 compression level 0-9 (default: 3)

# File Structure
```
/metadata           - Format version, package info, timestamps
/emitters           - Columnar emitter data (x, y, z, photons, psf_sigma_x, etc.)
/camera             - Camera type and calibration
/drift_correction   - Drift model coefficients (if provided)
/provenance         - Source file info
```

# Example
```julia
save_smld("results.h5", smld; source_file="/data/experiment.h5", drift_model=dm)
```
"""
function save_smld(filepath::String, smld::BasicSMLD{T,E};
                   source_file::Union{String,Nothing}=nothing,
                   drift_model=nothing,
                   compression::Int=3) where {T,E}

    n = length(smld.emitters)
    is_3d = E <: Emitter3DFit
    emitter_type_name = string(nameof(E))

    # Detect emitter type features
    has_sigma = n > 0 && hasproperty(smld.emitters[1], :σ)      # Emitter2DFitSigma
    has_sigma_xy = n > 0 && hasproperty(smld.emitters[1], :σx)  # Emitter2DFitSigmaXY
    has_pvalue = n > 0 && hasproperty(smld.emitters[1], :pvalue)

    h5open(filepath, "w") do fid
        # === /metadata group ===
        meta = create_group(fid, "metadata")
        meta["format_version"] = SMLD_FORMAT_VERSION
        meta["package_name"] = "SMLMAnalysis"
        meta["package_version"] = _get_package_version()
        meta["emitter_type"] = emitter_type_name
        meta["element_type"] = string(T)
        meta["is_3d"] = is_3d
        meta["has_psf_sigma"] = has_sigma
        meta["has_psf_sigma_xy"] = has_sigma_xy
        meta["has_pvalue"] = has_pvalue
        meta["save_timestamp"] = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")
        meta["n_frames"] = smld.n_frames
        meta["n_datasets"] = smld.n_datasets
        meta["n_emitters"] = n

        # === /emitters group ===
        em = create_group(fid, "emitters")

        if n > 0
            # Core fields (all emitter types)
            x = [e.x for e in smld.emitters]
            y = [e.y for e in smld.emitters]
            photons = [e.photons for e in smld.emitters]
            bg = [e.bg for e in smld.emitters]
            σ_x = [e.σ_x for e in smld.emitters]
            σ_y = [e.σ_y for e in smld.emitters]
            σ_photons = [e.σ_photons for e in smld.emitters]
            σ_bg = [e.σ_bg for e in smld.emitters]
            frame = Int32[e.frame for e in smld.emitters]
            dataset = Int32[e.dataset for e in smld.emitters]
            track_id = Int32[e.track_id for e in smld.emitters]
            id = Int32[e.id for e in smld.emitters]

            # Write core fields
            em["x", compress=compression] = x
            em["y", compress=compression] = y
            em["photons", compress=compression] = photons
            em["bg", compress=compression] = bg
            em["sigma_x", compress=compression] = σ_x
            em["sigma_y", compress=compression] = σ_y
            em["sigma_photons", compress=compression] = σ_photons
            em["sigma_bg", compress=compression] = σ_bg
            em["frame", compress=compression] = frame
            em["dataset", compress=compression] = dataset
            em["track_id", compress=compression] = track_id
            em["id", compress=compression] = id

            # 3D fields
            if is_3d
                z = [e.z for e in smld.emitters]
                σ_z = [e.σ_z for e in smld.emitters]
                em["z", compress=compression] = z
                em["sigma_z", compress=compression] = σ_z
            end

            # GaussMLE Emitter2DFitSigma fields (isotropic PSF)
            if has_sigma
                psf_sigma = [e.σ for e in smld.emitters]
                em["psf_sigma", compress=compression] = psf_sigma

                if hasproperty(smld.emitters[1], :σ_σ)
                    σ_sigma = [e.σ_σ for e in smld.emitters]
                    em["sigma_psf_sigma", compress=compression] = σ_sigma
                end
            end

            # GaussMLE Emitter2DFitSigmaXY fields (anisotropic PSF)
            if has_sigma_xy
                psf_sigma_x = [e.σx for e in smld.emitters]
                psf_sigma_y = [e.σy for e in smld.emitters]
                em["psf_sigma_x", compress=compression] = psf_sigma_x
                em["psf_sigma_y", compress=compression] = psf_sigma_y

                if hasproperty(smld.emitters[1], :σ_σx)
                    σ_sigma_x = [e.σ_σx for e in smld.emitters]
                    σ_sigma_y = [e.σ_σy for e in smld.emitters]
                    em["sigma_psf_sigma_x", compress=compression] = σ_sigma_x
                    em["sigma_psf_sigma_y", compress=compression] = σ_sigma_y
                end
            end

            # p-value (GaussMLE emitters)
            if has_pvalue
                pvalue = [e.pvalue for e in smld.emitters]
                em["pvalue", compress=compression] = pvalue
            end
        end

        # === /camera group ===
        cam = create_group(fid, "camera")
        camera = smld.camera

        cam["type"] = camera isa IdealCamera ? "IdealCamera" : "SCMOSCamera"
        cam["pixel_edges_x", compress=compression] = collect(camera.pixel_edges_x)
        cam["pixel_edges_y", compress=compression] = collect(camera.pixel_edges_y)

        if camera isa SCMOSCamera
            # Store calibration data (handle scalars vs arrays)
            if camera.offset isa AbstractArray
                cam["offset", compress=compression] = collect(camera.offset)
            else
                cam["offset"] = camera.offset
            end
            if camera.gain isa AbstractArray
                cam["gain", compress=compression] = collect(camera.gain)
            else
                cam["gain"] = camera.gain
            end
            if camera.readnoise isa AbstractArray
                cam["readnoise", compress=compression] = collect(camera.readnoise)
            else
                cam["readnoise"] = camera.readnoise
            end
            if camera.qe isa AbstractArray
                cam["qe", compress=compression] = collect(camera.qe)
            else
                cam["qe"] = camera.qe
            end
        end

        # === /provenance group ===
        prov = create_group(fid, "provenance")
        if source_file !== nothing
            prov["source_file"] = source_file
        end

        # === /drift_correction group (optional) ===
        if drift_model !== nothing
            dc = create_group(fid, "drift_correction")
            _save_drift_model!(dc, drift_model, compression)
        end

        # === /user_metadata group ===
        if !isempty(smld.metadata)
            user_meta = create_group(fid, "user_metadata")
            for (key, value) in smld.metadata
                try
                    # Only save HDF5-compatible types
                    if value isa Union{String, Number, AbstractArray{<:Number}}
                        user_meta[key] = value
                    elseif value isa AbstractArray{<:String}
                        user_meta[key] = collect(value)
                    end
                    # Skip unsupported types silently
                catch
                    # Skip values that can't be serialized
                end
            end
        end
    end

    return filepath
end

"""
    _save_drift_model!(group, dm::LegendrePolynomial, compression)

Save LegendrePolynomial drift model to HDF5 group.
"""
function _save_drift_model!(group, dm, compression::Int)
    # Generic handling - store type name
    group["model_type"] = string(typeof(dm).name.name)

    # Check for LegendrePolynomial-like structure
    if hasproperty(dm, :ndatasets) && hasproperty(dm, :intra) && hasproperty(dm, :inter)
        group["ndatasets"] = dm.ndatasets

        # Get parameters from first intra model
        if !isempty(dm.intra)
            first_intra = dm.intra[1]
            group["ndims"] = first_intra.ndims
            if !isempty(first_intra.dm)
                group["degree"] = first_intra.dm[1].degree
                group["n_frames"] = first_intra.dm[1].n_frames
            end
        end

        # Store intra coefficients: (ndatasets, ndims, degree)
        ndatasets = dm.ndatasets
        ndims = dm.intra[1].ndims
        degree = dm.intra[1].dm[1].degree

        intra_coeffs = zeros(ndatasets, ndims, degree)
        for d in 1:ndatasets
            for dim in 1:ndims
                intra_coeffs[d, dim, :] = dm.intra[d].dm[dim].coefficients
            end
        end
        group["intra_coefficients", compress=compression] = intra_coeffs

        # Store inter shifts: (ndatasets, ndims)
        inter_shifts = zeros(ndatasets, ndims)
        for d in 1:ndatasets
            for dim in 1:ndims
                inter_shifts[d, dim] = dm.inter[d].dm[dim]
            end
        end
        group["inter_shifts", compress=compression] = inter_shifts
    end
end

"""
    load_smld(filepath::String) -> BasicSMLD

Load BasicSMLD from HDF5 file.

Returns a BasicSMLD with the saved emitters, camera, and metadata.
Automatically reconstructs the correct emitter type (including GaussMLE types
with PSF width fields).

Drift correction info is stored in metadata["drift_correction"] if present.

# Example
```julia
smld = load_smld("results.h5")
```
"""
function load_smld(filepath::String)
    h5open(filepath, "r") do fid
        # Read metadata
        meta = fid["metadata"]
        format_version = read(meta["format_version"])
        emitter_type_str = read(meta["emitter_type"])
        element_type_str = read(meta["element_type"])
        is_3d = read(meta["is_3d"])
        n_frames = read(meta["n_frames"])
        n_datasets = read(meta["n_datasets"])
        n_emitters = read(meta["n_emitters"])

        # Feature flags (may not exist in v1.0 files)
        has_psf_sigma = haskey(meta, "has_psf_sigma") ? read(meta["has_psf_sigma"]) : false
        has_psf_sigma_xy = haskey(meta, "has_psf_sigma_xy") ? read(meta["has_psf_sigma_xy"]) : false
        has_pvalue = haskey(meta, "has_pvalue") ? read(meta["has_pvalue"]) : false

        # Determine numeric type
        T = element_type_str == "Float32" ? Float32 : Float64

        # Read emitters
        em = fid["emitters"]

        if n_emitters > 0
            # Core fields
            x = read(em["x"])
            y = read(em["y"])
            photons = read(em["photons"])
            bg = read(em["bg"])
            σ_x = read(em["sigma_x"])
            σ_y = read(em["sigma_y"])
            σ_photons = read(em["sigma_photons"])
            σ_bg = read(em["sigma_bg"])
            frame = read(em["frame"])
            dataset = read(em["dataset"])
            track_id = read(em["track_id"])
            id = read(em["id"])

            # Optional 3D fields
            z = is_3d ? read(em["z"]) : nothing
            σ_z = is_3d ? read(em["sigma_z"]) : nothing

            # Optional PSF sigma fields (check both metadata flag and dataset existence)
            psf_sigma = (has_psf_sigma || haskey(em, "psf_sigma")) ? read(em["psf_sigma"]) : nothing
            σ_sigma = haskey(em, "sigma_psf_sigma") ? read(em["sigma_psf_sigma"]) : nothing

            # Optional PSF sigma_xy fields
            psf_sigma_x = (has_psf_sigma_xy || haskey(em, "psf_sigma_x")) ? read(em["psf_sigma_x"]) : nothing
            psf_sigma_y = (has_psf_sigma_xy || haskey(em, "psf_sigma_y")) ? read(em["psf_sigma_y"]) : nothing
            σ_sigma_x = haskey(em, "sigma_psf_sigma_x") ? read(em["sigma_psf_sigma_x"]) : nothing
            σ_sigma_y = haskey(em, "sigma_psf_sigma_y") ? read(em["sigma_psf_sigma_y"]) : nothing

            # Optional p-value
            pvalue = (has_pvalue || haskey(em, "pvalue")) ? read(em["pvalue"]) : nothing

            # Construct emitters based on type
            emitters = _construct_emitters(
                emitter_type_str, T, n_emitters, is_3d,
                x, y, z, photons, bg,
                σ_x, σ_y, σ_z, σ_photons, σ_bg,
                psf_sigma, σ_sigma,
                psf_sigma_x, psf_sigma_y, σ_sigma_x, σ_sigma_y,
                pvalue,
                frame, dataset, track_id, id
            )
        else
            emitters = _empty_emitters(emitter_type_str, T, is_3d)
        end

        # Read camera
        cam_grp = fid["camera"]
        cam_type = read(cam_grp["type"])
        pixel_edges_x = read(cam_grp["pixel_edges_x"])
        pixel_edges_y = read(cam_grp["pixel_edges_y"])

        if cam_type == "IdealCamera"
            camera = IdealCamera(pixel_edges_x, pixel_edges_y)
        else
            # SCMOSCamera
            offset = read(cam_grp["offset"])
            gain = read(cam_grp["gain"])
            readnoise = read(cam_grp["readnoise"])
            qe = read(cam_grp["qe"])
            camera = SCMOSCamera(pixel_edges_x, pixel_edges_y;
                                 offset=offset, gain=gain,
                                 readnoise=readnoise, qe=qe)
        end

        # Build metadata dict
        metadata = Dict{String,Any}()

        # Read provenance
        if haskey(fid, "provenance")
            prov = fid["provenance"]
            if haskey(prov, "source_file")
                metadata["source_file"] = read(prov["source_file"])
            end
        end

        # Store format info
        metadata["smld_format_version"] = format_version
        metadata["saved_package_version"] = haskey(meta, "package_version") ?
            read(meta["package_version"]) : "unknown"
        metadata["save_timestamp"] = haskey(meta, "save_timestamp") ?
            read(meta["save_timestamp"]) : "unknown"

        # Read drift correction if present
        if haskey(fid, "drift_correction")
            dc = fid["drift_correction"]
            drift_info = Dict{String,Any}()
            drift_info["model_type"] = read(dc["model_type"])
            if haskey(dc, "degree")
                drift_info["degree"] = read(dc["degree"])
            end
            if haskey(dc, "n_frames")
                drift_info["n_frames"] = read(dc["n_frames"])
            end
            if haskey(dc, "intra_coefficients")
                drift_info["intra_coefficients"] = read(dc["intra_coefficients"])
            end
            if haskey(dc, "inter_shifts")
                drift_info["inter_shifts"] = read(dc["inter_shifts"])
            end
            metadata["drift_correction"] = drift_info
        end

        # Read user metadata
        if haskey(fid, "user_metadata")
            user_meta = fid["user_metadata"]
            for key in keys(user_meta)
                metadata[key] = read(user_meta[key])
            end
        end

        return BasicSMLD(emitters, camera, n_frames, n_datasets, metadata)
    end
end

"""
    _construct_emitters(...)

Internal: Construct emitters of the appropriate type based on emitter_type_str.
"""
function _construct_emitters(
    emitter_type_str::String, T::Type, n::Int, is_3d::Bool,
    x, y, z, photons, bg,
    σ_x, σ_y, σ_z, σ_photons, σ_bg,
    psf_sigma, σ_sigma,
    psf_sigma_x, psf_sigma_y, σ_sigma_x, σ_sigma_y,
    pvalue,
    frame, dataset, track_id, id
)
    # Try to get GaussMLE emitter types if available
    GaussMLEEmitterTypes = _get_gaussmle_emitter_types()

    if emitter_type_str == "Emitter2DFitSigmaXY" && psf_sigma_x !== nothing && GaussMLEEmitterTypes !== nothing
        # GaussMLE anisotropic PSF type
        Emitter2DFitSigmaXY = GaussMLEEmitterTypes.Emitter2DFitSigmaXY
        return [Emitter2DFitSigmaXY{T}(
            T(x[i]), T(y[i]),
            T(photons[i]), T(bg[i]),
            T(psf_sigma_x[i]), T(psf_sigma_y[i]),
            T(σ_x[i]), T(σ_y[i]), T(σ_photons[i]), T(σ_bg[i]),
            σ_sigma_x !== nothing ? T(σ_sigma_x[i]) : T(0),
            σ_sigma_y !== nothing ? T(σ_sigma_y[i]) : T(0),
            pvalue !== nothing ? T(pvalue[i]) : T(0),
            Int(frame[i]), Int(dataset[i]), Int(track_id[i]), Int(id[i])
        ) for i in 1:n]

    elseif emitter_type_str == "Emitter2DFitSigma" && psf_sigma !== nothing && GaussMLEEmitterTypes !== nothing
        # GaussMLE isotropic PSF type
        Emitter2DFitSigma = GaussMLEEmitterTypes.Emitter2DFitSigma
        return [Emitter2DFitSigma{T}(
            T(x[i]), T(y[i]),
            T(photons[i]), T(bg[i]),
            T(psf_sigma[i]),
            T(σ_x[i]), T(σ_y[i]), T(σ_photons[i]), T(σ_bg[i]),
            σ_sigma !== nothing ? T(σ_sigma[i]) : T(0),
            pvalue !== nothing ? T(pvalue[i]) : T(0),
            Int(frame[i]), Int(dataset[i]), Int(track_id[i]), Int(id[i])
        ) for i in 1:n]

    elseif is_3d
        # Standard 3D emitter
        return [Emitter3DFit{T}(
            T(x[i]), T(y[i]), T(z[i]),
            T(photons[i]), T(bg[i]),
            T(σ_x[i]), T(σ_y[i]), T(σ_z[i]),
            T(σ_photons[i]), T(σ_bg[i]);
            frame=Int(frame[i]),
            dataset=Int(dataset[i]),
            track_id=Int(track_id[i]),
            id=Int(id[i])
        ) for i in 1:n]

    else
        # Standard 2D emitter (fallback)
        return [Emitter2DFit{T}(
            T(x[i]), T(y[i]),
            T(photons[i]), T(bg[i]),
            T(σ_x[i]), T(σ_y[i]),
            T(σ_photons[i]), T(σ_bg[i]);
            frame=Int(frame[i]),
            dataset=Int(dataset[i]),
            track_id=Int(track_id[i]),
            id=Int(id[i])
        ) for i in 1:n]
    end
end

"""
    _empty_emitters(emitter_type_str, T, is_3d)

Internal: Create empty emitter vector of appropriate type.
"""
function _empty_emitters(emitter_type_str::String, T::Type, is_3d::Bool)
    GaussMLEEmitterTypes = _get_gaussmle_emitter_types()

    if emitter_type_str == "Emitter2DFitSigmaXY" && GaussMLEEmitterTypes !== nothing
        return GaussMLEEmitterTypes.Emitter2DFitSigmaXY{T}[]
    elseif emitter_type_str == "Emitter2DFitSigma" && GaussMLEEmitterTypes !== nothing
        return GaussMLEEmitterTypes.Emitter2DFitSigma{T}[]
    elseif is_3d
        return Emitter3DFit{T}[]
    else
        return Emitter2DFit{T}[]
    end
end

"""
    _get_gaussmle_emitter_types()

Internal: Try to get GaussMLE emitter types if the package is loaded.
Returns nothing if GaussMLE is not available.
"""
function _get_gaussmle_emitter_types()
    try
        # Check if GaussMLE is loaded
        if isdefined(Main, :GaussMLE)
            return Main.GaussMLE
        end
        # Check if loaded in SMLMAnalysis
        if isdefined(@__MODULE__, :GaussMLE)
            return GaussMLE
        end
        # Try to get from package extensions or parent module
        for m in values(Base.loaded_modules)
            if nameof(m) == :GaussMLE
                return m
            end
        end
    catch
    end
    return nothing
end

"""
    smld_info(filepath::String)

Print summary info about an SMLD HDF5 file without loading all data.
"""
function smld_info(filepath::String)
    h5open(filepath, "r") do fid
        meta = fid["metadata"]

        println("SMLD File: $filepath")
        println("  Format version: ", read(meta["format_version"]))
        println("  Package: ", read(meta["package_name"]), " v",
                haskey(meta, "package_version") ? read(meta["package_version"]) : "?")
        println("  Emitter type: ", read(meta["emitter_type"]))
        println("  Element type: ", read(meta["element_type"]))
        println("  Emitters: ", read(meta["n_emitters"]))
        println("  Frames: ", read(meta["n_frames"]))
        println("  Datasets: ", read(meta["n_datasets"]))

        # PSF fields info
        if haskey(meta, "has_psf_sigma") && read(meta["has_psf_sigma"])
            println("  PSF model: isotropic (σ)")
        elseif haskey(meta, "has_psf_sigma_xy") && read(meta["has_psf_sigma_xy"])
            println("  PSF model: anisotropic (σx, σy)")
        end

        if haskey(meta, "save_timestamp")
            println("  Saved: ", read(meta["save_timestamp"]))
        end

        if haskey(fid, "provenance") && haskey(fid["provenance"], "source_file")
            println("  Source: ", read(fid["provenance"]["source_file"]))
        end

        if haskey(fid, "drift_correction")
            dc = fid["drift_correction"]
            println("  Drift model: ", read(dc["model_type"]))
            if haskey(dc, "degree")
                println("    Degree: ", read(dc["degree"]))
            end
        end

        # Show available emitter fields
        if haskey(fid, "emitters")
            em = fid["emitters"]
            em_keys = sort(collect(keys(em)))
            println("  Emitter fields: ", join(em_keys, ", "))
        end
    end
end
