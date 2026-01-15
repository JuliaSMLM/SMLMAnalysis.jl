"""
    config.jl

Configuration types for SMLM analysis workflows.
"""

using TOML

# =============================================================================
# AnalysisConfig - Configuration for SMLM analysis workflow
# =============================================================================

"""
    AnalysisConfig

Configuration struct for SMLM analysis workflow. Fields are ordered by pipeline flow.
Use `true` to enable a step, `false` to skip it.

# Example
```julia
# All defaults
config = AnalysisConfig()

# Custom configuration
config = AnalysisConfig(
    detect_min_photons = 500.0,  # Detection threshold
    min_photons = 500.0,         # Filter threshold
    render = true,
    outdir = "output/"
)

# Save/load config
save_config("myconfig.toml", config)
config = load_config("myconfig.toml")
```
"""
Base.@kwdef struct AnalysisConfig
    # === Detection (always runs) ===
    boxsize::Int = 11
    overlap::Float64 = 2.0
    detect_min_photons::Float64 = 500.0  # Detection threshold in photons (auto-converted to DoG threshold)
    detect_gpu::Bool = true

    # === PSF sigma (used for both detection and fixed-sigma fitting) ===
    # Detection: DoG uses sigma_small = psf_sigma, sigma_large = 2*psf_sigma (in pixels)
    # Fitting: Only used when fit_model=:fixed (GaussianXYNB)
    psf_sigma::Float32 = 0.135f0    # PSF sigma in MICRONS (~135nm typical for TIRF)

    # === Fitting (always runs) ===
    # PSF model: :fixed (GaussianXYNB), :variable (GaussianXYNBS), :anisotropic (GaussianXYNBSXSY)
    fit_model::Symbol = :variable   # Default to variable sigma for better pvalue
    iterations::Int = 20
    fit_device::Symbol = :auto      # :cpu, :gpu, :auto

    # === Filtering (runs after fitting) ===
    filter::Bool = true
    min_photons::Union{Float64, Nothing} = 500.0
    max_precision::Union{Float64, Nothing} = 0.015  # Max localization precision in μm (15nm default)
    psf_sigma_mode_tolerance::Union{Float64, Nothing} = 0.10  # Keep PSF sigma within ±10% of mode (for variable sigma fits)
    min_pvalue::Union{Float64, Nothing} = 1e-3    # p-value threshold

    # === Frame Connection (runs after filtering) ===
    frameconnect::Bool = false      # OFF by default
    fc_maxframegap::Int = 5         # Max frame gap between connected locs
    fc_nsigmadev::Float64 = 5.0     # Sigma multiplier for preclustering distance
    fc_nnearestclusters::Int = 2    # Nearest clusters for density estimation
    fc_nmaxnn::Int = 2              # Max nearest neighbors for preclustering

    # === Uncertainty Calibration (runs after frame connection) ===
    calibrate_uncertainties::Bool = true  # Adjust uncertainties using k and σ_motion from frame connection

    # === Drift Correction (runs after uncertainty calibration) ===
    drift::Bool = true              # ON by default for DNA-PAINT
    drift_degree::Int = 2           # Polynomial degree (2 usually sufficient)
    drift_cost_fun::String = "Kdtree"  # "Kdtree" (fast) or "Entropy"
    drift_model::String = "Polynomial" # "Polynomial" or "LegendrePoly"
    dataset_indices::Union{Nothing, Vector{Int}} = nothing  # Map frame → dataset ID for multi-dataset drift correction

    # === Isolated Emitter Filter (runs after drift correction) ===
    filter_isolated::Bool = false   # OFF by default
    isolated_n_sigma::Float64 = 2.0 # Neighbor if dist < n_sigma * sqrt(σ_i² + σ_j²)
    isolated_min_neighbors::Union{Int, Symbol} = :auto  # :auto uses triangle method, or set Int manually

    # === Rendering ===
    render::Bool = true             # ON by default - primary output
    render_gaussian::Bool = true    # Gaussian blur @ 20x + inferno
    render_histogram::Bool = true   # Histogram @ 10x + time coloring
    render_circles::Bool = true     # Circles @ 50x + time coloring
    render_gaussian_zoom::Int = 20
    render_histogram_zoom::Int = 10
    render_circles_zoom::Int = 50
    render_time_colormap::Symbol = :turbo  # colormap for time-colored renders
    render_clip_percentile::Union{Float64, Symbol} = :auto  # :auto adapts to n_locs, or set 0.0-1.0

    # === BAGOL (Deep Learning) ===
    bagol::Bool = false             # OFF by default
    bagol_model::String = "default"

    # === Output ===
    outdir::Union{String, Nothing} = nothing
    save_figures::Bool = true
    save_smld::Bool = true
end

"""
    save_config(filepath, config::AnalysisConfig)

Save an AnalysisConfig to a TOML file.
"""
function save_config(filepath::String, config::AnalysisConfig)
    d = Dict{String, Any}()
    for field in fieldnames(AnalysisConfig)
        val = getfield(config, field)
        if val !== nothing
            # Convert Symbol to String for TOML compatibility
            d[string(field)] = val isa Symbol ? string(val) : val
        end
    end
    open(filepath, "w") do io
        TOML.print(io, d)
    end
    return filepath
end

"""
    load_config(filepath) -> AnalysisConfig

Load an AnalysisConfig from a TOML file.
"""
function load_config(filepath::String)
    d = TOML.parsefile(filepath)
    kwargs = Dict{Symbol, Any}()
    # Fields that should be Symbols
    symbol_fields = (:fit_device, :fit_model, :render_time_colormap, :render_clip_percentile, :isolated_min_neighbors)
    for (k, v) in d
        sym = Symbol(k)
        # Convert string back to Symbol for appropriate fields
        if sym in symbol_fields && v isa String
            kwargs[sym] = Symbol(v)
        else
            kwargs[sym] = v
        end
    end
    return AnalysisConfig(; kwargs...)
end

# =============================================================================
# AnalysisResult - Return type from analyze()
# =============================================================================

"""
    AnalysisResult

Result container from `analyze()`. Contains the final SMLD plus intermediate results.

# Fields
- `smld`: Final BasicSMLD (after all processing)
- `smld_raw`: Raw fitted SMLD (before filtering/frameconnect)
- `roi_batch`: Detection results
- `timings`: Dict of step timings in seconds
- `workflow`: SMLMWorkflow provenance tracking
"""
struct AnalysisResult{T}
    smld::BasicSMLD{T}
    smld_raw::BasicSMLD{T}
    roi_batch::ROIBatch
    timings::Dict{String, Float64}
    workflow::SMLMWorkflow
end

# Convenience accessors
Base.getproperty(r::AnalysisResult, s::Symbol) =
    s == :emitters ? r.smld.emitters : getfield(r, s)

function Base.show(io::IO, r::AnalysisResult)
    n_raw = length(r.smld_raw.emitters)
    n_final = length(r.smld.emitters)
    total_time = sum(values(r.timings))
    print(io, "AnalysisResult: $n_final localizations ($n_raw raw) in $(round(total_time, digits=2))s")
end
