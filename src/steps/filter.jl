"""
Filtering step - filters localizations by various criteria
"""

"""
    FilterConfig <: AbstractSMLMConfig

Quality-based filtering of localizations. All criteria use `(min, max)` tuples.

# Keywords
- `photons`: Photon count range, e.g. `(500.0, Inf)`
- `precision`: Localization precision range in microns, e.g. `(0.0, 0.007)`
- `pvalue`: Goodness-of-fit p-value range, e.g. `(1e-3, 1.0)`
- `psf_sigma`: PSF width filter. `:auto` uses mode ± 10%, or explicit `(min, max)` in microns

All filters default to `nothing` (disabled).
"""
@kwdef struct FilterConfig <: SMLMData.AbstractSMLMConfig
    # All filters use (min, max) tuples. Use -Inf/Inf for unbounded.
    photons::Union{Tuple{Float64, Float64}, Nothing} = nothing      # (min, max)
    precision::Union{Tuple{Float64, Float64}, Nothing} = nothing    # (min, max) in microns
    pvalue::Union{Tuple{Float64, Float64}, Nothing} = nothing       # (min, max)
    # PSF sigma: :auto (mode ± 10%), or (min, max) tuple in microns
    psf_sigma::Union{Symbol, Tuple{Float64, Float64}, Nothing} = nothing
end

"""
    filter_step(smld, cfg; smld_raw=nothing, outdir=nothing, step_number=0, verbose=Verbosity.STANDARD)

Filter localizations by quality criteria. Returns `(filtered_smld, info)`.

# Arguments
- `smld::BasicSMLD`: Input localizations
- `cfg::FilterConfig`: Filter criteria

# Keyword Arguments
- `smld_raw`: Original unfiltered SMLD for detailed output diagnostics
- `outdir`: Output directory (nothing to skip file output)
- `step_number`: Step number for output directory naming
- `verbose`: Verbosity level

# Returns
`(filtered_smld, (step_record, n_before, n_after))`
"""
function filter_step(smld::BasicSMLD, cfg::FilterConfig;
                     smld_raw::Union{BasicSMLD,Nothing}=nothing,
                     outdir::Union{String,Nothing}=nothing,
                     step_number::Int=0,
                     verbose::Int=Verbosity.STANDARD)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)

    v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg))" photons=cfg.photons precision=cfg.precision

    n_before = length(smld.emitters)
    t = @elapsed filtered = _filter_smld(smld, cfg)
    n_after = length(filtered.emitters)

    summary = Dict{Symbol,Any}(
        :n_before => n_before,
        :n_after => n_after,
        :acceptance => round(n_after / n_before, digits=3)
    )
    record = StepRecord(step_number, cfg, t, summary)

    if dir !== nothing
        _save_filter_outputs!(dir, cfg, v, t, n_before, n_after, smld_raw, filtered)
    end

    v >= Verbosity.PROGRESS && @info "  → $n_after / $n_before ($(round(t, digits=2))s)"
    (filtered, (step_record=record, n_before=n_before, n_after=n_after))
end

function _filter_smld(smld::BasicSMLD, cfg::FilterConfig)
    emitters = smld.emitters
    mask = trues(length(emitters))

    if cfg.photons !== nothing
        lo, hi = cfg.photons
        mask .&= [lo <= e.photons <= hi for e in emitters]
    end

    if cfg.precision !== nothing
        lo, hi = cfg.precision
        mask .&= [lo <= max(e.σ_x, e.σ_y) <= hi for e in emitters]
    end

    if cfg.pvalue !== nothing
        lo, hi = cfg.pvalue
        mask .&= [lo <= e.pvalue <= hi for e in emitters]
    end

    if cfg.psf_sigma !== nothing && length(emitters) > 0
        # Determine bounds: :auto calculates mode ± 10%, or use explicit (min, max)
        if hasproperty(emitters[1], :σ)
            lo, hi = _get_psf_sigma_bounds(cfg.psf_sigma, [e.σ for e in emitters])
            if lo > 0 && hi > 0
                mask .&= [lo <= e.σ <= hi for e in emitters]
            end
        elseif hasproperty(emitters[1], :σx) && hasproperty(emitters[1], :σy)
            lo_x, hi_x = _get_psf_sigma_bounds(cfg.psf_sigma, [e.σx for e in emitters])
            lo_y, hi_y = _get_psf_sigma_bounds(cfg.psf_sigma, [e.σy for e in emitters])
            if lo_x > 0 && hi_x > 0 && lo_y > 0 && hi_y > 0
                mask .&= [lo_x <= e.σx <= hi_x && lo_y <= e.σy <= hi_y for e in emitters]
            end
        end
    end

    filtered = emitters[mask]
    BasicSMLD(filtered, smld.camera, smld.n_frames, smld.n_datasets, smld.metadata)
end

"""
    _get_psf_sigma_bounds(range_spec, values) -> (lo, hi)

Calculate PSF sigma filter bounds.
- `:auto` → mode ± 10%
- `(min, max)` → explicit bounds in microns
"""
function _get_psf_sigma_bounds(range_spec, values::Vector)
    if range_spec === :auto
        mode = _calculate_mode(values)
        mode > 0 || return (0.0, 0.0)
        return (mode * 0.90, mode * 1.10)
    elseif range_spec isa Tuple{Float64, Float64}
        return range_spec
    else
        error("psf_sigma_range must be :auto or (min, max) tuple, got: $range_spec")
    end
end

function _save_filter_outputs!(dir::String, cfg::FilterConfig, v::Int, t::Float64,
                               n_before::Int, n_after::Int,
                               smld_raw::Union{BasicSMLD,Nothing}, smld_filtered::BasicSMLD)
    mkpath(dir)
    _save_config!(dir, cfg)

    if v >= Verbosity.STANDARD
        _write_filter_stats(dir, cfg, n_before, n_after, t)
    end

    if v >= Verbosity.DETAILED && smld_raw !== nothing
        _save_filter_detailed(dir, smld_raw, smld_filtered, cfg)
    end
end

function _write_filter_stats(dir, cfg, n_before, n_after, t)
    acceptance = n_after / n_before

    filepath = joinpath(dir, "stats.md")
    open(filepath, "w") do io
        println(io, "# Filter Statistics\n")
        println(io, "## Summary")
        println(io, "- **Input**: $n_before")
        println(io, "- **Output**: $n_after")
        println(io, "- **Acceptance**: $(round(100*acceptance, digits=1))%")
        println(io, "- **Time**: $(round(t*1000, digits=1))ms")
        println(io, "")
        println(io, "## Criteria Applied")
        if cfg.photons !== nothing
            lo, hi = cfg.photons
            println(io, "- photons: $(lo) - $(hi == Inf ? "∞" : hi)")
        end
        if cfg.precision !== nothing
            lo, hi = cfg.precision
            println(io, "- precision: $(round(lo*1000, digits=1)) - $(hi == Inf ? "∞" : round(hi*1000, digits=1)) nm")
        end
        if cfg.pvalue !== nothing
            lo, hi = cfg.pvalue
            println(io, "- pvalue: $(lo) - $(hi)")
        end
        if cfg.psf_sigma !== nothing
            if cfg.psf_sigma === :auto
                println(io, "- psf_sigma: :auto (mode ± 10%)")
            else
                lo, hi = cfg.psf_sigma
                println(io, "- psf_sigma: $(round(lo*1000, digits=1)) - $(round(hi*1000, digits=1)) nm")
            end
        end
    end
end

function _save_filter_detailed(dir, smld_raw, smld_filtered, cfg)
    # Show which criteria rejected what
    emitters = smld_raw.emitters
    n = length(emitters)

    filepath = joinpath(dir, "detailed_stats.md")
    open(filepath, "w") do io
        println(io, "# Filter Breakdown\n")
        println(io, "| Criterion | Pass | Fail | % Pass |")
        println(io, "|-----------|------|------|--------|")

        if cfg.photons !== nothing
            lo, hi = cfg.photons
            pass = sum(lo <= e.photons <= hi for e in emitters)
            hi_str = hi == Inf ? "∞" : string(hi)
            println(io, "| Photons ∈ [$lo, $hi_str] | $pass | $(n - pass) | $(round(100*pass/n, digits=1))% |")
        end

        if cfg.precision !== nothing
            lo, hi = cfg.precision
            pass = sum(lo <= max(e.σ_x, e.σ_y) <= hi for e in emitters)
            hi_str = hi == Inf ? "∞" : "$(round(hi*1000, digits=1))nm"
            println(io, "| Precision ∈ [$(round(lo*1000, digits=1))nm, $hi_str] | $pass | $(n - pass) | $(round(100*pass/n, digits=1))% |")
        end

        if cfg.pvalue !== nothing
            lo, hi = cfg.pvalue
            pass = sum(lo <= e.pvalue <= hi for e in emitters)
            println(io, "| P-value ∈ [$lo, $hi] | $pass | $(n - pass) | $(round(100*pass/n, digits=1))% |")
        end
    end
end
