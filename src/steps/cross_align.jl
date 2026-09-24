"""
Cross-channel alignment step — aligns multiple SMLDs to a common reference.

Dispatches on `CrossAlignConfig <: AbstractMultiTargetStep` operating on
`Vector{BasicSMLD}`. State-modifying: returns aligned SMLDs.
"""

"""
    CrossAlignConfig <: AbstractMultiTargetStep

Configuration for cross-channel alignment in the multi-target pipeline.

Wraps `SMLMDriftCorrection.align_smld`, which uses entropy-based or FFT
cross-correlation alignment. Every alignment parameter lives on the upstream
`AlignConfig` (re-exported), which is passed through unchanged — including its
`verbose` field.

# Fields
- `align::AlignConfig`: upstream alignment config (default:
  `AlignConfig(method=:entropy, maxn=100, histbinsize=0.05)`)

# Example
```julia
CrossAlignConfig()                                        # entropy (CC + entropy refinement)
CrossAlignConfig(align=AlignConfig(method=:fft))          # CC only
```
"""
@kwdef struct CrossAlignConfig <: AbstractMultiTargetStep
    align::SMLMDriftCorrection.AlignConfig =
        SMLMDriftCorrection.AlignConfig(method=:entropy, maxn=100, histbinsize=0.05)
end

step_name(::CrossAlignConfig) = "crossalign"

"""
    crossalign_step(smlds, cfg; outdir, step_number, verbose) -> (aligned_smlds, CrossAlignInfo)

Align multiple SMLDs to a common reference using entropy-based alignment.
"""
function crossalign_step(smlds::Vector{<:SMLMData.BasicSMLD}, cfg::CrossAlignConfig;
                         outdir::Union{String,Nothing}=nothing,
                         step_number::Int=0,
                         verbose::Int=Verbosity.STANDARD)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)

    v >= Verbosity.PROGRESS && @info "[$step_number] crossalign: $(cfg.align.method), $(length(smlds)) channels"

    local aligned_smlds, align_info
    t = @elapsed begin
        (aligned_smlds, align_info) = SMLMDriftCorrection.align_smld(smlds, cfg.align)
    end

    # Convert shifts to nm and compute max
    shifts_nm = [s .* 1000 for s in align_info.shifts]
    max_shift_nm = maximum(sqrt(sum(s .^ 2)) for s in shifts_nm)

    if dir !== nothing
        mkpath(dir)
        _save_config!(dir, cfg)
        _save_info!(dir, align_info)
        if v >= Verbosity.STANDARD
            _write_crossalign_stats(dir, cfg, align_info, shifts_nm, max_shift_nm, t)
        end
    end

    v >= Verbosity.PROGRESS && @info "  -> aligned $(length(smlds)) channels, max shift $(round(max_shift_nm, digits=1))nm ($(round(t, digits=2))s)"

    info = CrossAlignInfo(align_info, align_info.shifts, max_shift_nm, t)
    (aligned_smlds, info)
end

_step_summary(info::CrossAlignInfo) = Dict{Symbol,Any}(
    :max_shift_nm => round(info.max_shift_nm, digits=1),
    :n_channels => length(info.shifts),
    :method => info.align_info.method,
)

"""
    analyze(smlds::Vector{BasicSMLD}, cfg::CrossAlignConfig; kwargs...) -> (aligned_smlds, StepInfo)

Multi-target dispatch: cross-channel alignment. Modifies SMLDs.
"""
function analyze(smlds::Vector{<:SMLMData.BasicSMLD}, cfg::CrossAlignConfig;
                 outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD, kwargs...)
    t = @elapsed (aligned, ca_info) = crossalign_step(smlds, cfg;
        outdir=outdir, step_number=step_number, verbose=verbose)
    (aligned, StepInfo(step_number, cfg, t, _step_summary(ca_info); info=ca_info))
end

function _write_crossalign_stats(dir, cfg::CrossAlignConfig, align_info, shifts_nm, max_shift_nm, t)
    filepath = joinpath(dir, "stats.md")
    open(filepath, "w") do io
        println(io, "# Cross-Channel Alignment Statistics\n")
        println(io, "## Summary")
        println(io, "- **Method**: $(cfg.align.method)")
        println(io, "- **Channels**: $(length(align_info.shifts))")
        println(io, "- **Max shift**: $(round(max_shift_nm, digits=1)) nm")
        println(io, "- **Time**: $(round(t, digits=2))s")
        println(io)
        println(io, "## Per-Channel Shifts")
        println(io, "| Channel | X (nm) | Y (nm) | Magnitude (nm) |")
        println(io, "|---------|--------|--------|----------------|")
        for (i, s) in enumerate(shifts_nm)
            mag = sqrt(sum(s .^ 2))
            if length(s) >= 2
                println(io, "| $i | $(round(s[1], digits=1)) | $(round(s[2], digits=1)) | $(round(mag, digits=1)) |")
            else
                println(io, "| $i | $(round(s[1], digits=1)) | - | $(round(mag, digits=1)) |")
            end
        end
    end
end
