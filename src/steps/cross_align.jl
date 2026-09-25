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
  `AlignConfig()`, i.e. upstream's defaults)

# Example
```julia
CrossAlignConfig()                                        # entropy (CC + entropy refinement)
CrossAlignConfig(align=AlignConfig(method=:fft))          # CC only
```
"""
@kwdef struct CrossAlignConfig <: AbstractMultiTargetStep
    align::SMLMDriftCorrection.AlignConfig =
        SMLMDriftCorrection.AlignConfig()
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
        _align_edge_geometry!(aligned_smlds, smlds, align_info)
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

"""
    _align_edge_geometry!(aligned, smlds, info) -> aligned

`SMLMDriftCorrection.align_smld` only moves emitters — it never touches
`smld.metadata`. Edge classification (`steps/edgeclassify.jl`) stores cell-mask
geometry there (`"edge_outer_polygon"`, `"edge_cells"`) in the same μm frame as
the emitters, so saving a post-alignment SMLD would otherwise pair aligned
emitters with an unaligned mask. For every non-reference channel whose
metadata carries either key, replace it with a copy run through the exact
point transform `align_smld` applied to that channel's emitters — recovered
from the original (`smlds[i]`) and aligned (`aligned[i]`) emitter positions
themselves, not replayed from `info`'s diagnostic coefficients (for `:affine`,
those describe two sequential passes, and summing them is not the same as
composing them). `aligned[1]` is the caller's own object (`smlds[1]`,
untouched by `align_smld`) and is never mutated here — its correction is the
identity anyway.
"""
function _align_edge_geometry!(aligned::Vector{<:SMLMData.BasicSMLD},
                                smlds::Vector{<:SMLMData.BasicSMLD},
                                info::SMLMDriftCorrection.AlignInfo)
    for i in 2:length(aligned)
        md = aligned[i].metadata
        (haskey(md, "edge_outer_polygon") || haskey(md, "edge_cells")) || continue
        tf = _edge_point_transform(smlds[i], aligned[i], info, i)
        if tf === nothing
            @warn "Cross-align: dropping edge geometry for channel $i -- the aligned emitters could not be reproduced as an exact map of the original emitters, so the mask cannot be aligned with them" i
            delete!(md, "edge_outer_polygon")
            delete!(md, "edge_cells")
            continue
        end
        if haskey(md, "edge_outer_polygon")
            md["edge_outer_polygon"] = tf.(md["edge_outer_polygon"])
        end
        if haskey(md, "edge_cells")
            md["edge_cells"] = SMLMClustering.CellPolygon[
                SMLMClustering.CellPolygon(tf.(c.outer),
                    Vector{NTuple{2,Float64}}[tf.(h) for h in c.holes])
                for c in md["edge_cells"]
            ]
        end
    end
    aligned
end

# Point transform matching the emitter correction `align_smld` applied to
# channel `i`, exactly. For `:shift` this is just the recorded translation
# (see `correctdrift!` in SMLMDriftCorrection `intrainter.jl`). For `:affine`,
# `align_smld` composes two sequential affine passes (global shift, then
# affine — twice; see `_align_affine_fft` in `align.jl`) but its `diagnostic`
# only records each pass's own coefficients summed together, which is NOT the
# same map as composing the two passes. Rather than replay that (inexact)
# sum, recover the exact composed affine map by least-squares fit against the
# actual before/after emitter positions -- `aligned[i]` is a deepcopy of
# `smld` with the same emitter order, so this fit is exact up to floating-point
# roundoff whenever the true correction is affine. Returns `nothing` if it
# isn't (mismatched emitter counts, or a residual too large to be roundoff).
function _edge_point_transform(smld::SMLMData.BasicSMLD, aligned_smld::SMLMData.BasicSMLD,
                                info::SMLMDriftCorrection.AlignInfo, i::Int)
    if info.transform == :shift
        dx, dy = info.shifts[i][1], info.shifts[i][2]
        return p -> (p[1] - dx, p[2] - dy)
    else  # :affine
        n = length(smld.emitters)
        n == length(aligned_smld.emitters) || return nothing
        X = Matrix{Float64}(undef, n, 3)
        T = Matrix{Float64}(undef, n, 2)
        for (k, (e0, e1)) in enumerate(zip(smld.emitters, aligned_smld.emitters))
            X[k, 1] = e0.x
            X[k, 2] = e0.y
            X[k, 3] = 1.0
            T[k, 1] = e1.x
            T[k, 2] = e1.y
        end
        M = X \ T   # least-squares affine map: [x y 1] * M ≈ [x' y']
        maximum(abs, X * M .- T) <= 1e-6 || return nothing
        return p -> (p[1] * M[1, 1] + p[2] * M[2, 1] + M[3, 1],
                     p[1] * M[1, 2] + p[2] * M[2, 2] + M[3, 2])
    end
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
