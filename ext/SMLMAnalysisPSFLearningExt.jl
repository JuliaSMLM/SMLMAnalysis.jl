"""
    SMLMAnalysisPSFLearningExt

Package extension that adds the `psflearning` analyze() step when `PSFLearning`
is loaded alongside `SMLMAnalysis`.

PSFLearning is a weak dependency, so SMLMAnalysis's own dependency tree stays
free of its Reactant/Enzyme stack.

`psflearning` is a calibration step: it takes a bead z-stack and returns a
MicroscopePSFs PSF, not an SMLD, so it does not go in the `steps` vector. It
writes `psf.h5` with `PSFLearning.save_psf`. That file is the interface to the
fitters: `MicroscopePSFs.load_psf` reads it without PSFLearning loaded.
"""
module SMLMAnalysisPSFLearningExt

# SMLMAnalysis must load before PSFLearning. On a CUDA machine SMLMAnalysis loads
# CUDNN_jll (through SMLMBoxer), and CUDNN_jll fails to initialise once Reactant's
# bundled cuDNN is already loaded ("undefined symbol ... libcudnn_graph.so.9");
# the other order works.
using SMLMAnalysis: SMLMAnalysis, StepInfo, Verbosity, _save_config!, step_outdir
import SMLMAnalysis: _produces_smld, analyze, step_name
using SMLMData: SMLMData
using PSFLearning: PSFLearning

step_name(::PSFLearning.PSFLearningConfig) = "psflearning"
_produces_smld(::PSFLearning.PSFLearningConfig) = false

# The numeric fields of an upstream info struct, for StepInfo's summary.
function _numeric_summary(info)
    d = Dict{Symbol, Any}()
    for f in fieldnames(typeof(info))
        v = getfield(info, f)
        v isa Number && (d[f] = v)
    end
    return d
end

_as_info(x) = x isa SMLMData.AbstractSMLMInfo ? x : nothing

"""
    analyze(stack, cfg::PSFLearningConfig; z_positions) -> (psf, StepInfo)

Learn a PSF from a bead z-stack with `PSFLearning.learn_psf(z_positions, cfg; data = stack)`.
SMLMAnalysis does the detection upstream, so this step does not detect beads.

# Arguments
- `stack::AbstractArray{<:Real,3}`: bead ROI z-stack of shape `(N_z, M, M)`, with
  `M == cfg.roi_size`.
- `z_positions::AbstractVector`: axial position (µm) of each plane, length `N_z`.

Returns the learned PSF and writes `psf.h5` (MicroscopePSFs HDF5) into the step
directory when `outdir` is set.
"""
function analyze(
        stack::AbstractArray{<:Real, 3}, cfg::PSFLearning.PSFLearningConfig;
        z_positions::AbstractVector, outdir = nothing, step_number::Int = 0,
        verbose::Int = Verbosity.STANDARD, kwargs...
    )
    dir = step_outdir(outdir, step_number, cfg)
    verbose >= Verbosity.PROGRESS && @info "[$step_number] psflearning"
    zp = Vector{Float32}(z_positions)
    data = Float32.(stack)
    t = @elapsed ((psf, info) = PSFLearning.learn_psf(zp, cfg; data = data))
    if dir !== nothing
        mkpath(dir)
        _save_config!(dir, cfg)
        try
            PSFLearning.save_psf(joinpath(dir, "psf.h5"), info)
        catch err
            @warn "psflearning: save_psf failed" err
        end
    end
    verbose >= Verbosity.PROGRESS && @info "  → PSF learned ($(round(t, digits = 1))s)"
    return (psf, StepInfo(step_number, cfg, t, _numeric_summary(info); info = _as_info(info)))
end

end # module
