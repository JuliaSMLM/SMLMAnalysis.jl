"""
    SMLMAnalysisPSFLearningExt

Package extension that activates the `psflearning` analyze() step when
`PSFLearning` is loaded alongside `SMLMAnalysis`.

Kept as an extension (PSFLearning is a weakdep) so SMLMAnalysis's core dependency
tree stays free of the Reactant/Enzyme deep-learning stack.

`psflearning` is a calibration-phase step: it consumes bead calibration data and
produces a MicroscopePSFs PSF (ScalarPSF/VectorPSF), NOT an SMLD — so it does not
belong in the linear smld-threading `steps` vector. It writes `psf.h5`
(`save_psf`, io_version 1.0.0) which the deepfit trainer loads directly via
`MicroscopePSFs.load_psf` (format locked with @smlmdeepfit — zero converter).
"""
module SMLMAnalysisPSFLearningExt

using SMLMData
using PSFLearning
import SMLMAnalysis: analyze, step_name, _produces_smld,
    step_outdir, _save_config!, StepInfo, Verbosity

step_name(::PSFLearning.PSFLearningConfig) = "psflearning"
_produces_smld(::PSFLearning.PSFLearningConfig) = false

"""Defensive numeric-field summary for an upstream info struct."""
function _numeric_summary(info)
    d = Dict{Symbol,Any}()
    for f in fieldnames(typeof(info))
        v = getfield(info, f)
        v isa Number && (d[f] = v)
    end
    d
end

_as_info(x) = x isa SMLMData.AbstractSMLMInfo ? x : nothing

"""
    analyze(stack, cfg::PSFLearningConfig; z_positions) -> (psf, StepInfo)

Calibration-phase PSF learning. Wraps the detection-free core
`learn_psf(z_positions, cfg; data=stack)` (the clean boundary, since SMLMAnalysis
owns detection upstream — do not re-detect here).

# Arguments
- `stack::AbstractArray{<:Real,3}`: pre-extracted bead ROI z-stack, shape
  `(N_z, M, M)` with `M == cfg.roi_size`.
- `z_positions::AbstractVector`: axial position (µm) of each z-plane, length `N_z`.

Returns the learned MicroscopePSFs PSF (`ScalarPSF`/`VectorPSF`) and writes
`psf.h5` (canonical MicroscopePSFs HDF5, io_version 1.0.0) into the step dir for
the deepfit trainer to load via `MicroscopePSFs.load_psf`.
"""
function analyze(stack::AbstractArray{<:Real,3}, cfg::PSFLearning.PSFLearningConfig;
        z_positions::AbstractVector, outdir=nothing, step_number::Int=0,
        verbose::Int=Verbosity.STANDARD, kwargs...)
    dir = step_outdir(outdir, step_number, cfg)
    verbose >= Verbosity.PROGRESS && @info "[$step_number] psflearning"
    zp = Vector{Float32}(z_positions)
    data = Float32.(stack)
    t = @elapsed ((psf, info) = PSFLearning.learn_psf(zp, cfg; data=data))
    if dir !== nothing
        mkpath(dir)
        _save_config!(dir, cfg)
        try
            PSFLearning.save_psf(joinpath(dir, "psf.h5"), info)
        catch err
            @warn "psflearning: save_psf failed" err
        end
    end
    verbose >= Verbosity.PROGRESS && @info "  → PSF learned ($(round(t, digits=1))s)"
    (psf, StepInfo(step_number, cfg, t, _numeric_summary(info); info=_as_info(info)))
end

end # module
