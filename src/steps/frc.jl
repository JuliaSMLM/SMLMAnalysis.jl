# NOT YET WIRED IN. This file is deliberately absent from the include list in
# src/SMLMAnalysis.jl, so it is inert source and cannot affect loading or resolution.
#
# SMLMResolution is not yet registered in General. Adding it to [deps] before it is
# registered AND its 3-day new-package hold has cleared would fail AutoMerge for
# SMLMAnalysis itself, because AutoMerge resolves against the live registry.
#
# To land, once SMLMResolution is live:
#   1. Project.toml [deps]    SMLMResolution = "<uuid>"
#   2. Project.toml [compat]  SMLMResolution = "0.1"
#   3. src/SMLMAnalysis.jl    include("steps/frc.jl")  (next to the other steps)
#   4. src/SMLMAnalysis.jl    export FRCConfig, FRCInfo
#   5. docs/src/steps/frc.md  step page + entry in docs/src/steps/index.md
#   6. Verify: fast + SMLM_TEST_FULL tiers, and that FRC composes as a diagnostic
#      (the SMLD must come out of the pipeline unchanged).
#
# Written against SMLMResolution's final FRCConfig. Conventions follow
# src/steps/clustering.jl (wrapping an upstream-owned config) and src/steps/render.jl
# (pass-through diagnostic return).

# Upstream owns the config and info types; we alias and dispatch, never re-declare.
const FRCConfig = SMLMResolution.FRCConfig
const FRCInfo   = SMLMResolution.FRCInfo

# step_name is intentionally not defined: the generic fallback in src/types.jl strips
# "Config" and lowercases, giving FRCConfig -> "frc" and an output dir of NN_frc/.

"""
    analyze(smld, cfg::FRCConfig; kwargs...) -> (smld, StepInfo)

Estimate image resolution by Fourier Ring Correlation (Nieuwenhuizen et al. 2013,
doi:10.1038/nmeth.2488) via SMLMResolution.

**The SMLD is returned unchanged.** FRC is a diagnostic, like `RenderConfig` and the
clustering statistics steps; the resolution and its diagnostics live in the step's
`FRCInfo`, reachable as `stepinfo(info, :frc).info`.

Spurious-correlation (`Q`) correction is off by default. The available `:lowq`
estimator is documented upstream as biased in the structure-dominated regime — on
simulated data with a true `Q` of zero it returns several, inflating a ~46 nm
resolution to several microns. Pass a known `Q::Float64` if you have one measured
independently.
"""
function analyze(smld::BasicSMLD, cfg::FRCConfig;
                 outdir=nothing, step_number::Int=0,
                 verbose::Int=Verbosity.STANDARD, kwargs...)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)
    v >= Verbosity.PROGRESS &&
        @info "[$step_number] $(step_name(cfg))" n_locs=length(smld.emitters)

    t = @elapsed (analysis = SMLMResolution.compute_frc_resolution(smld, cfg))
    info = SMLMResolution.FRCInfo(analysis; elapsed_s = t)

    # frc_summary omits keys whose value was not computed, so every key present is
    # numeric. Required: _toml_value stringifies `nothing` to "nothing", which would
    # otherwise flip a provenance key's type between runs.
    summary = SMLMResolution.frc_summary(analysis)

    if dir !== nothing && v >= Verbosity.STANDARD
        mkpath(dir)
        _save_config!(dir, cfg)
        _save_info!(dir, info)
    end

    if v >= Verbosity.PROGRESS
        corr = get(summary, :resolution_corrected_nm, nothing)
        msg = corr === nothing ? "" : ", corrected $(round(corr, digits=1)) nm"
        @info "  → resolution $(round(summary[:resolution_nm], digits=1)) nm$msg"
    end

    # Pass-through: the SMLD continues down the pipeline untouched.
    (smld, StepInfo(step_number, cfg, t, summary; info=info))
end
