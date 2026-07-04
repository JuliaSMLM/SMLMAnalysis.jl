"""
BaGoL (Bayesian Grouping of Localizations) step for SMLMAnalysis.

Groups localizations into emitters using Bayesian inference via RJMCMC.
Requires frame-connected (and ideally calibrated) data.
"""

# BaGoLConfig is defined in SMLMBaGoL (upstream owns config, like DriftConfig/RenderConfig)
const BaGoLConfig = SMLMBaGoL.BaGoLConfig

"""
    bagol_step(smld, cfg; outdir=nothing, step_number=0, verbose=Verbosity.STANDARD)

Group localizations into emitters via BaGoL. Returns `(bagol_smld, BaGoLInfo)`.
"""
function bagol_step(smld::BasicSMLD, cfg::BaGoLConfig;
                    outdir::Union{String,Nothing}=nothing,
                    step_number::Int=0,
                    verbose::Int=Verbosity.STANDARD,
                    checkpoint::Int=Checkpoint.EXPENSIVE)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)

    v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg))" μ=cfg.μ shape=cfg.shape n_iterations=cfg.n_iterations

    n_locs_in = length(smld.emitters)

    # Config dispatch: BaGoLConfig fields map 1:1 to run_bagol kwargs
    # Override verbose from pipeline verbosity level
    cfg_run = BaGoLConfig(; (f => getfield(cfg, f) for f in fieldnames(BaGoLConfig))...,
                            verbose = v >= Verbosity.PROGRESS)
    # keep_se_finder (runtime kwarg, SMLMBaGoL v0.3.7-DEV): stash the full estimate_se_adjust
    # result on diagnostics.se_finder at render verbosity so the finder plot uses one finder run.
    bagol_smld, diagnostics = SMLMBaGoL.run_bagol(smld, cfg_run; keep_se_finder = v >= Verbosity.STANDARD)

    n_emitters = diagnostics.n_emitters
    compression = n_locs_in > 0 ? round(n_locs_in / max(1, n_emitters), digits=1) : 0.0

    # τ̂ the finder actually applied: diagnostics.se_adjust is the applied (τx,τy) μm (or nothing under :auto skip)
    tau_um = diagnostics.se_adjust === nothing ? 0.0 : Float64(diagnostics.se_adjust[1])
    info = BaGoLInfo(n_locs_in, n_emitters, compression,
                     diagnostics.final_μ, diagnostics.final_shape,
                     diagnostics.n_partitions, tau_um, diagnostics)

    # Diagnostic outputs via SMLMBaGoL report system
    if dir !== nothing && v >= Verbosity.STANDARD
        mkpath(dir)
        _save_config!(dir, cfg)
        _save_info!(dir, diagnostics)
        # compute_report wants locs WITHOUT FC link IDs (track_id is not a GT emitter
        # assignment). Operate on a copy so the shared pipeline smld is never mutated in place
        # (downstream steps / the caller may still read the original track_id).
        locs_for_report = deepcopy(smld)
        for e in locs_for_report.emitters
            e.track_id = 0
        end
        report = SMLMBaGoL.compute_report(bagol_smld, diagnostics; locs_smld=locs_for_report)
        SMLMBaGoL.write_report(report; output_dir=dir)
        SMLMBaGoL.plot_report(report; output_dir=dir)

        # BaGoL renders via upstream render_report (SMLMBaGoL v0.3.7-DEV): writes the render_<noun>
        # set (render_mapn / render_sr / render_circles / render_partitions) on ONE shared target;
        # se_adjust inflation handled inside render_report.
        se = diagnostics.se_adjust === nothing ? 0.0 : diagnostics.se_adjust
        SMLMBaGoL.render_report(smld, bagol_smld; output_dir=dir,
            partition_ids=diagnostics.partition_ids, se_adjust=se, zoom=50, prefix="render")
        # finder plot at the same (STANDARD) verbosity as the renders (one finder run via keep_se_finder)
        diagnostics.se_finder !== nothing && SMLMBaGoL.plot_se_adjust(diagnostics.se_finder; output_dir=dir)
    end

    if dir !== nothing && checkpoint >= Checkpoint.EXPENSIVE
        _save_step_smld(dir, bagol_smld; filename="smld_bagol.jld2")
    end

    v >= Verbosity.PROGRESS && @info "  → $n_emitters emitters from $n_locs_in locs ($(compression)x compression, $(round(diagnostics.final_μ, digits=1)) locs/emitter)"

    (bagol_smld, info)
end

_step_summary(info::BaGoLInfo) = Dict{Symbol,Any}(
    :n_locs_in => info.n_locs_in,
    :n_emitters => info.n_emitters,
    :compression => info.compression,
    :final_μ => round(info.final_μ, digits=2),
    :final_shape => round(info.final_shape, digits=2),
    :n_partitions => info.n_partitions,
    :tau_um => round(info.tau_um, digits=4),  # applied se_adjust τ̂ (μm) — finder-by-default (SMLMBaGoL v0.3.5)
    :se_adjust => (hasproperty(info.diagnostics, :se_adjust) ? info.diagnostics.se_adjust : nothing),  # applied (τx,τy) provenance (SMLMBaGoL v0.3.1-DEV)
)

# (BaGoL diagnostic renders + se_adjust inflation are produced by SMLMBaGoL.render_report,
#  called from bagol_step above — no local duplicate needed since SMLMBaGoL v0.3.7-DEV.)

"""
    analyze(smld, cfg::BaGoLConfig; kwargs...) -> (bagol_smld, StepInfo)

Group localizations into emitters via Bayesian inference (BaGoL).
"""
function analyze(smld::BasicSMLD, cfg::BaGoLConfig;
                 outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD,
                 checkpoint::Int=Checkpoint.EXPENSIVE, kwargs...)
    t = @elapsed (bagol_smld, bagol_info) = bagol_step(smld, cfg;
        outdir=outdir, step_number=step_number, verbose=verbose, checkpoint=checkpoint)
    (bagol_smld, StepInfo(step_number, cfg, t, _step_summary(bagol_info); info=bagol_info))
end
