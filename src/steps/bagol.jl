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
    bagol_smld, diagnostics = SMLMBaGoL.run_bagol(smld, cfg_run)

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

        # BaGoL-specific renders: partition circles + overlay
        _render_bagol_diagnostics(smld, bagol_smld, diagnostics, dir, cfg)
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

"""Copy of `smld` with per-loc σ inflated in quadrature by `se_adjust` (scalar τ or (τx,τy), µm)
so rendered ellipses show the σ BaGoL used for grouping (matches render_report). `nothing`→unchanged."""
function _se_inflated_smld(smld::BasicSMLD, se_adjust)
    se_adjust === nothing && return smld
    se_adjust isa Symbol && return smld   # :auto/unresolved sentinel (SMLMBaGoL v0.3.5) — nothing numeric to inflate
    τx, τy = se_adjust isa Number ? (se_adjust, se_adjust) : (se_adjust[1], se_adjust[2])
    (τx <= 0 && τy <= 0) && return smld
    out = deepcopy(smld)
    for e in out.emitters
        hasproperty(e, :σ_x) && (e.σ_x = sqrt(e.σ_x^2 + τx^2))
        hasproperty(e, :σ_y) && (e.σ_y = sqrt(e.σ_y^2 + τy^2))
    end
    out
end

"""
Render BaGoL-specific diagnostic images: partition ellipses and loc/emitter overlay.

Matches SMLMBaGoL.render_report naming: `circles.png` (locs + MAP-N overlay)
and `partitions.png` (partition-colored localizations).

These renders need the pre-BaGoL localizations, so they are produced inside the BaGoL
step rather than as separate pipeline steps.
"""
function _render_bagol_diagnostics(smld::BasicSMLD, bagol_smld::BasicSMLD,
                                   diagnostics::SMLMBaGoL.BaGoLDiagnostics,
                                   dir::String, cfg::BaGoLConfig)
    zoom = 50.0
    partition_ids = diagnostics.partition_ids

    # Shared target so both renders have identical bounds
    target = SMLMRender.create_target_from_smld(smld; zoom=zoom)

    # Circles overlay: white localizations (σ inflated by the APPLIED se_adjust = the grouping σ) + red MAP-N emitters.
    # Use diagnostics.se_adjust (the resolved (τx,τy) BaGoL actually applied), NOT cfg.se_adjust — under
    # finder-by-default (SMLMBaGoL v0.3.5, se_adjust=:auto) cfg.se_adjust is the :auto Symbol, not a number.
    try
        bg_smld = _se_inflated_smld(smld, hasproperty(diagnostics, :se_adjust) ? diagnostics.se_adjust : nothing)
        (bg_img, _) = SMLMRender.render(bg_smld; strategy=EllipseRender(),
            color=:white, target=target, clip_percentile=nothing)
        (fg_img, _) = SMLMRender.render(bagol_smld; strategy=EllipseRender(),
            color=:red, target=target, clip_percentile=nothing)
        combined = SMLMRender.compose(bg_img, fg_img; blend=:replace)
        SMLMRender.save_image(joinpath(dir, "circles.png"), combined)
    catch e
        @warn "BaGoL circles render failed" exception=e
    end

    # Partition-colored localizations (copy emitters so the shared input smld is never mutated;
    # we repurpose the dataset field as a partition label only for this render).
    if length(partition_ids) == length(smld.emitters)
        try
            part_emitters = deepcopy(smld.emitters)
            for (i, e) in enumerate(part_emitters)
                e.dataset = partition_ids[i]
            end
            part_smld = BasicSMLD(part_emitters, smld.camera, smld.n_frames, 1)
            SMLMRender.render(part_smld; strategy=EllipseRender(),
                color_by=:dataset, categorical=true, zoom=zoom,
                filename=joinpath(dir, "partitions.png"))
        catch e
            @warn "BaGoL partitions render failed" exception=e
        end
    end
end

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
