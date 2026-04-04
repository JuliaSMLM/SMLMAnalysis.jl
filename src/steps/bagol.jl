"""
BaGoL (Bayesian Grouping of Localizations) step for SMLMAnalysis.

Groups localizations into emitters using Bayesian inference via RJMCMC.
Requires frame-connected (and ideally calibrated) data.
"""

"""
    BaGoLConfig <: AbstractSMLMConfig

Configuration for BaGoL grouping step. Fields map directly to `SMLMBaGoL.run_bagol` kwargs.

# Keywords
- `μ::Float64`: Mean localizations per emitter (default: 10.0)
- `shape::Float64`: NegBin shape parameter — 1.0 = exponential/dSTORM, >1 = peaked/DNA-PAINT (default: 2.0)
- `learn_distribution::Union{Bool,Symbol}`: Count distribution learning — `true` = learn both μ and shape,
  `false` = fix both, `:mu` = learn μ only, `:shape` = learn shape only (default: `true`)
- `n_iterations::Int`: Total MCMC iterations (default: 10000)
- `burn_in::Int`: Burn-in iterations before recording (default: 2000)
- `sync_interval::Int`: Iterations between global μ/shape updates (default: 500)
- `partition_sigma::Float64`: DBSCAN threshold in sigma units (default: 3.0)
- `min_partition_size::Int`: Minimum locs per partition; smaller dropped as noise (default: 0)
- `max_partition_size::Int`: Maximum locs per partition; larger are split (default: 1000)
- `skip_partition_size::Int`: Skip partitions larger than this (default: typemax(Int))
- `posterior_pixel_size::Float64`: Rao-Blackwellized posterior image pixel size in μm; 0.0 to disable (default: 0.002)
"""
@kwdef struct BaGoLConfig <: SMLMData.AbstractSMLMConfig
    # Count model
    μ::Float64 = 10.0
    shape::Float64 = 2.0
    learn_distribution::Union{Bool, Symbol} = true

    # MCMC
    n_iterations::Int = 10000
    burn_in::Int = 2000
    sync_interval::Int = 500

    # Partitioning
    partition_sigma::Float64 = 3.0
    min_partition_size::Int = 0
    max_partition_size::Int = 1000
    skip_partition_size::Int = typemax(Int)

    # Posterior image
    posterior_pixel_size::Float64 = 0.002
end

"""
    bagol_step(smld, cfg; outdir=nothing, step_number=0, verbose=Verbosity.STANDARD)

Group localizations into emitters via BaGoL. Returns `(bagol_smld, BaGoLInfo)`.
"""
function bagol_step(smld::BasicSMLD, cfg::BaGoLConfig;
                    outdir::Union{String,Nothing}=nothing,
                    step_number::Int=0,
                    verbose::Int=Verbosity.STANDARD)
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)

    v >= Verbosity.PROGRESS && @info "[$step_number] $(step_name(cfg))" μ=cfg.μ shape=cfg.shape n_iterations=cfg.n_iterations

    n_locs_in = length(smld.emitters)

    bagol_smld, diagnostics = SMLMBaGoL.run_bagol(smld;
        μ = cfg.μ,
        shape = cfg.shape,
        learn_distribution = cfg.learn_distribution,
        n_iterations = cfg.n_iterations,
        burn_in = cfg.burn_in,
        sync_interval = cfg.sync_interval,
        partition_sigma = cfg.partition_sigma,
        min_partition_size = cfg.min_partition_size,
        max_partition_size = cfg.max_partition_size,
        skip_partition_size = cfg.skip_partition_size,
        posterior_pixel_size = cfg.posterior_pixel_size,
        verbose = v >= Verbosity.PROGRESS,
    )

    n_emitters = diagnostics.n_emitters
    compression = n_locs_in > 0 ? round(n_locs_in / max(1, n_emitters), digits=1) : 0.0

    info = BaGoLInfo(n_locs_in, n_emitters, compression,
                     diagnostics.final_μ, diagnostics.final_shape,
                     diagnostics.n_partitions, diagnostics)

    # Diagnostic outputs via SMLMBaGoL report system
    if dir !== nothing && v >= Verbosity.STANDARD
        mkpath(dir)
        _save_config!(dir, cfg)
        _save_info!(dir, diagnostics)
        # Zero out track_id before compute_report — FC link IDs are not GT emitter assignments
        for e in smld.emitters
            e.track_id = 0
        end
        report = SMLMBaGoL.compute_report(bagol_smld, diagnostics; locs_smld=smld)
        SMLMBaGoL.write_report(report; output_dir=dir)
        SMLMBaGoL.plot_report(report; output_dir=dir)

        # BaGoL-specific renders: partition circles + overlay
        _render_bagol_diagnostics(smld, bagol_smld, diagnostics, dir, cfg)
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
)

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

    # Circles overlay: white localizations + red MAP-N emitters
    try
        (bg_img, _) = SMLMRender.render(smld; strategy=EllipseRender(),
            color=:white, target=target, clip_percentile=nothing)
        (fg_img, _) = SMLMRender.render(bagol_smld; strategy=EllipseRender(),
            color=:red, target=target, clip_percentile=nothing)
        combined = SMLMRender.compose(bg_img, fg_img; blend=:replace)
        SMLMRender.save_image(joinpath(dir, "circles.png"), combined)
    catch e
        @warn "BaGoL circles render failed" exception=e
    end

    # Partition-colored localizations
    if length(partition_ids) == length(smld.emitters)
        try
            for (i, e) in enumerate(smld.emitters)
                e.dataset = partition_ids[i]
            end
            part_smld = BasicSMLD(smld.emitters, smld.camera, smld.n_frames, 1)
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
                 outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD, kwargs...)
    t = @elapsed (bagol_smld, bagol_info) = bagol_step(smld, cfg;
        outdir=outdir, step_number=step_number, verbose=verbose)
    (bagol_smld, StepInfo(step_number, cfg, t, _step_summary(bagol_info); info=bagol_info))
end
