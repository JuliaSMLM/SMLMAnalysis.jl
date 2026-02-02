"""
BaGoL (Bayesian Grouping of Localizations) step for SMLMAnalysis.

This step groups localizations into emitters using Bayesian inference.
Requires frame-connected and calibrated data (a.smld).
"""

using SMLMBaGoL

"""
    BaGoLConfig

Configuration for BaGoL grouping step.

# Fields
- `name`: Step name (default: "bagol")
- `verbose`: Verbosity level (default: Verbosity.STANDARD)
- `n_iterations`: Number of MCMC iterations (default: 10000)
- `burn_in`: Burn-in iterations before recording samples (default: 2000)
- `α`: Shape parameter for count distribution. Float64 or :auto (default: :auto)
- `learn_α`: Whether to update α during MCMC (default: true)
- `λ_K`: Poisson prior mean for emitter count K. Nothing = use default (n_locs/5) (default: nothing)
- `partition_threshold`: Use partitioning if n_locs > threshold (default: 500)
- `nsigma`: DBSCAN threshold in sigma units (default: 3.0)
- `min_partition_size`: Minimum locs per partition (smaller dropped as noise) (default: 0)
- `max_partition_size`: Maximum locs per partition (default: 1000)
- `render_zoom`: Zoom factor for renders (default: 50)
"""
Base.@kwdef struct BaGoLConfig <: StepConfig
    name::String = "bagol"
    verbose::Int = Verbosity.STANDARD
    # RJMCMC parameters
    n_iterations::Int = 10000
    burn_in::Int = 2000
    α::Union{Float64, Symbol} = :auto
    learn_α::Bool = true
    λ_K::Union{Float64, Nothing} = nothing  # Poisson prior mean for K
    # Partitioning
    partition_threshold::Int = 500
    nsigma::Float64 = 3.0
    min_partition_size::Int = 0  # Keep ALL clusters (was 10, dropped most emitters)
    max_partition_size::Int = 1000
    # Output
    render_zoom::Int = 50  # 50x for ellipse/circle plots (2nm effective pixel)
end

"""
    run_step!(a::Analysis, cfg::BaGoLConfig)

Run BaGoL grouping on calibrated, frame-connected localizations.

Requires:
- `a.smld` to contain calibrated localizations (from frameconnect + calibration)

Produces:
- `a.bagol_result`: BaGoLDiagnostics with posterior info
- Updates a.smld with grouped emitters (optional, controlled by replace_smld)
- Renders via SMLMRender if verbosity >= STANDARD
"""
function run_step!(a::Analysis, cfg::BaGoLConfig)
    a.smld === nothing && error("No smld - run detectfit, filter, and frameconnect first")

    a.step_counter += 1
    t0 = time()
    verbose = _get_verbose(a, cfg)

    if verbose >= Verbosity.PROGRESS
        @info "Step $(a.step_counter): $(cfg.name)"
    end

    n_locs = length(a.smld.emitters)

    # Compute λ_K if not specified (default: n_locs / 5)
    λ_K = cfg.λ_K === nothing ? Float64(n_locs) / 5.0 : cfg.λ_K

    # Run BaGoL - returns (BasicSMLD, BaGoLDiagnostics)
    bagol_smld, diagnostics = SMLMBaGoL.run_bagol(
        a.smld;
        n_iterations = cfg.n_iterations,
        burn_in = cfg.burn_in,
        α = cfg.α,
        learn_α = cfg.learn_α,
        λ_K = λ_K,
        partition_threshold = cfg.partition_threshold,
        nsigma = cfg.nsigma,
        min_partition_size = cfg.min_partition_size,
        max_partition_size = cfg.max_partition_size,
        verbose = verbose >= Verbosity.PROGRESS
    )

    # Store results for QC and downstream use
    a.bagol_result = diagnostics
    a.bagol_smld = bagol_smld
    n_emitters = diagnostics.n_emitters

    # Record timing and summary
    t = time() - t0
    summary = Dict{Symbol, Any}(
        :n_locs_in => n_locs,
        :n_emitters => n_emitters,
        :compression => n_locs > 0 ? round(n_locs / max(1, n_emitters), digits=1) : 0.0
    )
    # diagnostics is the BaGoL info struct (tuple-pattern)
    _record!(a, cfg, t, summary; info=diagnostics)

    # Save outputs using SMLMRender
    if verbose >= Verbosity.STANDARD
        _save_bagol_outputs!(a, cfg, bagol_smld, diagnostics, verbose)
    end

    _checkpoint!(a)

    if verbose >= Verbosity.PROGRESS
        @info "  BaGoL: $(n_locs) locs → $(n_emitters) emitters ($(summary[:compression])x compression)"
    end

    a
end

"""
Save BaGoL output renders and statistics using SMLMRender.

Renders generated:
- `overlay.png`: PRIMARY QC - Input (gray) + Output (red) ellipses showing σ
- `input_gaussian.png`: Input SR image (inferno colormap)
- `bagol_gaussian.png`: Output SR image (viridis colormap)
- `k_posterior.png`: P(K|data) histogram
- `stats.md`: Summary with precision improvement metrics
"""
function _save_bagol_outputs!(a::Analysis, cfg::BaGoLConfig, bagol_smld::BasicSMLD, diagnostics, verbose::Int)
    stepdir = _stepdir(a, cfg)
    stepdir === nothing && return
    mkpath(stepdir)

    # Save config
    _save_config!(stepdir, cfg)

    zoom = cfg.render_zoom

    # 1. PRIMARY QC: Input/Output Overlay with ellipses showing σ
    # Gray ellipses = input locs (many, small σ)
    # Red ellipses = output emitters (few, even smaller σ from combining)
    # Tuple-pattern: render returns (image, RenderInfo), ignore here
    if diagnostics.n_emitters > 0
        _ = SMLMRender.render([a.smld, bagol_smld];
            colors = [:gray, :red],
            strategy = EllipseRender(),
            zoom = zoom,
            filename = joinpath(stepdir, "overlay.png")
        )
    end

    # 2. Gaussian renders for detail/publication
    _ = SMLMRender.render(a.smld;
        strategy = GaussianRender(),
        zoom = zoom,
        colormap = :inferno,
        filename = joinpath(stepdir, "input_gaussian.png")
    )

    if diagnostics.n_emitters > 0
        _ = SMLMRender.render(bagol_smld;
            strategy = GaussianRender(),
            zoom = zoom,
            colormap = :viridis,
            filename = joinpath(stepdir, "bagol_gaussian.png")
        )
    end

    # 3. K posterior histogram
    _plot_k_posterior(diagnostics, joinpath(stepdir, "k_posterior.png"))

    # 4. Save summary stats with precision comparison
    _save_bagol_stats!(stepdir, a.smld, bagol_smld, diagnostics)
end

"""
Save BaGoL statistics including precision improvement metrics.
"""
function _save_bagol_stats!(stepdir::String, input_smld::BasicSMLD, bagol_smld::BasicSMLD, diagnostics)
    open(joinpath(stepdir, "stats.md"), "w") do io
        println(io, "# BaGoL Results\n")

        n_in = length(input_smld.emitters)
        n_out = diagnostics.n_emitters

        println(io, "## Grouping")
        println(io, "- Input localizations: $n_in")
        println(io, "- Output emitters: $n_out")
        println(io, "- Compression: $(round(n_in / max(1, n_out), digits=1))x")
        println(io, "")

        println(io, "## Model Parameters")
        println(io, "- Final μ: $(round(diagnostics.final_μ, digits=2)) (mean locs/emitter)")
        println(io, "- Final α: $(round(diagnostics.final_α, digits=2))")
        println(io, "")

        # Precision improvement comparison
        println(io, "## Precision")
        if n_in > 0
            σ_in = mean(sqrt(e.σ_x^2 + e.σ_y^2) for e in input_smld.emitters)
            println(io, "- Mean input σ: $(round(σ_in * 1000, digits=1)) nm")
        end
        if n_out > 0 && !isempty(bagol_smld.emitters)
            σ_out = mean(sqrt(e.σ_x^2 + e.σ_y^2) for e in bagol_smld.emitters)
            println(io, "- Mean output σ: $(round(σ_out * 1000, digits=1)) nm")
            if n_in > 0
                improvement = σ_in / σ_out
                println(io, "- Precision improvement: $(round(improvement, digits=1))x")
            end
        end
    end
end

"""
Plot K posterior histogram.
"""
function _plot_k_posterior(diagnostics, filepath::String)
    fig = Figure(size=(600, 400))
    ax = Axis(fig[1, 1],
              xlabel="K (number of emitters)",
              ylabel="Count",
              title="Posterior P(K|data)")

    posterior_k = diagnostics.posterior_k
    ks = 0:(length(posterior_k) - 1)
    barplot!(ax, ks, posterior_k, color=:steelblue)

    # Mark MAP
    map_k = argmax(posterior_k) - 1
    vlines!(ax, [map_k], color=:red, linestyle=:dash, label="MAP K = $map_k")
    axislegend(ax, position=:rt)

    save(filepath, fig)
end
