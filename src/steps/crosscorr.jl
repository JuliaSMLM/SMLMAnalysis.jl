"""
Cross-correlation g(r) step — pair cross-correlation between two channels.

Dispatches on `CrossCorrConfig <: AbstractMultiTargetStep` operating on
`Vector{BasicSMLD}`. Pass-through: SMLDs are not modified.

Computes the pair cross-correlation function g(r) in coordinate space using
KDTree range queries with Ripley's isotropic edge correction.
- g(r) > 1: spatial clustering / co-localization
- g(r) = 1: random (CSR)
- g(r) < 1: exclusion / anti-correlation
"""

"""
    CrossCorrConfig <: AbstractMultiTargetStep

Configuration for pair cross-correlation g(r) between two channels.

# Fields
- `r_max`: Maximum correlation distance in μm (default: 1.0)
- `dr`: Radial bin width in μm (default: 0.01)
- `edge_correction`: Apply Ripley's isotropic edge correction (default: true)
- `channels`: Which channels to correlate, as 1-based indices (default: (1, 2))
"""
@kwdef struct CrossCorrConfig <: AbstractMultiTargetStep
    r_max::Float64 = 1.0
    dr::Float64 = 0.01
    edge_correction::Bool = true
    channels::Tuple{Int,Int} = (1, 2)
end

step_name(::CrossCorrConfig) = "crosscorr"

"""
    crosscorr_step(smlds, cfg; outdir, step_number, verbose, labels) -> (smlds, CrossCorrInfo)

Compute pair cross-correlation g(r) between two channels. SMLDs pass through unmodified.
"""
function crosscorr_step(smlds::Vector{<:SMLMData.BasicSMLD}, cfg::CrossCorrConfig;
                         outdir::Union{String,Nothing}=nothing,
                         step_number::Int=0,
                         verbose::Int=Verbosity.STANDARD,
                         labels::Vector{Symbol}=Symbol[])
    v = verbose
    dir = step_outdir(outdir, step_number, cfg)

    ch_a, ch_b = cfg.channels
    (1 <= ch_a <= length(smlds) && 1 <= ch_b <= length(smlds)) ||
        error("CrossCorrConfig channels ($ch_a, $ch_b) out of range for $(length(smlds)) channels")
    ch_a != ch_b || error("CrossCorrConfig channels must be different (got $ch_a, $ch_b)")

    label_a = length(labels) >= ch_a ? labels[ch_a] : Symbol("Ch$ch_a")
    label_b = length(labels) >= ch_b ? labels[ch_b] : Symbol("Ch$ch_b")

    v >= Verbosity.PROGRESS && @info "[$step_number] crosscorr: $label_a × $label_b, r_max=$(cfg.r_max)μm, dr=$(cfg.dr)μm"

    smld_a = smlds[ch_a]
    smld_b = smlds[ch_b]

    local r, g, area
    t = @elapsed begin
        (r, g, area) = _compute_crosscorr(smld_a, smld_b, cfg)
    end

    n_a = length(smld_a.emitters)
    n_b = length(smld_b.emitters)

    if dir !== nothing
        mkpath(dir)
        _save_config!(dir, cfg)

        cc_info = CrossCorrInfo(r, g, n_a, n_b, area, cfg.r_max, cfg.dr, label_a, label_b, t)
        _save_info!(dir, cc_info)

        _write_crosscorr_csv(dir, r, g)

        if v >= Verbosity.STANDARD
            _write_crosscorr_stats(dir, cfg, cc_info)
            _save_crosscorr_plot(dir, r, g, label_a, label_b)
        end
    end

    v >= Verbosity.PROGRESS && @info "  -> g(r) computed: $(length(r)) bins, peak=$(round(maximum(g; init=0.0), digits=2)) ($(round(t, digits=2))s)"

    info = CrossCorrInfo(r, g, n_a, n_b, area, cfg.r_max, cfg.dr, label_a, label_b, t)
    (smlds, info)
end

_step_summary(info::CrossCorrInfo) = Dict{Symbol,Any}(
    :n_a => info.n_a,
    :n_b => info.n_b,
    :n_bins => length(info.r),
    :peak_g => round(maximum(info.g; init=0.0), digits=3),
    :channel_a => info.channel_a,
    :channel_b => info.channel_b,
)

"""
    analyze(smlds::Vector{BasicSMLD}, cfg::CrossCorrConfig; kwargs...) -> (smlds, StepInfo)

Multi-target dispatch: pair cross-correlation. SMLDs pass through.
"""
function analyze(smlds::Vector{<:SMLMData.BasicSMLD}, cfg::CrossCorrConfig;
                 outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD,
                 labels::Vector{Symbol}=Symbol[], kwargs...)
    t = @elapsed (smlds, cc_info) = crosscorr_step(smlds, cfg;
        outdir=outdir, step_number=step_number, verbose=verbose, labels=labels)
    (smlds, StepInfo(step_number, cfg, t, _step_summary(cc_info); info=cc_info))
end

# ============================================================
# Core algorithm
# ============================================================

"""
    _compute_crosscorr(smld_a, smld_b, cfg) -> (r, g, area)

Compute pair cross-correlation g(r) between two SMLDs using KDTree + Ripley's edge correction.
"""
function _compute_crosscorr(smld_a::SMLMData.BasicSMLD, smld_b::SMLMData.BasicSMLD, cfg::CrossCorrConfig)
    em_a = smld_a.emitters
    em_b = smld_b.emitters
    n_a = length(em_a)
    n_b = length(em_b)

    # FOV bounds from camera (use channel A's camera; should be the same)
    cam = smld_a.camera
    x_min = first(cam.pixel_edges_x)
    x_max = last(cam.pixel_edges_x)
    y_min = first(cam.pixel_edges_y)
    y_max = last(cam.pixel_edges_y)
    area = (x_max - x_min) * (y_max - y_min)
    bounds = (x_min, x_max, y_min, y_max)

    # Bin edges
    n_bins = ceil(Int, cfg.r_max / cfg.dr)
    r_edges = range(0.0, cfg.r_max, length=n_bins + 1)
    r_centers = [(r_edges[i] + r_edges[i+1]) / 2 for i in 1:n_bins]

    # Handle empty channels
    if n_a == 0 || n_b == 0
        return (collect(r_centers), ones(n_bins), area)
    end

    # Build KDTree from channel B
    coords_b = zeros(2, n_b)
    for i in 1:n_b
        coords_b[1, i] = em_b[i].x
        coords_b[2, i] = em_b[i].y
    end
    tree_b = KDTree(coords_b)

    # Accumulate weighted pair counts per bin
    counts = zeros(Float64, n_bins)
    density_b = n_b / area

    for i in 1:n_a
        xa, ya = em_a[i].x, em_a[i].y
        point = [xa, ya]
        idxs = inrange(tree_b, point, cfg.r_max)

        for j in idxs
            xb, yb = em_b[j].x, em_b[j].y
            dist = sqrt((xa - xb)^2 + (ya - yb)^2)
            dist == 0.0 && continue

            bin = floor(Int, dist / cfg.dr) + 1
            bin > n_bins && continue

            if cfg.edge_correction
                w = _ripley_edge_weight(xa, ya, dist, bounds)
                counts[bin] += 1.0 / w
            else
                counts[bin] += 1.0
            end
        end
    end

    # Normalize: g(r) = counts(r) / (n_A * density_B * annulus_area)
    g = zeros(Float64, n_bins)
    for k in 1:n_bins
        r_inner = r_edges[k]
        r_outer = r_edges[k+1]
        annulus_area = pi * (r_outer^2 - r_inner^2)
        expected = n_a * density_b * annulus_area
        g[k] = expected > 0 ? counts[k] / expected : 1.0
    end

    (collect(r_centers), g, area)
end

"""
    _ripley_edge_weight(x, y, r, bounds) -> Float64

Ripley's isotropic edge correction: fraction of circle circumference at (x,y)
with radius r that falls inside the rectangular FOV defined by bounds.

Returns a weight in (0, 1]. Points near boundaries get lower weights,
which are inverted to up-weight their contributions.
"""
function _ripley_edge_weight(x::Real, y::Real, r::Real, bounds::NTuple{4,<:Real})
    x_min, x_max, y_min, y_max = bounds

    # Distances to each boundary
    d_left   = x - x_min
    d_right  = x_max - x
    d_bottom = y - y_min
    d_top    = y_max - y

    # If circle fits entirely inside FOV, weight = 1
    if d_left >= r && d_right >= r && d_bottom >= r && d_top >= r
        return 1.0
    end

    # Compute the angle subtended by the arc outside each boundary
    # For a boundary at distance d < r, the arc outside subtends 2*acos(d/r)
    theta_outside = 0.0

    for d in (d_left, d_right, d_bottom, d_top)
        if d < r
            theta_outside += 2.0 * acos(clamp(d / r, -1.0, 1.0))
        end
    end

    # Corner corrections: if circle extends past two adjacent boundaries,
    # the corner region is double-counted. Check each of the 4 corners.
    # A corner at (cx, cy) contributes if both dx and dy < r AND
    # the corner is within the circle (sqrt(dx^2 + dy^2) < r)
    corners = ((d_left, d_bottom), (d_left, d_top),
               (d_right, d_bottom), (d_right, d_top))

    for (dx, dy) in corners
        if dx < r && dy < r && dx^2 + dy^2 < r^2
            # The corner overlap angle needs to be added back
            # (it was subtracted twice, once for each boundary)
            theta_outside -= acos(clamp(dx / r, -1.0, 1.0)) + acos(clamp(dy / r, -1.0, 1.0)) - pi/2
        end
    end

    theta_inside = 2pi - clamp(theta_outside, 0.0, 2pi)
    weight = theta_inside / (2pi)
    return max(weight, 0.01)  # Floor to avoid division by near-zero
end

# ============================================================
# Output helpers
# ============================================================

function _write_crosscorr_csv(dir::String, r::Vector{Float64}, g::Vector{Float64})
    filepath = joinpath(dir, "crosscorr_gr.csv")
    open(filepath, "w") do io
        println(io, "r,g")
        for i in eachindex(r)
            println(io, "$(r[i]),$(g[i])")
        end
    end
end

function _write_crosscorr_stats(dir::String, cfg::CrossCorrConfig, info::CrossCorrInfo)
    filepath = joinpath(dir, "stats.md")
    peak_g = maximum(info.g; init=0.0)
    peak_r = length(info.g) > 0 ? info.r[argmax(info.g)] : 0.0

    open(filepath, "w") do io
        println(io, "# Cross-Correlation g(r) Statistics\n")
        println(io, "## Summary")
        println(io, "- **Channel A**: $(info.channel_a) ($(info.n_a) localizations)")
        println(io, "- **Channel B**: $(info.channel_b) ($(info.n_b) localizations)")
        println(io, "- **FOV area**: $(round(info.area, digits=2)) μm²")
        println(io, "- **r_max**: $(cfg.r_max) μm")
        println(io, "- **dr**: $(cfg.dr) μm")
        println(io, "- **Edge correction**: $(cfg.edge_correction ? "Ripley's isotropic" : "none")")
        println(io, "- **Peak g(r)**: $(round(peak_g, digits=3)) at r = $(round(peak_r, digits=4)) μm")
        println(io, "- **Time**: $(round(info.elapsed_s, digits=2))s")
    end
end

function _save_crosscorr_plot(dir::String, r::Vector{Float64}, g::Vector{Float64},
                               label_a::Symbol, label_b::Symbol)
    isempty(r) && return

    peak_g = maximum(g)
    peak_r = r[argmax(g)]

    fig = Figure(size=(700, 450))
    ax = Axis(fig[1, 1],
        xlabel="r (μm)",
        ylabel="g(r)",
        title="Cross-correlation: $label_a × $label_b",
    )

    lines!(ax, r, g, color=:steelblue, linewidth=2)
    hlines!(ax, [1.0], color=:gray40, linestyle=:dash, linewidth=1, label="CSR (g=1)")

    # Annotate peak if above CSR
    if peak_g > 1.05
        scatter!(ax, [peak_r], [peak_g], color=:red, markersize=8)
        text!(ax, peak_r, peak_g,
            text="  g=$(round(peak_g, digits=2))",
            fontsize=11, color=:red, align=(:left, :bottom))
    end

    axislegend(ax, position=:rt, framevisible=false)
    save(joinpath(dir, "crosscorr_gr.png"), fig)
end
