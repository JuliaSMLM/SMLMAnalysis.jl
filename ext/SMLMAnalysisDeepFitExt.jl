"""
    SMLMAnalysisDeepFitExt

Package extension that activates the `deepfit_training` and `deepfit_inference`
analyze() steps when `SMLMDeepFit` is loaded alongside `SMLMAnalysis`.

Kept as an extension (SMLMDeepFit is a weakdep) so SMLMAnalysis's core dependency
tree stays free of the Reactant/Enzyme/Lux deep-learning stack — downstream
consumers that do not `using SMLMDeepFit` never resolve it.

Step model:
  - `deepfit_training`  — calibration phase, standalone `analyze(cfg::TrainConfig)`.
    Self-generates simulated data from a PSF and trains a DECODE U-Net. Produces a
    model artifact (`TrainResult.model_path`), NOT an SMLD — does not belong in the
    linear smld-threading `steps` vector.
  - `deepfit_inference` — localization step, `analyze(movie, cfg::DeepFitConfig)`.
    Drops into the AnalysisConfig.steps vector in place of DetectFitConfig: raw
    `[H,W,T]` movie -> `BasicSMLD` of `Emitter3DFit`. The camera is injected from
    the pipeline (`_prepare_step`), replacing SMLMDeepFit's placeholder camera.
"""
module SMLMAnalysisDeepFitExt

using SMLMData
using SMLMDeepFit
import SMLMRender         # host dep: SR render of inferred localizations
import CairoMakie         # host dep: overlay / training-curve figures (qualified — dodges save/attributes clashes)
import SMLMAnalysis: analyze, step_name, _produces_smld, _prepare_step,
    step_outdir, _save_config!, _save_step_smld, _save_loc_per_frame, _save_info!,
    _save_box_overlay, StepInfo, Verbosity, Checkpoint

# ------------------------------------------------------------
# helpers
# ------------------------------------------------------------

"""Reconstruct an immutable @kwdef config with one field replaced (positional ctor)."""
function _set_field(cfg::T, field::Symbol, val) where {T}
    field in fieldnames(T) || return cfg
    T((f === field ? val : getfield(cfg, f) for f in fieldnames(T))...)
end

"""Defensive numeric-field summary for an upstream info struct."""
function _numeric_summary(info)
    d = Dict{Symbol,Any}()
    for f in fieldnames(typeof(info))
        v = getfield(info, f)
        v isa Number && (d[f] = v)
    end
    d
end

"""StepInfo's typed `info` slot requires <: AbstractSMLMInfo; guard the pre-reparent window."""
_as_info(x) = x isa SMLMData.AbstractSMLMInfo ? x : nothing

# ------------------------------------------------------------
# diagnostics (best-effort: a plotting failure never fails the step)
# ------------------------------------------------------------

"""Pixel size (µm) from an IdealCamera's uniform pixel edges; placeholder-camera fallback."""
function _deepfit_pixelsize(cam)
    if cam !== nothing && hasproperty(cam, :pixel_edges_x) && length(cam.pixel_edges_x) >= 2
        return Float64(cam.pixel_edges_x[2] - cam.pixel_edges_x[1])
    end
    0.1   # SMLMDeepFit placeholder camera default
end

"""DECODE training loss + accuracy curves from TrainInfo (its own LossPlot/AccuracyPlot ship empty)."""
function _deepfit_training_curves(dir, info)
    (hasproperty(info, :train_losses) && !isempty(info.train_losses)) || return
    ep = 1:length(info.train_losses)
    fig = CairoMakie.Figure(size = (940, 380))
    ax1 = CairoMakie.Axis(fig[1, 1], xlabel = "logged step", ylabel = "loss", title = "DECODE training loss")
    CairoMakie.lines!(ax1, ep, info.train_losses, label = "train")
    (hasproperty(info, :test_losses) && length(info.test_losses) == length(ep)) &&
        CairoMakie.lines!(ax1, ep, info.test_losses, label = "test")
    CairoMakie.axislegend(ax1)
    ax2 = CairoMakie.Axis(fig[1, 2], xlabel = "logged step", ylabel = "accuracy / efficiency", title = "DECODE accuracy")
    hasproperty(info, :train_accuracies) && !isempty(info.train_accuracies) &&
        CairoMakie.lines!(ax2, ep, info.train_accuracies, label = "train")
    (hasproperty(info, :test_accuracies) && length(info.test_accuracies) == length(ep)) &&
        CairoMakie.lines!(ax2, ep, info.test_accuracies, label = "test")
    CairoMakie.axislegend(ax2)
    CairoMakie.save(joinpath(dir, "loss_accuracy.png"), fig)
end

"""detectfit-style overlay: inferred localizations BOXED on sample movie frames, via the shared
`_save_box_overlay` core — so deepfit_inference reads as the same figure family as detectfit/filter
(boxes on the raw data), not bare × marks. Box is centered on each localization (µm→px via camera)."""
function _deepfit_inference_overlay(dir, movie, smld, cam; box_size = 9)
    em = [e for e in smld.emitters if 1 <= e.frame <= size(movie, 3)]
    isempty(em) && return
    ps = _deepfit_pixelsize(cam)
    xc = Float64[e.x / ps - box_size / 2 for e in em]   # box corner (px), centered on the localization
    yc = Float64[e.y / ps - box_size / 2 for e in em]
    fr = Int[e.frame for e in em]
    colors = fill(:red, length(em))
    _save_box_overlay(dir, "inference_overlay.png", movie, xc, yc, fr, Float64(box_size), colors;
                      title_prefix = "frame",
                      suptitle = "deepfit_inference: localizations (boxed) on movie frames")
end

"""Super-resolution Gaussian render of the inferred localizations (sparse on sparse data — clip is standard)."""
function _deepfit_sr_render(dir, smld; zoom = 20, clip_percentile = 0.99)
    isempty(smld.emitters) && return
    (img, _) = SMLMRender.render(smld; strategy = SMLMRender.GaussianRender(),
                                 zoom = zoom, clip_percentile = clip_percentile)
    SMLMRender.save_image(joinpath(dir, "inferred_sr.png"), img)
end

# ============================================================
# deepfit_training  (calibration phase; standalone, f(cfg))
# ============================================================
step_name(::SMLMDeepFit.TrainConfig) = "deepfit_training"
_produces_smld(::SMLMDeepFit.TrainConfig) = false

function analyze(cfg::SMLMDeepFit.TrainConfig;
        outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD, kwargs...)
    dir = step_outdir(outdir, step_number, cfg)
    verbose >= Verbosity.PROGRESS && @info "[$step_number] deepfit_training"
    t = @elapsed ((result, info) = SMLMDeepFit.train(cfg))
    if dir !== nothing
        mkpath(dir)
        _save_config!(dir, cfg)
        try
            _save_info!(dir, info)
            _deepfit_training_curves(dir, info)
        catch err
            @warn "deepfit_training: diagnostics failed (model still saved)" err
        end
    end
    summary = _numeric_summary(info)
    hasproperty(result, :model_path) && (summary[:model_path] = result.model_path)
    verbose >= Verbosity.PROGRESS && @info "  → trained ($(round(t, digits=1))s)"
    (result, StepInfo(step_number, cfg, t, summary; info=_as_info(info)))
end

# ============================================================
# deepfit_inference  (pipeline localization step; images -> BasicSMLD)
# ============================================================
step_name(::SMLMDeepFit.DeepFitConfig) = "deepfit_inference"

"""Inject the pipeline-level camera into the config (mirrors DetectFitConfig)."""
_prepare_step(cfg::SMLMDeepFit.DeepFitConfig, camera::SMLMData.AbstractCamera) =
    _set_field(cfg, :camera, camera)

# Single dataset: a raw [H,W,T] movie.
function analyze(movie::AbstractArray{<:Real,3}, cfg::SMLMDeepFit.DeepFitConfig;
        outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD,
        checkpoint::Int=Checkpoint.EXPENSIVE, kwargs...)
    _analyze_deepfit([movie], cfg, outdir, step_number, verbose, checkpoint)
end

# Multi-dataset: the pipeline's normalized state, a Vector of [H,W,T] movies.
function analyze(movies::AbstractVector, cfg::SMLMDeepFit.DeepFitConfig;
        outdir=nothing, step_number::Int=0, verbose::Int=Verbosity.STANDARD,
        checkpoint::Int=Checkpoint.EXPENSIVE, kwargs...)
    _analyze_deepfit(movies, cfg, outdir, step_number, verbose, checkpoint)
end

function _analyze_deepfit(movies, cfg, outdir, step_number, verbose, checkpoint)
    dir = step_outdir(outdir, step_number, cfg)
    verbose >= Verbosity.PROGRESS && @info "[$step_number] deepfit_inference" n_datasets=length(movies)
    smld = nothing
    info = nothing
    t = @elapsed begin
        per   = [SMLMDeepFit.deepfit(m, cfg) for m in movies]
        smlds = SMLMData.BasicSMLD[p[1] for p in per]
        info  = per[1][2]
        smld  = _combine_datasets(smlds)
    end
    if dir !== nothing
        mkpath(dir)
        _save_config!(dir, cfg)
        try
            _save_info!(dir, info)
            _deepfit_inference_overlay(dir, movies[1], smld, cfg.camera)
            _deepfit_sr_render(dir, smld)
            _save_loc_per_frame(dir, smld; filename="localizations_per_frame.png",
                                title="DeepFit localizations per frame")
        catch err
            @warn "deepfit_inference: diagnostics failed (localizations still saved)" err
        end
    end
    checkpoint >= Checkpoint.EXPENSIVE && _save_step_smld(dir, smld; filename="smld_deepfit.jld2")
    verbose >= Verbosity.PROGRESS && @info "  → $(length(smld.emitters)) localizations ($(round(t, digits=1))s)"
    (smld, StepInfo(step_number, cfg, t, _numeric_summary(info); info=_as_info(info)))
end

# Single dataset is the validated path. Multi-dataset combine needs an
# Emitter3DFit-aware dataset retag (SMLMData side) — wire after single-dataset
# reconstruction is green rather than ship an untested emitter reconstruction.
function _combine_datasets(smlds::Vector{<:SMLMData.BasicSMLD})
    length(smlds) == 1 && return smlds[1]
    error("deepfit_inference: multi-dataset combine not yet wired " *
          "(need Emitter3DFit dataset retag); run datasets individually for now.")
end

end # module
