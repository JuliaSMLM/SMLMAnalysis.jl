"""
    SMLMAnalysisDeepFitExt

Package extension that adds the `deepfit_training` and `deepfit_inference`
analyze() steps when `SMLMDeepFit` is loaded alongside `SMLMAnalysis`.

SMLMDeepFit is a weak dependency, so SMLMAnalysis's own dependency tree stays
free of its Reactant/Enzyme/Lux stack.

Steps:
  - `deepfit_training`: a calibration step, run standalone as `analyze(cfg::TrainConfig)`.
    It simulates training data from a PSF file and trains a DECODE U-Net. It returns
    the model (`TrainResult.model_path`), not an SMLD, so it does not go in the
    `steps` vector.
  - `deepfit_inference`: a localization step, `analyze(movie, cfg::DeepFitConfig)`.
    It takes DetectFitConfig's place in `AnalysisConfig.steps`: a raw `[H,W,T]`
    movie in, a `BasicSMLD` of `Emitter3DFit` out. The pipeline's camera replaces
    SMLMDeepFit's placeholder camera (`_prepare_step`).
"""
module SMLMAnalysisDeepFitExt

using SMLMData: SMLMData
using SMLMDeepFit: SMLMDeepFit
using SMLMRender: SMLMRender
using CairoMakie: CairoMakie
using SMLMAnalysis: SMLMAnalysis, Checkpoint, StepInfo, Verbosity, _save_box_overlay,
    _save_config!, _save_info!, _save_loc_per_frame, _save_step_smld, step_outdir
import SMLMAnalysis: _prepare_step, _produces_smld, analyze, step_name

# ------------------------------------------------------------
# helpers
# ------------------------------------------------------------

# Rebuild an immutable @kwdef config with one field replaced (positional constructor).
function _set_field(cfg::T, field::Symbol, val) where {T}
    field in fieldnames(T) || return cfg
    return T((f === field ? val : getfield(cfg, f) for f in fieldnames(T))...)
end

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

# ------------------------------------------------------------
# diagnostics (a plotting failure never fails the step)
# ------------------------------------------------------------

# Pixel size (µm) from the camera's uniform pixel edges; SMLMDeepFit's placeholder default otherwise.
function _deepfit_pixelsize(cam)
    if cam !== nothing && hasproperty(cam, :pixel_edges_x) && length(cam.pixel_edges_x) >= 2
        return Float64(cam.pixel_edges_x[2] - cam.pixel_edges_x[1])
    end
    return 0.1
end

# DECODE training loss and accuracy curves from TrainInfo.
function _deepfit_training_curves(dir, info)
    (hasproperty(info, :train_losses) && !isempty(info.train_losses)) || return
    ep = 1:length(info.train_losses)
    fig = CairoMakie.Figure(size = (940, 380))
    ax1 = CairoMakie.Axis(fig[1, 1], xlabel = "logged step", ylabel = "loss", title = "DECODE training loss")
    CairoMakie.lines!(ax1, ep, info.train_losses, label = "train")
    (hasproperty(info, :test_losses) && length(info.test_losses) == length(ep)) &&
        CairoMakie.lines!(ax1, ep, info.test_losses, label = "test")
    CairoMakie.axislegend(ax1)
    ax2 = CairoMakie.Axis(
        fig[1, 2], xlabel = "logged step", ylabel = "accuracy / efficiency", title = "DECODE accuracy"
    )
    hasproperty(info, :train_accuracies) && !isempty(info.train_accuracies) &&
        CairoMakie.lines!(ax2, ep, info.train_accuracies, label = "train")
    (hasproperty(info, :test_accuracies) && length(info.test_accuracies) == length(ep)) &&
        CairoMakie.lines!(ax2, ep, info.test_accuracies, label = "test")
    CairoMakie.axislegend(ax2)
    return CairoMakie.save(joinpath(dir, "loss_accuracy.png"), fig)
end

# The detectfit-style figure: a box centred on each localization, on sample movie frames.
function _deepfit_inference_overlay(dir, movie, smld, cam; box_size = 9)
    em = [e for e in smld.emitters if 1 <= e.frame <= size(movie, 3)]
    isempty(em) && return
    ps = _deepfit_pixelsize(cam)
    xc = Float64[e.x / ps - box_size / 2 for e in em]
    yc = Float64[e.y / ps - box_size / 2 for e in em]
    fr = Int[e.frame for e in em]
    colors = fill(:red, length(em))
    return _save_box_overlay(
        dir, "inference_overlay.png", movie, xc, yc, fr, Float64(box_size), colors;
        title_prefix = "frame",
        suptitle = "deepfit_inference: localizations (boxed) on movie frames"
    )
end

# Gaussian render of the inferred localizations.
function _deepfit_sr_render(dir, smld; zoom = 20, clip_percentile = 0.99)
    isempty(smld.emitters) && return
    (img, _) = SMLMRender.render(
        smld; strategy = SMLMRender.GaussianRender(), zoom = zoom, clip_percentile = clip_percentile
    )
    return SMLMRender.save_image(joinpath(dir, "inferred_sr.png"), img)
end

# ============================================================
# deepfit_training (calibration step; standalone, analyze(cfg))
# ============================================================
step_name(::SMLMDeepFit.TrainConfig) = "deepfit_training"
_produces_smld(::SMLMDeepFit.TrainConfig) = false

function analyze(
        cfg::SMLMDeepFit.TrainConfig;
        outdir = nothing, step_number::Int = 0, verbose::Int = Verbosity.STANDARD, kwargs...
    )
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
    verbose >= Verbosity.PROGRESS && @info "  → trained ($(round(t, digits = 1))s)"
    return (result, StepInfo(step_number, cfg, t, summary; info = _as_info(info)))
end

# ============================================================
# deepfit_inference (pipeline localization step; images -> BasicSMLD)
# ============================================================
step_name(::SMLMDeepFit.DeepFitConfig) = "deepfit_inference"

# Inject the pipeline's camera into the config, as for DetectFitConfig.
_prepare_step(cfg::SMLMDeepFit.DeepFitConfig, camera::SMLMData.AbstractCamera) =
    _set_field(cfg, :camera, camera)

# One dataset: a raw [H,W,T] movie.
function analyze(
        movie::AbstractArray{<:Real, 3}, cfg::SMLMDeepFit.DeepFitConfig;
        outdir = nothing, step_number::Int = 0, verbose::Int = Verbosity.STANDARD,
        checkpoint::Int = Checkpoint.EXPENSIVE, kwargs...
    )
    return _analyze_deepfit([movie], cfg, outdir, step_number, verbose, checkpoint)
end

# Several datasets: the pipeline's normalized state, a Vector of [H,W,T] movies.
function analyze(
        movies::AbstractVector, cfg::SMLMDeepFit.DeepFitConfig;
        outdir = nothing, step_number::Int = 0, verbose::Int = Verbosity.STANDARD,
        checkpoint::Int = Checkpoint.EXPENSIVE, kwargs...
    )
    return _analyze_deepfit(movies, cfg, outdir, step_number, verbose, checkpoint)
end

function _analyze_deepfit(movies, cfg, outdir, step_number, verbose, checkpoint)
    dir = step_outdir(outdir, step_number, cfg)
    verbose >= Verbosity.PROGRESS && @info "[$step_number] deepfit_inference" n_datasets = length(movies)
    smld = nothing
    info = nothing
    t = @elapsed begin
        per = [SMLMDeepFit.deepfit(m, cfg) for m in movies]
        smlds = SMLMData.BasicSMLD[p[1] for p in per]
        info = per[1][2]
        smld = _combine_datasets(smlds)
    end
    if dir !== nothing
        mkpath(dir)
        _save_config!(dir, cfg)
        try
            _save_info!(dir, info)
            _deepfit_inference_overlay(dir, movies[1], smld, cfg.camera)
            _deepfit_sr_render(dir, smld)
            _save_loc_per_frame(
                dir, smld; filename = "localizations_per_frame.png",
                title = "DeepFit localizations per frame"
            )
        catch err
            @warn "deepfit_inference: diagnostics failed (localizations still saved)" err
        end
    end
    checkpoint >= Checkpoint.EXPENSIVE && _save_step_smld(dir, smld; filename = "smld_deepfit.h5")
    verbose >= Verbosity.PROGRESS && @info "  → $(length(smld.emitters)) localizations ($(round(t, digits = 1))s)"
    return (smld, StepInfo(step_number, cfg, t, _numeric_summary(info); info = _as_info(info)))
end

# Only the single-dataset path is validated. Combining datasets needs an
# Emitter3DFit-aware dataset retag (SMLMData side), so it throws until then.
function _combine_datasets(smlds::Vector{<:SMLMData.BasicSMLD})
    length(smlds) == 1 && return smlds[1]
    return error(
        "deepfit_inference: combining datasets is not implemented yet " *
            "(needs an Emitter3DFit dataset retag); run each dataset on its own."
    )
end

end # module
