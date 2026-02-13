"""
Checkpoint persistence using JLD2.

Provides save/load for SMLD data using fast columnar format.
In the functional pipeline, users save intermediate SMLDs directly:

```julia
# Save after expensive detectfit
(smld, info) = analyze(image_stacks, DetectFitConfig(camera=cam, boxer=BoxerConfig(boxsize=9)))
save_smld("output/detectfit.h5", smld)

# Resume later
smld = load_smld("output/detectfit.h5")
(smld, info) = analyze(smld, FilterConfig(photons=(500.0, Inf)))
```

For full pipeline state save/restore (including drift model, connected SMLD, etc.),
use `save_pipeline_state` / `load_pipeline_state`.
"""

using JLD2

# ============================================================
# Columnar SMLD serialization (fast for millions of emitters)
# ============================================================

"""
Convert SMLD emitters to columnar format for fast JLD2 serialization.
Returns a NamedTuple with arrays for each field.
"""
function _smld_to_columnar(smld::Union{BasicSMLD, Nothing})
    smld === nothing && return nothing
    isempty(smld.emitters) && return (
        emitter_type = Nothing,
        n = 0,
        camera = smld.camera,
        n_frames = smld.n_frames,
        n_datasets = smld.n_datasets,
        metadata = smld.metadata
    )

    e = smld.emitters
    n = length(e)
    T = typeof(e[1])
    fields = fieldnames(T)

    # Build columnar data dict
    cols = Dict{Symbol, Any}()
    for f in fields
        vals = [getfield(em, f) for em in e]
        cols[f] = vals
    end

    (
        emitter_type = T,
        n = n,
        columns = cols,
        camera = smld.camera,
        n_frames = smld.n_frames,
        n_datasets = smld.n_datasets,
        metadata = smld.metadata
    )
end

"""
Reconstruct SMLD from columnar format.
"""
function _columnar_to_smld(cols)
    cols === nothing && return nothing
    cols.n == 0 && return BasicSMLD(
        AbstractEmitter[],
        cols.camera,
        cols.n_frames,
        cols.n_datasets,
        cols.metadata
    )

    T = cols.emitter_type
    n = cols.n
    fields = fieldnames(T)

    # Reconstruct emitters
    emitters = Vector{T}(undef, n)
    for i in 1:n
        args = [cols.columns[f][i] for f in fields]
        emitters[i] = T(args...)
    end

    BasicSMLD(emitters, cols.camera, cols.n_frames, cols.n_datasets, cols.metadata)
end

# ============================================================
# Pipeline state save/load
# ============================================================

"""
    save_pipeline_state(path::String, result::AnalysisResult;
                        smld_raw=nothing, step_records=StepRecord[], camera=nothing)

Save pipeline state to JLD2 for cross-session resume.

# Example
```julia
(result, info) = analyze(image_stacks, config)
save_pipeline_state("output/checkpoint.jld2", result;
    step_records=info.step_records, camera=config.camera)
```
"""
function save_pipeline_state(path::String, result::AnalysisResult;
                             smld_raw::Union{BasicSMLD,Nothing}=nothing,
                             step_records::Vector{StepRecord}=StepRecord[],
                             camera::Union{SMLMData.AbstractCamera,Nothing}=nothing)
    mkpath(dirname(path))

    smld_cols = _smld_to_columnar(result.smld)
    smld_raw_cols = _smld_to_columnar(smld_raw)
    smld_connected_cols = _smld_to_columnar(result.smld_connected)

    jldsave(path;
        smld_cols = smld_cols,
        smld_raw_cols = smld_raw_cols,
        smld_connected_cols = smld_connected_cols,
        drift_model = result.drift_model,
        step_records = step_records,
        camera = camera,
        checkpoint_version = 7
    )

    path
end

"""
    load_pipeline_state(path::String) -> NamedTuple

Load pipeline state from JLD2.

Returns a NamedTuple with fields: smld, smld_raw, smld_connected, drift_model,
step_records, camera.

# Example
```julia
state = load_pipeline_state("output/checkpoint.jld2")
smld = state.smld
# Continue pipeline from here
(smld, info) = analyze(smld, FilterConfig(photons=(500.0, Inf)))
```
"""
function load_pipeline_state(path::String)
    !isfile(path) && error("Checkpoint not found: $path")

    jldopen(path, "r") do file
        version = haskey(file, "checkpoint_version") ? file["checkpoint_version"] : 1

        # Load SMLD based on version
        if version >= 3
            smld = _columnar_to_smld(file["smld_cols"])
            smld_raw = _columnar_to_smld(file["smld_raw_cols"])
            smld_connected = _columnar_to_smld(file["smld_connected_cols"])
        else
            smld = file["smld"]
            smld_raw = file["smld_raw"]
            smld_connected = file["smld_connected"]
        end

        camera = haskey(file, "camera") ? file["camera"] : nothing
        step_records = if version >= 7 && haskey(file, "step_records")
            file["step_records"]
        elseif haskey(file, "steps")
            file["steps"]  # Legacy format
        else
            StepRecord[]
        end

        (
            smld = smld,
            smld_raw = smld_raw,
            smld_connected = smld_connected,
            drift_model = file["drift_model"],
            step_records = step_records,
            camera = camera
        )
    end
end
