"""
Checkpoint persistence using JLD2.

Saves analysis state after each step to enable cross-session resume.
Checkpoints are stored in outdir/checkpoints/step_NNN_name.jld2

Saved state:
- smld, smld_raw, smld_connected, bagol_smld (as columnar arrays for fast serialization)
- drift_model, bagol_result
- roi_datasets (dataset index for each ROI)
- step_history (Vector{StepRecord})
- camera
- image_source (path or :memory)
- n_datasets, n_frames_per_dataset

Not saved (reconstructable):
- images (reload from path)
- roi_batch (can reconstruct from smld if needed)

Versions:
- v3: Uses columnar format for emitters (100-1000x faster serialization)
- v4: Added bagol_result (BaGoLDiagnostics)
- v5: Added bagol_smld (grouped emitters from BaGoL)
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
# Checkpoint directory management
# ============================================================

function _checkpoint_dir(a::Analysis)
    a.outdir === nothing && return nothing
    joinpath(a.outdir, "checkpoints")
end

function _checkpoint_path(a::Analysis, step::Int, name::String)
    dir = _checkpoint_dir(a)
    dir === nothing && return nothing
    joinpath(dir, "step_$(lpad(step, 3, '0'))_$(name).jld2")
end

function _checkpoint_path(a::Analysis)
    isempty(a.steps) && return nothing
    _checkpoint_path(a, a.step_counter, a.steps[end].name)
end

# ============================================================
# Save checkpoint to disk
# ============================================================

"""
    _save_checkpoint!(a::Analysis)

Save current analysis state to disk. Called automatically after each step
when `a.checkpoint == true` and `a.outdir` is set.

Uses columnar format for SMLD data (v3) for fast serialization of millions of emitters.
"""
function _save_checkpoint!(a::Analysis)
    path = _checkpoint_path(a)
    path === nothing && return nothing

    # Ensure directory exists
    mkpath(dirname(path))

    # Determine image source
    image_source = a.data.path !== nothing ? a.data.path : :memory

    # Convert SMLD to columnar format for fast serialization
    smld_cols = _smld_to_columnar(a.smld)
    smld_raw_cols = _smld_to_columnar(a.smld_raw)
    smld_connected_cols = _smld_to_columnar(a.smld_connected)
    bagol_smld_cols = _smld_to_columnar(a.bagol_smld)

    # Save state (exclude roi_batch - can reconstruct from smld if needed)
    jldsave(path;
        smld_cols = smld_cols,
        smld_raw_cols = smld_raw_cols,
        smld_connected_cols = smld_connected_cols,
        bagol_smld_cols = bagol_smld_cols,
        drift_model = a.drift_model,
        bagol_result = a.bagol_result,
        roi_datasets = a.roi_datasets,
        steps = a.steps,
        step_counter = a.step_counter,
        camera = a.camera,
        image_source = image_source,
        verbose = a.verbose,
        n_datasets = a.n_datasets,
        n_frames_per_dataset = a.n_frames_per_dataset,
        checkpoint_version = 5  # v5: added bagol_smld
    )

    path
end

# ============================================================
# Load checkpoint from disk
# ============================================================

"""
    _load_checkpoint(path::String) -> NamedTuple

Load checkpoint data from disk. Returns a NamedTuple with all saved fields.
Handles v1, v2 (direct SMLD), and v3 (columnar SMLD) formats.
"""
function _load_checkpoint(path::String)
    !isfile(path) && error("Checkpoint not found: $path")

    jldopen(path, "r") do file
        # Check version for backward compatibility
        version = haskey(file, "checkpoint_version") ? file["checkpoint_version"] : 1

        # Load SMLD based on version
        if version >= 3
            # v3+: columnar format
            smld = _columnar_to_smld(file["smld_cols"])
            smld_raw = _columnar_to_smld(file["smld_raw_cols"])
            smld_connected = _columnar_to_smld(file["smld_connected_cols"])
            bagol_smld = version >= 5 && haskey(file, "bagol_smld_cols") ? _columnar_to_smld(file["bagol_smld_cols"]) : nothing
        else
            # v1/v2: direct SMLD objects
            smld = file["smld"]
            smld_raw = file["smld_raw"]
            smld_connected = file["smld_connected"]
            bagol_smld = nothing
        end

        (
            smld = smld,
            smld_raw = smld_raw,
            smld_connected = smld_connected,
            bagol_smld = bagol_smld,
            drift_model = file["drift_model"],
            bagol_result = version >= 4 && haskey(file, "bagol_result") ? file["bagol_result"] : nothing,
            roi_datasets = version >= 2 && haskey(file, "roi_datasets") ? file["roi_datasets"] : nothing,
            steps = file["steps"],
            step_counter = file["step_counter"],
            camera = file["camera"],
            image_source = file["image_source"],
            verbose = file["verbose"],
            n_datasets = version >= 2 && haskey(file, "n_datasets") ? file["n_datasets"] : 1,
            n_frames_per_dataset = version >= 2 && haskey(file, "n_frames_per_dataset") ? file["n_frames_per_dataset"] : nothing
        )
    end
end

"""
    _find_checkpoint(a::Analysis, step::Int) -> Union{String, Nothing}

Find checkpoint file for a given step number.
"""
function _find_checkpoint(a::Analysis, step::Int)
    dir = _checkpoint_dir(a)
    dir === nothing && return nothing
    !isdir(dir) && return nothing

    # Look for step_NNN_*.jld2
    pattern = "step_$(lpad(step, 3, '0'))_"
    for f in readdir(dir)
        if startswith(f, pattern) && endswith(f, ".jld2")
            return joinpath(dir, f)
        end
    end
    nothing
end

"""
    _list_checkpoints(outdir::String) -> Vector{Tuple{Int, String, String}}

List all checkpoints in outdir/checkpoints/. Returns vector of (step, name, path) tuples.
"""
function _list_checkpoints(outdir::String)
    dir = joinpath(outdir, "checkpoints")
    !isdir(dir) && return Tuple{Int, String, String}[]

    checkpoints = Tuple{Int, String, String}[]
    for f in readdir(dir)
        m = match(r"step_(\d+)_(.+)\.jld2$", f)
        if m !== nothing
            step = parse(Int, m.captures[1])
            name = m.captures[2]
            push!(checkpoints, (step, name, joinpath(dir, f)))
        end
    end
    sort!(checkpoints, by=x->x[1])
    checkpoints
end

# ============================================================
# Resume analysis from disk
# ============================================================

"""
    resume_analysis(outdir::String; images=nothing, step::Union{Int,Nothing}=nothing) -> Analysis

Resume an analysis from a checkpoint directory.

# Arguments
- `outdir`: Path to the analysis output directory (containing checkpoints/)
- `images`: Optional image array. If not provided and original was from a file,
            images will be reloaded from that path.
- `step`: Optional step number to resume from. Default: latest checkpoint.

# Returns
Analysis object ready to continue with additional steps.

# Example
```julia
# Resume from latest checkpoint
a = resume_analysis("output/")

# Resume from specific step
a = resume_analysis("output/"; step=2)

# Resume with new images (e.g., subset)
a = resume_analysis("output/"; images=new_images)
```
"""
function resume_analysis(outdir::String; images=nothing, step::Union{Int,Nothing}=nothing)
    checkpoints = _list_checkpoints(outdir)
    isempty(checkpoints) && error("No checkpoints found in $outdir/checkpoints/")

    # Find the checkpoint to load
    if step === nothing
        # Latest checkpoint
        _, _, path = last(checkpoints)
    else
        # Specific step
        idx = findfirst(c -> c[1] == step, checkpoints)
        idx === nothing && error("No checkpoint at step $step. Available: $(first.(checkpoints))")
        _, _, path = checkpoints[idx]
    end

    # Load checkpoint data
    cp = _load_checkpoint(path)

    # Resolve images
    if images === nothing
        if cp.image_source === :memory
            error("Original images were in-memory. Provide `images` argument to resume.")
        end
        data = DataSource(cp.image_source)
    else
        data = DataSource(images)
    end

    # Determine n_frames_per_dataset: from checkpoint or compute from images
    n_frames_per_dataset = cp.n_frames_per_dataset
    if n_frames_per_dataset === nothing && images !== nothing
        total_frames = size(images, 3)
        n_frames_per_dataset = div(total_frames, cp.n_datasets)
    elseif n_frames_per_dataset === nothing
        # Fallback: try to get from smld metadata if available
        n_frames_per_dataset = cp.smld !== nothing ? cp.smld.n_frames : 1
    end

    # Reconstruct Analysis
    a = Analysis(
        data,
        cp.camera,
        cp.n_datasets,
        n_frames_per_dataset,
        nothing,  # roi_batch - not saved
        cp.roi_datasets,
        cp.smld_raw,
        cp.smld,
        cp.smld_connected,
        cp.drift_model,
        cp.bagol_result,
        cp.bagol_smld,
        Dict{Int, AnalysisCheckpoint}(),  # Will be populated from disk on demand
        cp.steps,
        cp.step_counter,
        outdir,
        cp.verbose,
        true  # Enable checkpointing since we're resuming
    )

    @info "Resumed analysis from step $(cp.step_counter): $(cp.steps[end].name)"
    a
end
