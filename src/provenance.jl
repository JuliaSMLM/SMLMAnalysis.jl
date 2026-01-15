"""
    provenance.jl

Types and functions for tracking SMLM analysis workflows for reproducibility.
"""

using Dates

"""
    ProcessingStep

Records a single step in an SMLM analysis pipeline.

# Fields
- `step_name::String`: Descriptive name of the processing step
- `function_name::Symbol`: The actual function called
- `parameters::Dict{Symbol,Any}`: Parameters used in this step
- `timestamp::DateTime`: When this step was executed
- `package::Symbol`: Which package provided the function
- `output_summary::String`: Brief description of output (type, size, etc.)
"""
struct ProcessingStep
    step_name::String
    function_name::Symbol
    parameters::Dict{Symbol,Any}
    timestamp::DateTime
    package::Symbol
    output_summary::String
end

"""
    SMLMWorkflow

Tracks an SMLM analysis pipeline for reproducibility and provenance.

Records each processing step with parameters, timestamps, and data flow.
Enables complete reconstruction of analysis parameters for publications
and reproducible research.

# Fields
- `steps::Vector{ProcessingStep}`: Ordered list of processing steps
- `created::DateTime`: When the workflow was created
- `description::String`: User description of the analysis
- `metadata::Dict{String,Any}`: Additional workflow-level metadata

# Usage Patterns

## Pattern 1: Automatic tracking (returned from workflows)
```julia
result, workflow = localization_workflow(images, camera)
println(workflow)  # See what was done
save_workflow(workflow, "analysis.json")
```

## Pattern 2: Manual step tracking
```julia
workflow = SMLMWorkflow("My analysis")
add_step!(workflow, "Simulation", :simulate,
          Dict(:density => 1.0, :psf => 0.13), :SMLMSim,
          "BasicSMLD with 1000 emitters")
```

## Pattern 3: Pass workflow through pipeline
```julia
wf = SMLMWorkflow("Complete pipeline")
result1 = simulate(params; workflow=wf)
result2 = detect(result1; workflow=wf)
# wf now contains complete history
```
"""
mutable struct SMLMWorkflow
    steps::Vector{ProcessingStep}
    created::DateTime
    description::String
    metadata::Dict{String,Any}

    function SMLMWorkflow(description::String="SMLM Analysis";
                         metadata::Dict{String,Any}=Dict{String,Any}())
        new(ProcessingStep[], now(), description, metadata)
    end
end

"""
    add_step!(workflow::SMLMWorkflow, step_name, function_name, parameters,
              package, output_summary) → ProcessingStep

Add a processing step to the workflow with current timestamp.

# Arguments
- `workflow`: The workflow to update
- `step_name`: Human-readable step name (e.g., "Particle Detection")
- `function_name`: Symbol of function called (e.g., :getboxes)
- `parameters`: Dict of parameters used
- `package`: Symbol of source package (e.g., :SMLMBoxer)
- `output_summary`: Brief description of output

# Returns
The created ProcessingStep
"""
function add_step!(workflow::SMLMWorkflow, step_name::String, function_name::Symbol,
                   parameters::Dict{Symbol,Any}, package::Symbol, output_summary::String)
    step = ProcessingStep(step_name, function_name, parameters, now(), package, output_summary)
    push!(workflow.steps, step)
    return step
end

"""
    add_step!(workflow::Union{SMLMWorkflow,Nothing}, ...) → ProcessingStep

Convenience method that does nothing if workflow is `nothing`.
Enables optional workflow tracking without conditional code.
"""
function add_step!(workflow::Nothing, step_name::String, function_name::Symbol,
                   parameters::Dict{Symbol,Any}, package::Symbol, output_summary::String)
    return nothing
end

# Pretty printing for workflow
function Base.show(io::IO, workflow::SMLMWorkflow)
    println(io, "SMLMWorkflow: $(workflow.description)")
    println(io, "Created: $(Dates.format(workflow.created, "yyyy-mm-dd HH:MM:SS"))")
    println(io, "Steps: $(length(workflow.steps))")
    for (i, step) in enumerate(workflow.steps)
        elapsed = if i == 1
            ""
        else
            dt = (step.timestamp - workflow.steps[i-1].timestamp).value / 1000
            " (+$(round(dt, digits=2))s)"
        end
        println(io, "  $i. $(step.step_name)$elapsed")
        println(io, "     $(step.package).$(step.function_name)")
        if !isempty(step.parameters)
            println(io, "     Parameters: ", join(["$k=$v" for (k,v) in step.parameters], ", "))
        end
        println(io, "     → $(step.output_summary)")
    end
end

# Detailed printing showing all parameters
function Base.show(io::IO, ::MIME"text/plain", workflow::SMLMWorkflow)
    show(io, workflow)
end

# Pretty printing for single step
function Base.show(io::IO, step::ProcessingStep)
    print(io, "ProcessingStep: $(step.step_name) ($(step.package).$(step.function_name))")
end

"""
    summarize_output(obj) → String

Create a brief summary of an output object for workflow tracking.
Specialized methods for common SMLM types.
"""
summarize_output(obj) = string(typeof(obj))

# Specializations for common types
summarize_output(arr::AbstractArray) = "$(typeof(arr)) of size $(size(arr))"
summarize_output(nt::NamedTuple) = "NamedTuple with keys $(keys(nt))"
summarize_output(::Nothing) = "nothing"
