"""
Run remaining analyses:
1. 1ch sxsy model (sigma already done)
2. 2ch sxsy with polarization analysis
"""

# Modify run_modes for 1ch to skip sigma (already completed)
# Then run the standard workflow

using SMLMAnalysis
using SMLMDriftCorrection
using CairoMakie
using Statistics
using Printf
using Dates
using FileIO
using Images: RGB

# Files
file_1ch = "/mnt/nas/adapt/projects/smart-microscope/data/DNA paint ruler/2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-1ch--2025-10-23_11-29-25.h5"
file_2ch = "/mnt/nas/adapt/projects/smart-microscope/data/DNA paint ruler/2025-10-23/20R-ruler-0.1exp-TIRF-onlyZFocusLockDT1-2ch--2025-10-23_12-04-53.h5"

# Run 1ch sxsy only
println("="^80)
println("Running 1ch sxsy analysis...")
println("="^80)

# Set ARGS and include the workflow with modified run_modes
# Actually, let's just run the workflow directly for 2ch first since that's more critical

println("\n\n")
println("="^80)
println("Running 2ch sxsy analysis with polarization...")
println("="^80)

# Run 2ch
empty!(ARGS)
push!(ARGS, file_2ch)
include("workflow_standard.jl")

println("\n\n")
println("="^80)
println("Running 1ch sxsy analysis...")
println("="^80)

# For 1ch, we need to modify run_modes - do this by running workflow_1ch_sxsy_only.jl
# Or just accept that sigma will re-run (it's fast)
empty!(ARGS)
push!(ARGS, file_1ch)
include("workflow_standard.jl")
