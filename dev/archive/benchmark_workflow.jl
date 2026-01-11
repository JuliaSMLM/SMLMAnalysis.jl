"""
Comprehensive SMLM Workflow Benchmark

Profiles the complete pipeline: Simulation → Detection → Fitting
Analyzes:
- Wall clock time per stage
- Memory allocations
- GPU memory usage and utilization
- Throughput metrics
- Bottleneck identification

Run with: julia --project=. dev/benchmark_workflow.jl
"""

using SMLMAnalysis
using SMLMData
using SMLMSim
using SMLMBoxer
using GaussMLE
using MicroscopePSFs
using Statistics
using Printf
using CUDA

# =============================================================================
# GPU Memory Monitoring Utilities
# =============================================================================

"""Wait for GPU tasks to complete and return available memory"""
function gpu_sync_and_memory()
    if CUDA.functional()
        CUDA.synchronize()
        GC.gc(false)  # Quick GC
        CUDA.reclaim()  # Reclaim GPU memory
        return CUDA.available_memory() / 1e9  # GB
    end
    return 0.0
end

"""Get GPU memory info as named tuple"""
function gpu_memory_info()
    if !CUDA.functional()
        return (available=0.0, total=0.0, used=0.0, percent_used=0.0)
    end
    CUDA.synchronize()
    total = CUDA.total_memory() / 1e9
    avail = CUDA.available_memory() / 1e9
    used = total - avail
    pct = 100 * used / total
    return (available=avail, total=total, used=used, percent_used=pct)
end

"""Wait for GPU memory to be available (blocking until threshold is met)"""
function wait_for_gpu_memory(required_gb::Float64=1.0; timeout_sec::Float64=30.0, poll_interval::Float64=0.1)
    if !CUDA.functional()
        return true
    end

    start_time = time()
    while true
        CUDA.synchronize()
        GC.gc(false)
        CUDA.reclaim()

        avail = CUDA.available_memory() / 1e9
        if avail >= required_gb
            return true
        end

        if time() - start_time > timeout_sec
            @warn "Timeout waiting for GPU memory. Available: $(round(avail, digits=2)) GB, Required: $required_gb GB"
            return false
        end

        sleep(poll_interval)
    end
end

# =============================================================================
# Benchmark Result Types
# =============================================================================

struct StageResult
    name::String
    time_sec::Float64
    alloc_mb::Float64
    gpu_mem_before_gb::Float64
    gpu_mem_after_gb::Float64
    output_count::Int
    throughput::Float64  # items/sec
end

function Base.show(io::IO, r::StageResult)
    @printf(io, "%-20s %7.3fs  %8.1f MB  %6.1f → %6.1f GB GPU  %7d items  %8.0f/s",
            r.name, r.time_sec, r.alloc_mb,
            r.gpu_mem_before_gb, r.gpu_mem_after_gb,
            r.output_count, r.throughput)
end

# =============================================================================
# Benchmark Functions
# =============================================================================

"""Benchmark a single stage with GPU sync and memory tracking"""
function benchmark_stage(name::String, f::Function; wait_for_memory::Bool=true, required_gb::Float64=2.0)
    if wait_for_memory
        wait_for_gpu_memory(required_gb)
    end

    mem_before = gpu_memory_info()

    # Run the function with timing and allocation tracking
    stats = @timed result = f()

    CUDA.functional() && CUDA.synchronize()
    mem_after = gpu_memory_info()

    # Extract output count from result
    output_count = if result isa ROIBatch
        length(result)
    elseif result isa BasicSMLD
        length(result.emitters)
    elseif result isa AbstractArray
        size(result, ndims(result))  # Last dimension is typically count
    elseif result isa Tuple
        # For simulation results
        if length(result) >= 2 && result[2] isa BasicSMLD
            length(result[2].emitters)
        else
            0
        end
    else
        0
    end

    throughput = output_count > 0 ? output_count / stats.time : 0.0

    return StageResult(
        name,
        stats.time,
        stats.bytes / 1e6,
        mem_before.used,
        mem_after.used,
        output_count,
        throughput
    ), result
end

"""Run benchmark with multiple trials and warmup"""
function run_benchmark(config; n_trials::Int=3, warmup::Bool=true)
    println("="^80)
    println("SMLM WORKFLOW BENCHMARK")
    println("="^80)

    # Print GPU info
    if CUDA.functional()
        dev = CUDA.device()
        println("\nGPU: $(CUDA.name(dev))")
        println("     $(round(CUDA.total_memory()/1e9, digits=1)) GB total memory")
    else
        println("\nNo GPU available - running on CPU")
    end

    println("\nConfiguration:")
    println("  Image size: $(config.nx) × $(config.ny) pixels")
    println("  Frames: $(config.nframes)")
    println("  Density: $(config.density) emitters/μm²")
    println("  Pixel size: $(config.pixel_size) μm")
    println("  PSF σ: $(config.psf_sigma) μm")
    println("  Box size: $(config.boxsize)")
    println("  Trials: $n_trials $(warmup ? "(+ warmup)" : "")")

    # Setup camera and simulation parameters
    camera = IdealCamera(config.nx, config.ny, config.pixel_size)

    sim_params = StaticSMLMParams(
        density = config.density,
        σ_psf = config.psf_sigma,
        nframes = config.nframes,
        framerate = 50.0,
        ndatasets = 1,
        ndims = 2
    )

    pattern = Nmer2D(n=1, d=0.0)  # Single point emitters
    molecule = GenericFluor(photons=config.photons_per_sec, k_off=20.0, k_on=0.05)
    psf = GaussianPSF(config.psf_sigma)

    # Fitter setup
    fitter = GaussMLEFitter(
        psf_model = GaussianXYNB(Float32(config.psf_sigma)),
        iterations = config.fit_iterations,
        device = config.use_gpu ? :gpu : :cpu,
        batch_size = config.fit_batch_size
    )

    all_results = Dict{String, Vector{StageResult}}()

    n_runs = warmup ? n_trials + 1 : n_trials

    for trial in 1:n_runs
        is_warmup = warmup && trial == 1
        trial_label = is_warmup ? "Warmup" : "Trial $(trial - (warmup ? 1 : 0))"

        println("\n" * "-"^60)
        println(trial_label)
        println("-"^60)

        trial_results = StageResult[]

        # Wait for GPU before starting
        wait_for_gpu_memory(4.0)

        # Stage 1: Simulation
        result_sim, sim_result = benchmark_stage("1. Simulate", () -> begin
            simulate(sim_params; pattern=pattern, molecule=molecule, camera=camera)
        end)
        _, smld_true, smld_noisy = sim_result
        push!(trial_results, result_sim)
        println("  ", result_sim)

        # Stage 2: Image Generation
        # Use support=1.0 (1μm radius) for faster PSF integration
        result_gen, images = benchmark_stage("2. Gen Images", () -> begin
            gen_images(smld_noisy, psf; bg=config.background, poisson_noise=true, threaded=true, support=1.0)
        end)
        push!(trial_results, result_gen)
        println("  ", result_gen)

        # Stage 3: Detection
        result_detect, roi_batch = benchmark_stage("3. Detection", () -> begin
            getboxes(images, camera;
                boxsize=config.boxsize,
                overlap=2.0,
                sigma_small=config.psf_sigma / config.pixel_size,
                sigma_large=2.0 * config.psf_sigma / config.pixel_size,
                minval=config.detection_threshold,
                use_gpu=config.use_gpu)
        end)
        push!(trial_results, result_detect)
        println("  ", result_detect)

        # Stage 4: Fitting
        if length(roi_batch) > 0
            result_fit, smld_fitted = benchmark_stage("4. Fitting", () -> begin
                fit(fitter, roi_batch)
            end)
            push!(trial_results, result_fit)
            println("  ", result_fit)
        else
            println("  4. Fitting: SKIPPED (no detections)")
        end

        # Only record non-warmup trials
        if !is_warmup
            for r in trial_results
                if !haskey(all_results, r.name)
                    all_results[r.name] = StageResult[]
                end
                push!(all_results[r.name], r)
            end
        end

        # Force cleanup between trials
        images = nothing
        roi_batch = nothing
        smld_fitted = nothing
        GC.gc()
        CUDA.functional() && CUDA.reclaim()
    end

    return all_results, config
end

"""Analyze benchmark results and print summary"""
function analyze_results(all_results, config)
    println("\n" * "="^80)
    println("BENCHMARK SUMMARY")
    println("="^80)

    # Compute statistics for each stage
    println("\nStage Statistics ($(length(first(values(all_results)))) trials):")
    println("-"^80)
    @printf("%-20s %10s %10s %12s %12s\n", "Stage", "Mean Time", "Std Dev", "Mean Alloc", "Mean Thput")
    println("-"^80)

    total_time = 0.0
    total_alloc = 0.0
    bottleneck_stage = ""
    bottleneck_time = 0.0

    for stage_name in sort(collect(keys(all_results)))
        results = all_results[stage_name]
        times = [r.time_sec for r in results]
        allocs = [r.alloc_mb for r in results]
        thputs = [r.throughput for r in results]

        mean_time = mean(times)
        std_time = std(times)
        mean_alloc = mean(allocs)
        mean_thput = mean(thputs)

        total_time += mean_time
        total_alloc += mean_alloc

        if mean_time > bottleneck_time
            bottleneck_time = mean_time
            bottleneck_stage = stage_name
        end

        @printf("%-20s %9.3fs %9.3fs %10.1f MB %10.0f/s\n",
                stage_name, mean_time, std_time, mean_alloc, mean_thput)
    end

    println("-"^80)
    @printf("%-20s %9.3fs            %10.1f MB\n", "TOTAL", total_time, total_alloc)

    # Bottleneck Analysis
    println("\n" * "="^80)
    println("BOTTLENECK ANALYSIS")
    println("="^80)

    for stage_name in sort(collect(keys(all_results)))
        results = all_results[stage_name]
        mean_time = mean([r.time_sec for r in results])
        pct = 100 * mean_time / total_time
        bar_len = round(Int, pct / 2)
        bar = "█"^bar_len * "░"^(50 - bar_len)
        @printf("%-20s %5.1f%% %s\n", stage_name, pct, bar)
    end

    println("\nPrimary Bottleneck: $bottleneck_stage ($(round(100*bottleneck_time/total_time, digits=1))% of total time)")

    # Throughput Analysis
    println("\n" * "="^80)
    println("THROUGHPUT ANALYSIS")
    println("="^80)

    # Expected vs actual counts
    fov_area = config.nx * config.pixel_size * config.ny * config.pixel_size  # μm²
    expected_patterns_per_frame = config.density * fov_area
    expected_emitters_total = expected_patterns_per_frame * config.nframes * 0.3  # ~30% on per frame

    if haskey(all_results, "3. Detection")
        detect_results = all_results["3. Detection"]
        mean_detections = mean([r.output_count for r in detect_results])
        detection_rate = mean_detections / expected_emitters_total
        @printf("  Expected emitter events: ~%.0f\n", expected_emitters_total)
        @printf("  Mean detections:         %.0f (%.1f%% detection rate)\n",
                mean_detections, 100*detection_rate)
    end

    if haskey(all_results, "4. Fitting")
        fit_results = all_results["4. Fitting"]
        mean_fits = mean([r.output_count for r in fit_results])
        mean_thput = mean([r.throughput for r in fit_results])
        @printf("  Mean fits:               %.0f\n", mean_fits)
        @printf("  Fitting throughput:      %.0f fits/s\n", mean_thput)
    end

    # Per-frame analysis
    frames_per_sec = config.nframes / total_time
    @printf("\n  Overall pipeline: %.1f frames/sec\n", frames_per_sec)

    return (bottleneck_stage=bottleneck_stage, total_time=total_time, total_alloc=total_alloc)
end

"""Generate performance recommendations based on benchmark results"""
function generate_recommendations(all_results, config, analysis)
    println("\n" * "="^80)
    println("PERFORMANCE RECOMMENDATIONS")
    println("="^80)

    recommendations = String[]

    # Analyze image generation
    if haskey(all_results, "2. Gen Images")
        gen_results = all_results["2. Gen Images"]
        mean_time = mean([r.time_sec for r in gen_results])
        mean_alloc = mean([r.alloc_mb for r in gen_results])

        if mean_time > 1.0
            push!(recommendations, """
1. IMAGE GENERATION BOTTLENECK ($(round(mean_time, digits=2))s)
   - Current: Sequential pixel integration with PSF sampling
   - Issue: gen_images loops over frames and emitters on CPU

   RECOMMENDATIONS:
   a) Use PSF support parameter to limit integration radius:
      gen_images(smld, psf; support=1.0)  # 1μm radius instead of Inf

   b) Consider GPU-accelerated PSF rendering:
      - Pre-compute PSF lookup table
      - Use texture memory for PSF interpolation
      - Parallelize across emitters and pixels

   c) For simulation benchmarks, pre-generate images once and reuse

   EXPECTED SPEEDUP: 5-20x with finite support, 50-100x with GPU rendering
""")
        end
    end

    # Analyze detection
    if haskey(all_results, "3. Detection")
        detect_results = all_results["3. Detection"]
        mean_time = mean([r.time_sec for r in detect_results])
        mean_alloc = mean([r.alloc_mb for r in detect_results])
        mean_rois = mean([r.output_count for r in detect_results])

        # Memory per frame
        frame_size_mb = config.nx * config.ny * 4 / 1e6  # Float32
        total_stack_mb = frame_size_mb * config.nframes
        alloc_ratio = mean_alloc / total_stack_mb

        if alloc_ratio > 5
            push!(recommendations, """
2. DETECTION MEMORY OVERHEAD ($(round(mean_alloc, digits=0)) MB for $(round(total_stack_mb, digits=0)) MB input)
   - Allocation ratio: $(round(alloc_ratio, digits=1))x input size
   - This suggests excessive intermediate allocations

   RECOMMENDATIONS:
   a) Process frames in batches to reduce peak memory:
      - Current: Full stack convolution
      - Better: Stream frames through GPU in smaller batches

   b) Reuse convolution buffers:
      - Pre-allocate DoG filter outputs
      - Use in-place operations where possible

   c) For GPU: Use pinned memory for host-device transfers

   EXPECTED IMPROVEMENT: 2-4x memory reduction
""")
        end

        if mean_time > 0.5 && config.use_gpu
            thput = mean_rois / mean_time
            push!(recommendations, """
3. DETECTION THROUGHPUT ($(round(thput, digits=0)) ROIs/s)
   - Current implementation uses cuDNN for convolution
   - Local max finding may have suboptimal GPU utilization

   RECOMMENDATIONS:
   a) Fuse DoG filter and local max into single kernel:
      - Avoid intermediate array allocations
      - Reduce GPU memory bandwidth

   b) Use shared memory for local max neighborhood:
      - Current: Multiple global memory reads per pixel
      - Better: Load tile to shared memory once

   c) Async overlap between batches:
      - Pipeline: batch N detection while batch N-1 fitting runs

   EXPECTED SPEEDUP: 2-5x with kernel fusion
""")
        end
    end

    # Analyze fitting
    if haskey(all_results, "4. Fitting")
        fit_results = all_results["4. Fitting"]
        mean_time = mean([r.time_sec for r in fit_results])
        mean_fits = mean([r.output_count for r in fit_results])
        mean_thput = mean([r.throughput for r in fit_results])

        # Theoretical max throughput estimate
        # Modern GPU can do ~1M simple operations per ROI
        # With 11x11 ROI, 20 iterations, ~6 params: ~25k flops/fit
        # RTX 3090 theoretical: 35 TFLOPS → ~1.4M fits/s
        # Practical limit (memory bound): ~500k fits/s
        theoretical_max = 500_000
        efficiency = mean_thput / theoretical_max

        if efficiency < 0.1
            push!(recommendations, """
4. FITTING EFFICIENCY ($(round(100*efficiency, digits=1))% of theoretical max)
   - Current: $(round(mean_thput, digits=0)) fits/s
   - Theoretical GPU limit: ~500k fits/s

   BOTTLENECK ANALYSIS:
   a) Memory transfer overhead:
      - Each batch copies data to GPU and results back
      - With batch_size=$(config.fit_batch_size), this is $(round(mean_fits/config.fit_batch_size, digits=0)) transfers

   b) Kernel launch overhead:
      - Small batches = many kernel launches
      - Consider larger batch sizes if memory allows

   c) Fisher matrix inversion per-ROI:
      - Cholesky decomposition is sequential per-thread
      - Consider batched LAPACK for CPU path

   RECOMMENDATIONS:
   a) Increase batch size to reduce transfer overhead:
      fitter = GaussMLEFitter(batch_size=50_000)  # or 100_000

   b) Use pinned memory for async transfers:
      - Pre-allocate pinned buffers
      - Overlap transfer with computation

   c) Consider streaming multiple batches:
      - While batch N computes, transfer batch N+1

   d) For large datasets, keep data on GPU:
      - Avoid round-trip for detection → fitting
      - Allocate ROI extraction on GPU

   EXPECTED SPEEDUP: 2-10x with async streaming
""")
        end

        # Check for p-value computation overhead
        if mean_fits > 0
            pvalue_overhead = mean_fits * 0.00001  # Rough estimate: chi2 CDF is ~10μs per call
            if pvalue_overhead / mean_time > 0.1
                push!(recommendations, """
5. P-VALUE COMPUTATION OVERHEAD
   - Chi-squared CDF computed per fit on CPU
   - For $(round(mean_fits, digits=0)) fits, this adds ~$(round(pvalue_overhead*1000, digits=0))ms

   RECOMMENDATIONS:
   a) Batch p-value computation with vectorized operations
   b) Use lookup table for common chi2 values
   c) Compute p-values lazily (only when filtering)

   EXPECTED SPEEDUP: Minor (5-10% of fitting time)
""")
            end
        end
    end

    # Pipeline-level recommendations
    push!(recommendations, """
6. PIPELINE OPTIMIZATION OPPORTUNITIES

   a) END-TO-END GPU PIPELINE:
      - Current: Data moves CPU ↔ GPU at each stage
      - Better: Keep data on GPU throughout detection → fitting

      Implementation:
      ```julia
      # Allocate GPU buffers once
      d_images = CuArray(images)
      d_rois = preallocate_rois(n_expected, boxsize)

      # Detection outputs directly to GPU ROI buffer
      n_rois = detect_to_gpu!(d_rois, d_images, params)

      # Fitting reads from GPU ROI buffer
      smld = fit_from_gpu(fitter, d_rois, n_rois)
      ```

   b) ASYNC PROCESSING:
      - Stream frames through pipeline
      - Frame N detection while Frame N-1 fits
      - Requires careful synchronization

   c) MEMORY POOLING:
      - Pre-allocate all GPU arrays at startup
      - Reuse buffers across pipeline stages
      - Eliminates allocation overhead

   d) REDUCE PRECISION WHERE POSSIBLE:
      - Detection: Float16 sufficient for DoG
      - PSF lookup: Int8 or Float16 texture
      - Only fitting needs Float32 precision
""")

    # Print all recommendations
    for rec in recommendations
        println(rec)
    end

    return recommendations
end

# =============================================================================
# Main Benchmark Execution
# =============================================================================

function main()
    # Default benchmark configuration
    config = (
        # Image dimensions
        nx = 256,
        ny = 256,
        nframes = 100,

        # Physics
        pixel_size = 0.1f0,      # μm
        psf_sigma = 0.13f0,      # μm
        density = 2.0,           # emitters/μm² (higher for more detections)
        photons_per_sec = 15000.0,  # Higher photon rate
        background = 5.0,        # Lower background

        # Detection
        boxsize = 11,
        detection_threshold = 5.0,  # Lower threshold for actual detections

        # Fitting
        fit_iterations = 20,
        fit_batch_size = 10_000,

        # Hardware
        use_gpu = CUDA.functional()
    )

    # Run benchmark
    all_results, config = run_benchmark(config; n_trials=3, warmup=true)

    # Analyze results
    analysis = analyze_results(all_results, config)

    # Generate recommendations
    generate_recommendations(all_results, config, analysis)

    println("\n" * "="^80)
    println("Benchmark complete")
    println("="^80)

    return all_results, config, analysis
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
