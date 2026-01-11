# SMLMAnalysis Performance Report

**Date:** December 2025
**Hardware:** NVIDIA RTX 4090 (25.3 GB), Linux
**Packages:** SMLMAnalysis, SMLMBoxer, GaussMLE, SMLMSim

## Executive Summary

The SMLMAnalysis pipeline achieves **excellent fitting performance** (~900k fits/s on GPU) but is limited by **image generation** (CPU-bound) and **detection memory scaling** for large datasets.

| Stage | Performance | Status |
|-------|-------------|--------|
| Simulation | 80k emitters/s | Excellent |
| Image Generation | 20-40 frames/s | **BOTTLENECK** |
| Detection (small) | 2500 fps | Excellent |
| Detection (large) | 180 fps | Good |
| Fitting | 800k-960k fits/s | Excellent |

## Benchmark Results

### Detection Performance (GPU)

| Size | Frames | Time | Throughput | Allocations |
|------|--------|------|------------|-------------|
| 256×256×100 | 100 | 0.032s | 2506 fps | 38 MB |
| 512×512×500 | 500 | 0.493s | 961 fps | 588 MB |
| 1024×1024×1000 | 1000 | 5.47s | **181 fps** | 8533 MB |

**Observation:** Large images (>4 GB) trigger memory batching, reducing throughput by 5x.

### Fitting Performance (GPU)

| ROIs | Time | Throughput | Allocations |
|------|------|------------|-------------|
| 1,000 | 0.008s | 129k/s | 0.7 MB |
| 10,000 | 0.011s | 936k/s | 6.7 MB |
| 50,000 | 0.052s | 963k/s | 33 MB |
| 100,000 | 0.121s | 825k/s | 67 MB |
| 200,000 | 0.255s | 785k/s | 133 MB |
| 500,000 | 0.577s | 867k/s | 332 MB |

**Peak performance:** ~960k fits/s at optimal batch sizes (10k-50k ROIs).

### Full Pipeline (Simulation → Detection → Fitting)

| Stage | Mean Time | % Total | Allocs |
|-------|-----------|---------|--------|
| Simulate | 0.005s | 0.2% | 4.6 MB |
| Gen Images | 2.56s | **97.6%** | 1248 MB |
| Detection | 0.054s | 2.1% | 55 MB |
| Fitting | 0.005s | 0.2% | 0.2 MB |

**Primary Bottleneck:** Image generation is 98% of pipeline time.

## Bottleneck Analysis

### 1. Image Generation (CRITICAL)

**Current implementation:** `SMLMSim.gen_images` uses CPU-based pixel integration with PSF sampling. Each emitter requires integrating over all pixels within its support region.

**Root cause:**
- Sequential frame-by-frame processing
- PSF integration computed per-pixel, per-emitter
- No GPU acceleration

**Impact:** 2.5s for 100 frames of 256×256 (even with `support=1.0`)

### 2. Detection Memory Scaling

**Current implementation:** SMLMBoxer batches frames when total memory exceeds GPU capacity.

**Root cause:**
- `n_copies = 6` memory multiplier for safety
- Each batch requires full round-trip CPU↔GPU
- Batch coordination overhead

**Impact:** 5x slowdown when image stack exceeds GPU memory

### 3. CPU↔GPU Data Transfers

**Current implementation:** Each stage operates independently:
1. Detection: CPU images → GPU → CPU ROIs
2. Fitting: CPU ROIs → GPU → CPU results

**Root cause:**
- No shared GPU memory pool
- Defensive copying for safety

**Impact:** Minor for small datasets, significant at scale

## Recommendations

### HIGH PRIORITY

#### 1. GPU-Accelerated Image Generation

**Current:** ~40 fps (CPU)
**Target:** 1000+ fps (GPU)
**Expected speedup:** 25-50x

```julia
# Proposed approach: GPU texture-based PSF rendering
function gen_images_gpu(smld, psf; kwargs...)
    # Pre-compute PSF lookup table on GPU
    d_psf_lut = precompute_psf_texture(psf, sampling=4)

    # Allocate output on GPU
    d_images = CUDA.zeros(Float32, height, width, n_frames)

    # Render all emitters in parallel
    # Each thread handles one emitter, writes to shared image
    render_emitters_kernel!(d_images, emitters, d_psf_lut)

    return Array(d_images)
end
```

**Implementation notes:**
- Use texture memory for PSF lookup (hardware interpolation)
- Atomic adds for overlapping emitters
- Process all frames in single kernel launch

#### 2. End-to-End GPU Pipeline

**Current:** Multiple CPU↔GPU round trips
**Target:** Data stays on GPU throughout

```julia
struct GPUPipeline
    d_images::CuArray{Float32,3}
    d_rois::CuArray{Float32,3}
    d_results::CuArray{Float32,2}
end

function process_gpu!(pipeline, fitter)
    # Detection directly to GPU ROI buffer
    n_rois = detect_to_gpu!(pipeline.d_rois, pipeline.d_images)

    # Fitting from GPU ROI buffer
    fit_from_gpu!(pipeline.d_results, fitter, pipeline.d_rois, n_rois)

    # Single transfer at end
    return Array(@view pipeline.d_results[:, 1:n_rois])
end
```

### MEDIUM PRIORITY

#### 3. Detection Streaming

**Current:** Full stack loaded, then batched
**Target:** Stream frames through GPU

```julia
function getboxes_streaming(images, camera; chunk_size=100, kwargs...)
    # Process in overlapping chunks
    for chunk_start in 1:chunk_size:n_frames
        chunk = @view images[:, :, chunk_start:chunk_end]

        # Async transfer next chunk while processing current
        @async transfer_next_chunk!(d_buffer_next, images, next_chunk)

        # Process current chunk
        roi_batch_chunk = detect_gpu!(d_buffer, chunk)
        append!(all_rois, roi_batch_chunk)
    end
end
```

#### 4. Memory Pool Reuse

**Current:** Allocate/free per operation
**Target:** Pre-allocated buffer pool

```julia
struct MemoryPool
    buffers::Dict{Tuple{Type,Dims}, Vector{CuArray}}
end

function allocate!(pool, T, dims)
    key = (T, dims)
    if !isempty(pool.buffers[key])
        return pop!(pool.buffers[key])
    end
    return CuArray{T}(undef, dims)
end

function release!(pool, arr)
    key = (eltype(arr), size(arr))
    push!(pool.buffers[key], arr)
end
```

### LOWER PRIORITY

#### 5. Precision Reduction

Detection can use Float16 for initial filtering:
- DoG convolution in Float16 → 2x memory bandwidth
- Local max finding in Float16
- Only convert to Float32 for ROI extraction

#### 6. P-value Computation Optimization

Current: Per-fit chi-squared CDF on CPU
Target: Batched or lookup-table based

```julia
# Lookup table approach
const PVALUE_LUT = precompute_chi2_cdf(df_range, precision=1e-4)

function pvalues_lookup(log_likelihoods, df)
    χ² = -2 .* log_likelihoods
    return lookup.(Ref(PVALUE_LUT), χ², df)
end
```

## Implementation Roadmap

### Phase 1: Quick Wins (1-2 days)
- [ ] Add `support` parameter guidance to workflow documentation
- [ ] Increase default `batch_size` in GaussMLEFitter to 50,000
- [ ] Add memory estimation before processing

### Phase 2: GPU Image Generation (1-2 weeks)
- [ ] Implement PSF texture lookup
- [ ] Create GPU render kernel
- [ ] Benchmark against CPU version
- [ ] Integrate into `gen_images`

### Phase 3: Pipeline Integration (1 week)
- [ ] Create `GPUPipeline` struct
- [ ] Implement shared memory buffers
- [ ] Add streaming detection option

### Phase 4: Polish (ongoing)
- [ ] Memory pool system
- [ ] Async chunk processing
- [ ] Float16 detection path

## Appendix: Profiling Commands

```julia
# Profile detection
using Profile
@profile roi_batch = getboxes(images, camera; use_gpu=true)
Profile.print(sortby=:overhead)

# GPU profiling
CUDA.@profile roi_batch = getboxes(images, camera; use_gpu=true)

# Memory tracking
using CUDA
CUDA.memory_status()
```
