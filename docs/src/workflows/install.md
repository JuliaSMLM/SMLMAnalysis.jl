```@meta
CurrentModule = SMLMAnalysis
```

# Installation & Setup

## Requirements

- **Julia 1.12** or newer.
- A **CUDA-capable GPU** is recommended. The [fitting](@ref "Detection & Fitting")
  step uses GaussMLE, which is GPU-accelerated via CUDA; see the GaussMLE
  documentation for GPU/CPU options.
- Start Julia with **multiple threads** — several steps (drift correction, frame
  connection) are threaded:

  ```bash
  julia -t auto --project=.
  ```

## Installing

Once SMLMAnalysis is registered in the Julia General registry, installation will
be the usual one-liner:

```julia
using Pkg
Pkg.add("SMLMAnalysis")
```

!!! note "Registration status"
    SMLMAnalysis is being prepared for registration; every JuliaSMLM dependency
    is already in the General registry. Until SMLMAnalysis itself is registered,
    install from source as below. (This page will switch to `Pkg.add` once the
    registration is live.)

### From source (current)

Clone the repository and instantiate its environment:

```bash
git clone https://github.com/JuliaSMLM/SMLMAnalysis.jl
cd SMLMAnalysis.jl
julia -t auto --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then:

```julia
using SMLMAnalysis
```

`using SMLMAnalysis` re-exports the key ecosystem types and verbs (cameras,
emitter types, the step configs, `cluster`, `render`, `run_bagol`,
`frameconnect`, …), so for most work no further imports are needed.

## Verifying the install

A quick end-to-end check on simulated data (runs as-is; the simulation verbs are
re-exported from SMLMSim):

```julia
using SMLMAnalysis

# Simulate a small 8-mer sample and synthesize its camera frames
cam = IdealCamera(64, 64, 0.1)                          # 64×64 px, 100 nm pixels
sim = StaticSMLMConfig(density = 2.0, σ_psf = 0.13, nframes = 200, ndatasets = 1)
(_, si) = simulate(sim; pattern  = Nmer2D(n = 8, d = 0.05),
                        molecule = GenericFluor(photons = 5.0e4, k_off = 20.0, k_on = 0.05),
                        camera   = cam)
(images, _) = gen_images(si.smld_model, SMLMAnalysis.MicroscopePSFs.GaussianPSF(0.13);
                         dataset = 1, bg = 20.0, poisson_noise = true)

# Run a minimal pipeline on it
config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(boxer = BoxerConfig(boxsize = 7, psf_sigma = 0.13)),
        RenderConfig(zoom = 10),
    ],
)
(result, info) = analyze(images, config)
@show length(result.smld.emitters)
```

`examples/loading_data.jl` in the repository is the runnable version of this check.

See [Getting Started](@ref Tutorial) for a full simulated walkthrough, and
[Running a Pipeline](@ref) for the config-driven and step-by-step styles.

## Building the documentation

The documentation builds from a checkout of the repository. The `docs/`
environment sources SMLMAnalysis from that checkout itself (the `[sources]` entry
in `docs/Project.toml`, honoured on Julia ≥ 1.11), so no `Pkg.develop` step is
needed — and running one would rewrite the tracked `docs/Project.toml` with an
absolute path:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```
