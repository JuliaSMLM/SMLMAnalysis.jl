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
    SMLMAnalysis is being prepared for registration. Until it and its
    JuliaSMLM dependencies are registered, install from source as below. (This
    page will switch to `Pkg.add` once the registration is live.)

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

## Optional extensions

Two capabilities load only when their packages are present, as Julia package
**extensions** (weak dependencies) — you do not pay their load cost otherwise:

| Extension | Activated by | Adds |
|-----------|--------------|------|
| `SMLMAnalysisPSFLearningExt` | `PSFLearning` | PSF-learning step |
| `SMLMAnalysisDeepFitExt` | `SMLMDeepFit` | deep-learning detection/fitting steps |

Install the corresponding package alongside SMLMAnalysis to enable the extension;
the relevant step configs become available automatically.

## Verifying the install

A quick end-to-end check on simulated data:

```julia
using SMLMAnalysis

# Simulate a small dataset and its image stack (SMLMSim, re-exported)
# then run a minimal pipeline:
config = AnalysisConfig(
    camera = cam,
    steps = [
        DetectFitConfig(boxer = BoxerConfig(boxsize = 9, psf_sigma = 0.130)),
        RenderConfig(zoom = 10),
    ],
)
(result, info) = analyze(images, config)
@show length(result.smld.emitters)
```

See [Getting Started](@ref Tutorial) for a full simulated walkthrough, and
[Running a Pipeline](@ref) for the config-driven and step-by-step styles.

## Building the documentation

The documentation is built from the **main checkout** of the repository (its
`docs/` environment anchors its path-sourced JuliaSMLM dependencies relative to
that checkout):

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```
