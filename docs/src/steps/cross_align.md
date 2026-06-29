```@meta
CurrentModule = SMLMAnalysis
```

# Cross-Alignment

In a multi-color experiment each channel is acquired separately, so a residual
shift between channels — from chromatic offset, stage settling, or registration
error — leaves the colors misregistered in the overlay. The cross-alignment step
estimates that shift and applies it, **registering every channel to a common
reference frame** so structures that should coincide actually overlap. It is a
**multi-channel** step: it operates on a `Vector{BasicSMLD}` (one SMLD per
channel), is selected by a `CrossAlignConfig`, and is a native SMLMAnalysis step
built on [SMLMDriftCorrection](https://github.com/JuliaSMLM/SMLMDriftCorrection.jl)'s
`align_smld`.

```julia
analyze(smlds, CrossAlignConfig(method = :entropy))   # → (aligned_smlds, StepInfo)
```

## When to use / prerequisites

- Use in a multi-channel pipeline, after each channel has been independently
  detected, fitted, filtered, and (ideally) drift-corrected. It is normally
  listed in a [`MultiTargetConfig`](@ref)`.steps` and dispatched on the
  resulting `Vector{BasicSMLD}` — see [Multi-Channel](@ref
  "Multi-Channel Analysis").
- Alignment is data-driven: it registers channels by the structures they share,
  so it needs **spatially correlated content** across channels (overlapping or
  co-localized features). It corrects a global translation, not chromatic
  distortion or rotation.
- A common pattern is to bracket this step between two [Composite Render](@ref)
  steps to see the overlay before and after alignment.

## Inputs, returns & artifacts

- **Input:** the channels as a `Vector{BasicSMLD}`.
- **Returns:** `(aligned_smlds, StepInfo)`. The returned vector holds the same
  channels with the per-channel shift applied to their coordinates; the step is
  **state-modifying** for the downstream multi-target steps. The per-channel
  shifts and the upstream alignment info live on `StepInfo.info` (a
  [`CrossAlignInfo`](@ref)).
- **Artifacts** (when `outdir` is set), written under `NN_crossalign/`:
  `config.toml`, `info.toml` (the upstream `AlignInfo` fields), and — at
  `STANDARD` verbosity or above — `stats.md` with the method, channel count, max
  shift, and a per-channel table of X/Y/magnitude shifts in nm.

## Concept

Aligning channels and correcting drift are the same registration problem with a
shorter trajectory: here the "trajectory" is a single rigid translation per
channel rather than a time-varying path. With `method = :entropy` the step seeds
the offset with a histogram cross-correlation and then refines it by
**minimizing the entropy** of the combined localization cloud — the shift that
makes shared structures pile up most tightly (Cnossen et al. 2021), the same
fiducial-free, redundancy-based idea used for [Drift Correction](@ref) (Wester
et al. 2021). With `method = :fft` it stops at the cross-correlation peak (no
entropy refinement), which is faster but less precise. For the algorithm in full
see the SMLMDriftCorrection documentation.

## Configuration

`CrossAlignConfig` is an SMLMAnalysis multi-target config whose fields map onto
the upstream `AlignConfig`:

| field | typical/default | meaning |
|---|---|---|
| `method` | `:entropy` | `:entropy` (cross-correlation seed + entropy refinement) or `:fft` (cross-correlation only) |
| `maxn` | `100` | maximum neighbors used in the entropy calculation |
| `histbinsize` | `0.05` | histogram bin size (µm) for the cross-correlation seed |

```julia
# Two-channel overlay: render, align, render again
mt = MultiTargetConfig(
    labels = [:ch1, :ch2],
    steps = [
        CompositeRenderConfig(zoom = 20.0, strategy = GaussianRender()),
        CrossAlignConfig(method = :entropy, histbinsize = 0.05),
        CompositeRenderConfig(zoom = 20.0, strategy = GaussianRender()),  # post-alignment
    ],
    outdir = "output/cell1/",
)

# Or call the step directly on a channel vector
(aligned, info) = analyze(smlds, CrossAlignConfig(method = :entropy))
shifts = info.info.shifts          # Vector{Vector{Float64}}, one (x, y) per channel, in µm
```

## Output & interpretation

The step's `StepInfo.summary` reports the headline numbers:

| field | meaning |
|---|---|
| `max_shift_nm` | largest per-channel shift magnitude applied (nm) |
| `n_channels` | number of channels aligned |
| `method` | the alignment method used (`:entropy` / `:fft`) |

`StepInfo.info` (a [`CrossAlignInfo`](@ref)) carries the full detail:
`shifts` (the applied per-channel `(x, y)` offsets, **in µm**), `max_shift_nm`,
the upstream `align_info`, and `elapsed_s`.

Sanity checks: shifts of a few tens to a couple hundred nm are typical for
chromatic/registration offset, and the post-alignment composite render should
show the colors snapping into register. A shift much larger than the expected
inter-channel offset usually means the channels lack enough shared structure to
register against.

## Notes & caveats

- **Translation only.** Cross-alignment removes a global X/Y shift per channel;
  it does not correct rotation, scaling, or field-dependent chromatic
  distortion. For those, apply an upstream geometric transform first.
- **Shared structure is required.** Channels with no co-localized or overlapping
  features give an ill-defined offset; the cross-correlation seed will lock onto
  noise.
- **`:fft` vs `:entropy`.** Use `:fft` for a fast first pass; use `:entropy`
  (the default) when you need the tighter, sub-bin offset.
- **Units.** `shifts` are in µm (matching the SMLD coordinates); the summary and
  `stats.md` report magnitudes in nm.

## References

- J. Cnossen, T. J. Cui, C. Joo, C. Smith. "Drift correction in localization
  microscopy using entropy minimization." *Optics Express* **29**, 27961 (2021).
  [doi:10.1364/OE.426620](https://doi.org/10.1364/OE.426620)
- M. J. Wester, *et al.* "Robust, fiducial-free drift correction for
  super-resolution imaging." *Scientific Reports* **11**, 23672 (2021).
  [doi:10.1038/s41598-021-02850-7](https://doi.org/10.1038/s41598-021-02850-7)

See the [SMLMDriftCorrection documentation](https://github.com/JuliaSMLM/SMLMDriftCorrection.jl)
for the alignment algorithm in full and all configuration options.
