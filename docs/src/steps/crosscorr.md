```@meta
CurrentModule = SMLMAnalysis
```

# Cross-Correlation

Two proteins imaged in separate channels may be distributed independently, may
co-cluster, or may exclude one another. The cross-correlation step quantifies
this by computing the **pair cross-correlation function** ``g(r)`` between two
channels: the density of channel-B localizations at distance ``r`` from a
channel-A localization, normalized by what complete spatial randomness would
give. It is selected by a `CrossCorrConfig` and is a **native SMLMAnalysis**
step built on the KD-tree range queries of
[NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl).

```julia
analyze(smlds, CrossCorrConfig(r_max = 0.5))   # → (smlds, StepInfo)
```

This is a multi-channel step: it operates on a `Vector{BasicSMLD}` (see
[Multi-Channel](@ref "Multi-Channel Analysis")). It is **read-only** — it
computes a curve and writes it out, but the SMLDs pass through unmodified.

## When to use / prerequisites

- Run on **two channels** of localizations once both are reconstructed (after
  [detection/fitting](@ref "Detection & Fitting") and the usual
  [quality filter](@ref "Quality Filter") on each channel).
- The two channels must live in a **common coordinate frame**. If they came
  from different optical paths, register them first with
  [Cross-Alignment](@ref); otherwise a real co-localization can be masked by a
  residual channel-to-channel offset.
- Both channels need enough localizations to populate the radial bins — sparse
  data gives a noisy ``g(r)``, especially at small ``r``.

## Inputs, returns & artifacts

- **Input:** the channel vector `smlds`. The two channels are picked by the
  1-based `channels` index tuple; they must be distinct and in range (the step
  errors otherwise). The FOV area used for normalization is taken from channel
  A's camera, which assumes both channels share a camera.
- **Returns:** `(smlds, StepInfo)` — `smlds` is returned unchanged. The computed
  curve lives on `StepInfo.info`, a [`CrossCorrInfo`](@ref) carrying `r`, `g`,
  `n_a`, `n_b`, `area`, and the channel labels.
- **Artifacts** (when `outdir` is set): `crosscorr_gr.csv` (the `r,g` table,
  always written), and at `STANDARD` verbosity `crosscorr_gr.png` (the curve with
  the CSR=1 line and the peak annotated) and `stats.md` (counts, FOV area, peak
  ``g(r)`` and its radius). The config is also dumped for provenance. Nothing is
  checkpointed — the SMLDs are unchanged.

## Concept

For each localization in channel A, the step finds every channel-B localization
within `r_max` via a KD-tree `inrange` query, bins the pair distances into
shells of width `dr`, and normalizes each shell by the count expected under
complete spatial randomness (CSR):

```math
g(r) = \frac{\text{observed B–A pairs in shell } r}
            {n_A \,\rho_B \,(\pi(r_\text{outer}^2 - r_\text{inner}^2))}
```

where ``\rho_B = n_B / \text{area}`` is the mean channel-B density. The
denominator is the number of pairs a uniform, independent channel B would
contribute to that annulus. Reading the result:

- **``g(r) > 1``** — co-clustering: more B near A than chance, at separation ``r``.
- **``g(r) \approx 1``** — the two channels are independent (the CSR baseline).
- **``g(r) < 1``** — exclusion / anti-correlation: the species avoid each other.

Because this is a *cross*-correlation between distinct channels, there is no
self-pair spike at ``r \to 0`` from repeated blinks of one molecule (zero-distance
pairs are skipped); a small-``r`` rise is genuine co-localization blurred by
localization precision.

**Edge correction.** A localization near the FOV edge sees a clipped annulus, so
naive counts undercount pairs at large ``r``. With `edge_correction = true` each
pair is weighted by the inverse of the fraction of its circle of radius ``r``
that lies inside the rectangular FOV — Ripley's isotropic correction (Ripley
1977) — up-weighting near-boundary points to compensate.

## Configuration

| field | typical / default | meaning |
|---|---|---|
| `r_max` | `1.0` | maximum separation in μm over which ``g(r)`` is computed |
| `dr` | `0.01` | radial bin width in μm; sets the curve's resolution |
| `edge_correction` | `true` | apply Ripley's isotropic edge correction for the rectangular FOV |
| `channels` | `(1, 2)` | 1-based indices `(A, B)` of the two channels to correlate; must differ and be in range |

```julia
# g(r) between channels 1 and 2 out to 500 nm, 10 nm bins
(smlds, info) = analyze(smlds,
    CrossCorrConfig(r_max = 0.5, dr = 0.01, channels = (1, 2)))

gr = info.info          # CrossCorrInfo
gr.r                    # bin centers (µm)
gr.g                    # g(r) values
```

## Output & interpretation

The step's `StepInfo.summary` reports the headline numbers:

| field | meaning |
|---|---|
| `n_a`, `n_b` | localization counts in the two channels |
| `n_bins` | number of radial bins (`≈ r_max / dr`) |
| `peak_g` | the maximum of ``g(r)`` — the strongest correlation observed |
| `channel_a`, `channel_b` | the channel labels used (from the pipeline's `labels`, else `Ch1`, `Ch2`) |

Sanity checks: a `peak_g` near 1 with a flat curve means the channels are
spatially independent; a `peak_g` well above 1 at small ``r`` indicates
co-clustering, and the radius at which ``g(r)`` decays back to ~1 estimates the
co-cluster length scale. Confirm the curve does **not** diverge as ``r``
approaches `r_max` — runaway tails usually mean the edge correction is off or the
FOV/area is mis-set. Very small ``r`` bins are intrinsically noisy (few pairs);
shrink `dr` only when both channels are dense.

## Notes & caveats

- **Register channels first.** ``g(r)`` measures separation, so an uncorrected
  inter-channel offset shifts and smears the peak. Run [Cross-Alignment](@ref)
  upstream.
- **It does not relabel emitters.** This is an analysis-only step; downstream
  steps see the same SMLDs. Use it alongside [Composite Render](@ref) to *see*
  the overlap you are quantifying.
- **Shared camera assumed.** The normalization area comes from channel A's
  camera; mismatched cameras or cropped channels will bias ``\rho_B`` and hence
  the baseline.
- **Asymmetry.** Pairs are counted from A to B; for well-sampled channels
  ``g_{AB}(r)`` and ``g_{BA}(r)`` agree, but with very different counts the
  finite-sampling noise differs — check both if in doubt by swapping `channels`.

## References

- B. D. Ripley. "Modelling spatial patterns." *Journal of the Royal Statistical
  Society B* **39**, 172–212 (1977).
- S. Sengupta, T. Jovanovic-Talisman, D. Skoko, M. Renz, S. L. Veatch, J.
  Lippincott-Schwartz. "Probing protein heterogeneity in the plasma membrane
  using PALM and pair correlation analysis." *Nature Methods* **8**, 969–975
  (2011). [doi:10.1038/nmeth.1704](https://doi.org/10.1038/nmeth.1704)

The KD-tree range queries are provided by
[NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl).
