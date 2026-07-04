```@meta
CurrentModule = SMLMAnalysis
```

# Intensity Filter

When two fluorophores are active at once inside a single diffraction-limited
spot, the fitter sees one blob and returns a single localization with a corrupted
position and an abnormally **high photon count**. The intensity filter rejects
these multi-emitter events with a statistical test on brightness. It is a step
**native to SMLMAnalysis**, selected by an `IntensityFilterConfig`.

```julia
analyze(smld, IntensityFilterConfig(cutoff = 0.01))   # → (filtered_smld, StepInfo)
```

## When to use / prerequisites

- Run on a fitted `BasicSMLD` whose emitters carry `photons` (i.e. after
  [Detection & Fitting](@ref "Detection & Fitting")).
- Most valuable on **denser data**, where simultaneous activations within one PSF
  are common (e.g. high-density dSTORM/PAINT). On very sparse data it has little
  to reject.
- It needs enough localizations to estimate the excitation field — see
  `min_bin_count` below.

## Inputs, returns & artifacts

- **Input** — the current `smld`.
- **Returns** — `(filtered_smld, StepInfo)`; the `StepInfo.info` is an
  `IntensityFilterInfo`.
- **Artifacts** (when `outdir` is set) — diagnostic figures of the fitted
  excitation field and the per-localization p-value distribution, plus the usual
  `config.toml` / `info.toml`.

## Concept

A single emitter does not emit a fixed number of photons — the count is Poisson,
and its mean varies across the field of view with the excitation intensity. The
filter accounts for both:

1. **Estimate the expected single-emitter rate, spatially.** The field of view is
   binned (`n_bins` per axis); in each bin the `rate_percentile` (default 95th)
   of the photon distribution is taken as the upper bound of *single-emitter*
   emission there — high enough to sit above typical single-emitter counts but
   below true doubles. A smooth field model is then fit across bins: a
   spatially-varying Gaussian beam (`field_mode = :gaussian`) or a single global
   rate (`:uniform`).
2. **Test each localization.** Given its position's expected rate λ, a
   **Poisson upper-tail test** asks how improbable the observed photon count is
   under single-emitter emission. Localizations with p-value below `cutoff` are
   rejected as multi-emitter events.

Optionally (`estimate_p2`), the filter also estimates **p₂**, the fraction of
double-emitter events, by decomposing the bright tail of the photon distribution
— a useful diagnostic of how crowded the data is.

## Configuration

| field | default | meaning |
|---|---|---|
| `cutoff` | `0.01` | p-value cutoff; localizations with p < cutoff are rejected |
| `field_mode` | `:gaussian` | excitation-field model: `:gaussian` (spatially-varying beam) or `:uniform` (one global rate) |
| `n_bins` | `10` | spatial grid bins per axis for field estimation |
| `min_bin_count` | `30` | minimum localizations per bin to estimate its rate |
| `rate_percentile` | `0.95` | per-bin percentile taken as the single-emitter rate; higher = more permissive |
| `estimate_p2` | `true` | also estimate the double-emitter fraction p₂ |
| `p2_tail_threshold`, `p2_n_bins` | `1.0`, `200` | p₂ tail-decomposition parameters |

```julia
(filtered, info) = analyze(smld,
    IntensityFilterConfig(cutoff = 0.01, field_mode = :gaussian, rate_percentile = 0.95))

info.info.p2_estimate    # estimated double-emitter fraction (or nothing)
```

## Output & interpretation

`IntensityFilterInfo` reports:

| field | meaning |
|---|---|
| `n_before`, `n_after` | localizations in and surviving |
| `field_mode` | the field model actually used |
| `lambda_max_global` | peak expected single-emitter rate (photons) |
| `field_fit_r2` | goodness of the field fit |
| `p2_estimate` | estimated fraction of double-emitter events (or `nothing`) |
| `p2_tail_obs`, `p2_tail_f2` | observed vs. modeled bright-tail mass |

Sanity checks: `field_fit_r2` should be high for a clean `:gaussian` fit — if it
is poor, the beam assumption may not hold and `:uniform` is the safer choice. A
large `p2_estimate` signals crowded data where multi-emitter rejection matters
most; if it is near zero, the filter is mostly a no-op and can be dropped.

## Notes & caveats

- **Native, in-house method** — the per-bin-percentile + smooth-field +
  Poisson-tail approach is implemented in SMLMAnalysis; there is no separate
  upstream package for it.
- **Field estimation needs data.** With too few localizations per bin (below
  `min_bin_count`) the field estimate is unreliable; reduce `n_bins` or fall back
  to `field_mode = :uniform`.
- **Repeatable.** Like the other filters, it can appear more than once in a
  pipeline, though once is typical.

## References

This is a method native to SMLMAnalysis; there is no separate primary
publication. See the [API Reference](@ref) for the full `IntensityFilterConfig`
and `IntensityFilterInfo` field lists.
