```@meta
CurrentModule = SMLMAnalysis
```

# Quality Filter

Single-molecule localization produces a population of fits of wildly varying
quality — bright clean blinks alongside dim, mis-fit, or multi-emitter events.
The quality-filter step keeps only the localizations that pass per-property
threshold ranges, discarding the rest. It is selected by a `FilterConfig` and is
a **native SMLMAnalysis step** (no upstream package): the thresholding logic
lives in this package.

```julia
analyze(smld, FilterConfig(photons = (500.0, Inf)))   # → (filtered_smld, StepInfo)
```

## When to use / prerequisites

- Run on a `BasicSMLD` of localizations — typically right after
  [detection/fitting](@ref "Detection & Fitting") to cut weak and spurious fits
  before they reach later steps.
- The step is **repeatable** (see [The Pipeline Model](@ref)): a common pattern
  is a *coarse* filter early (e.g. a loose photon floor) and a *tighter* filter
  after [frame connection](@ref "Frame Connection"), once linked blinks have
  sharper precision and goodness-of-fit estimates.
- Filtering on `precision`, `psf_sigma`, `z`, or `sigma_z` requires emitters that
  actually carry those fields; the `z` / `sigma_z` filters are **no-ops on 2D
  SMLDs** (emitters without `z` / `σ_z`).

## Inputs, returns & artifacts

- **Input:** the current `smld`.
- **Returns:** `(filtered_smld, StepInfo)`. The filter's own
  [`FilterInfo`](@ref) (`n_before`, `n_after`, `elapsed_s`) is on `StepInfo.info`.
- **Artifacts** (when `outdir` is set, at `STANDARD` verbosity): `stats.md`
  (counts, acceptance %, criteria applied), `fit_quality.png` (photon,
  background, precision, p-value, PSF-width — and for 3D fits z and σ_z —
  distributions with the threshold/rejected regions drawn), a `fit_overlay.png`
  of sampled ROIs colored by pass/fail (drawn only when a detect-fit sample cache
  is present; gracefully skipped for a standalone filter), and a post-filter
  `localizations_per_frame.png`. At `DETAILED` verbosity a `detailed_stats.md`
  gives per-criterion pass/fail counts. A filtered SMLD is checkpointed to
  `smld_filtered.jld2` at `Checkpoint.ALL`.

## Concept

Each criterion is a closed interval. An emitter is kept only if **every**
enabled criterion passes — the filters combine as a logical AND over the
population. All criteria are `(min, max)` tuples (use `-Inf` / `Inf` for an
open bound), and **every filter defaults to `nothing`, i.e. disabled**, so you
pay only for the criteria you set.

- `photons` — total fitted photons; a floor removes dim, poorly-localized blinks.
- `precision` — lateral localization precision in microns, tested as
  `max(σ_x, σ_y)`. These σ are the *exact* CRLB fit uncertainties from the MLE
  Fisher information (Smith et al. 2010; Huang et al. 2013 for sCMOS), not an
  analytical approximation; a precision ceiling is the most direct way to bound
  reconstruction resolution.
- `pvalue` — goodness-of-fit p-value from the fit's log-likelihood-ratio statistic
  (Huang et al. 2011); a floor rejects fits whose residuals are inconsistent with
  the PSF model (often multi-emitter or contaminated ROIs).
- `psf_sigma` — fitted PSF width. Real single molecules cluster tightly around
  one width; out-of-band widths flag defocus or blends. Pass an explicit
  `(min, max)` in microns, or `:auto` to bound the population **mode ± 10 %**
  (computed per axis for anisotropic `σx`/`σy` models).
- `z`, `sigma_z` — axial position and axial precision (microns) for 3D fits;
  the axial analogs of `photons`-window and `precision`. The `z` window is
  handy to drop localizations railed to the z-grid edge.

## Configuration

All fields are `(min, max)` tuples (microns where a length), `nothing` to
disable; `psf_sigma` also accepts the `:auto` symbol.

| field | typical | meaning |
|---|---|---|
| `photons` | `(500.0, Inf)` | photon-count window; floor removes dim fits |
| `precision` | `(0.0, 0.007)` | lateral precision `max(σ_x, σ_y)` window, µm (≤7 nm here) |
| `pvalue` | `(1e-3, 1.0)` | goodness-of-fit p-value window; floor rejects bad fits |
| `psf_sigma` | `:auto` | PSF-width window, µm — or `:auto` for mode ± 10 % |
| `z` | `(-0.4, 0.4)` | axial-position window, µm (3D only; no-op on 2D) |
| `sigma_z` | `(0.0, 0.030)` | axial-precision window, µm (3D only; no-op on 2D) |

See the [API Reference](@ref) for the complete field list.

```julia
# Coarse cut after detection…
(smld, _) = analyze(smld, FilterConfig(photons = (500.0, Inf)))

# …then a tighter pass after frame connection
(smld, info) = analyze(smld, FilterConfig(
    photons   = (1000.0, Inf),
    precision = (0.0, 0.007),     # ≤ 7 nm lateral
    pvalue    = (1e-3, 1.0),
    psf_sigma = :auto,            # keep widths within mode ± 10 %
))

info.info.n_after / info.info.n_before    # acceptance fraction
```

## Output & interpretation

`StepInfo.summary` reports the headline numbers:

| field | meaning |
|---|---|
| `n_before` | localizations entering the filter |
| `n_after` | localizations kept |
| `acceptance` | `n_after / n_before`, rounded to 3 digits |

Sanity checks: inspect `fit_quality.png` and confirm the threshold lines sit on
the *shoulders* of each distribution, not through its peak — a filter that
removes most of the data is usually mis-set (e.g. a precision ceiling tighter
than the fits can achieve, or a photon floor above the median). For 3D data a
spike at the edge of the `z` panel signals grid-railed degenerate fits that the
`z` window should remove. The per-criterion `detailed_stats.md` (DETAILED
verbosity) shows which single criterion is responsible for most of the loss.

## Notes & caveats

- **Order matters for downstream precision.** Filtering on `precision` /
  `psf_sigma` before frame connection uses single-blink CRLB values; the same
  filter after connection acts on the combined, higher-precision tracks. Decide
  which population you mean to threshold.
- **`:auto` follows the data.** The mode ± 10 % band adapts to each dataset's
  fitted-width distribution, so the *same* `FilterConfig(psf_sigma = :auto)`
  yields different bounds on different acquisitions — a feature for batch runs,
  but check the bounds drawn in `fit_quality.png` if a dataset is unusual.
- **3D filters silently no-op on 2D.** `z` / `sigma_z` are skipped when emitters
  lack those fields, so leaving them set in a shared config is harmless.
- **Disabled by default.** An empty `FilterConfig()` passes everything through
  unchanged.

## References

- C. S. Smith, N. Joseph, B. Rieger, K. A. Lidke. "Fast, single-molecule
  localization that achieves theoretically minimum uncertainty." *Nature Methods*
  **7**, 373–375 (2010).
  [doi:10.1038/nmeth.1449](https://doi.org/10.1038/nmeth.1449) — the exact CRLB
  (from the MLE Fisher information) behind the `σ` values the `precision` criterion
  bounds.
- F. Huang, T. M. P. Hartwich, F. E. Rivera-Molina, *et al.* "Video-rate nanoscopy
  using sCMOS camera-specific single-molecule localization algorithms." *Nature
  Methods* **10**, 653–658 (2013).
  [doi:10.1038/nmeth.2488](https://doi.org/10.1038/nmeth.2488) — extends the exact
  CRLB to per-pixel sCMOS noise.
- F. Huang, S. L. Schwartz, J. M. Byars, K. A. Lidke. "Simultaneous
  multiple-emitter fitting for single molecule super-resolution imaging."
  *Biomedical Optics Express* **2**, 1377–1393 (2011).
  [doi:10.1364/BOE.2.001377](https://doi.org/10.1364/BOE.2.001377) — the
  log-likelihood-ratio goodness-of-fit test behind the `pvalue` criterion.

This is a native SMLMAnalysis step; see the [API Reference](@ref) for the full
`FilterConfig` field list.
