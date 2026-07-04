```@meta
CurrentModule = SMLMAnalysis
```

# Bayesian Grouping (BaGoL)

A single fluorophore blinks many times during an acquisition, so it appears not
as one point but as a scatter of localizations around its true position. Bayesian
Grouping of Localizations (BaGoL) groups those repeated blinks back into the
individual emitters that produced them, yielding emitter positions **more precise
than any single localization** — and, as a byproduct, a count of how many emitters
are actually present. The step is selected by a `BaGoLConfig` and backed by
[SMLMBaGoL](https://github.com/JuliaSMLM/SMLMBaGoL.jl), which runs a reversible-jump
MCMC (RJMCMC) sampler over the number and positions of emitters.

```julia
analyze(smld, BaGoLConfig(μ = 10.0, shape = 2.0))   # → (grouped_smld, StepInfo)
```

## When to use / prerequisites

Use BaGoL for counting and for maximum positional precision — typically as an
alternative or supplement to a final [Rendering](@ref "Rendering"): group first, then render
the grouped emitters. It is normally the last analysis step.

BaGoL makes two assumptions, both of which dictate preprocessing:

- **Each localization is a real observation of a single emitter.** Run BaGoL on
  cleaned localizations — apply a [Quality Filter](@ref) and an
  [Intensity Filter](@ref) first to drop spurious detections and multi-emitter
  fits (where one localization stands for two or more emitters), which BaGoL
  cannot model.
- **The reported precision `σ` is correct.** BaGoL groups by each localization's
  fitted `σ_x`, `σ_y`, so it expects [Frame Connection](@ref) (and ideally the
  uncertainty calibration configured inside `FrameConnectConfig`) to have run, and
  it defaults to `se_adjust = :auto` to repair any residual `σ` underestimate (see
  below). Underestimated `σ` causes BaGoL to over-split one emitter into several.

## Inputs, returns & artifacts

- **Input:** the current `smld` of 2D fitted localizations (`Emitter2DFit`), each
  carrying its uncertainty `σ_x`, `σ_y`. Positions and uncertainties are in microns.
- **Returns:** `(grouped_smld, StepInfo)`. `grouped_smld.emitters` are the **grouped
  emitters** (`Emitter2DFit`, each with a posterior position uncertainty), not the
  input localizations. The full `SMLMBaGoL.BaGoLDiagnostics` is on
  `StepInfo.info.diagnostics`.
- **Artifacts** (when `outdir` is set, at `STANDARD` verbosity, in `NN_bagol/`):
  `config.toml` and `info.toml`; a metrics report and posterior diagnostic figures
  (upstream `write_report` / `plot_report`); a render set written by upstream
  `render_report` — `render_mapn`, `render_sr`, `render_circles`, and
  `render_partitions` (at `zoom = 50`); and, when `se_adjust = :auto`, a finder
  diagnostic plot (`plot_se_adjust`). The grouped SMLD is checkpointed to
  `smld_bagol.jld2` at `Checkpoint.EXPENSIVE`.

## Concept

BaGoL treats the localizations as measurements from a mixture of emitters whose
number *and* positions are both unknown. Its RJMCMC moves split, merge, create, and
remove emitters, so the chain explores models with different emitter counts rather
than fixing one. Instead of committing to a single grouping it builds a full
posterior over both the emitter count `K` and the positions; for downstream use the
posterior is summarized by the **MAP-N** estimate — the most probable count with a
representative grouping at that count — pooling each emitter's localizations to a
super-resolved position. For the priors, the collapsed sampler, and the MAP-N
estimators see the SMLMBaGoL documentation.

## Configuration

`BaGoLConfig` is re-exported from SMLMBaGoL (the "upstream owns the config" idiom);
its fields map 1:1 onto `run_bagol` kwargs. The fields you set most often:

| field | typical / default | meaning |
|---|---|---|
| `μ` | `10.0` | NegBin mean localizations (blinks) per emitter (typed `\mu`+Tab) |
| `shape` | `2.0` | count-distribution shape (`1` ≈ dSTORM/exponential, `>1` ≈ peaked/DNA-PAINT) |
| `learn_distribution` | `true` | learn count params: `true` (both), `false` (fix both), `:mu`, or `:shape` |
| `n_iterations` | `4000` | total MCMC iterations |
| `burn_in` | `2000` | iterations discarded before accumulating |
| `se_adjust` | `:auto` | `:auto` estimates excess uncertainty `τ` and adds it in quadrature (`σ²+τ²`); `0.0` = none; a number = manual `τ` (µm) |
| `partition_sigma` | `3.0` | precision-weighted DBSCAN partition threshold in σ units (`Inf` = no partitioning) |

See the SMLMBaGoL documentation for the complete field list (model choices,
split/merge tuning, partitioning) and defaults.

```julia
# Group cleaned, frame-connected localizations into emitters
(grouped, info) = analyze(smld, BaGoLConfig(μ = 10.0, shape = 2.0,
                                            n_iterations = 4000, se_adjust = :auto))

grouped.emitters            # the grouped emitters (super-resolved positions)
info.info.n_emitters        # MAP-N emitter count
info.info.diagnostics       # full BaGoLDiagnostics (posterior_k, acceptance_rates, …)
```

## Output & interpretation

The step's `StepInfo.summary` reports the headline numbers:

| field | meaning |
|---|---|
| `n_locs_in` | input localization count |
| `n_emitters` | MAP-N number of grouped emitters |
| `compression` | `n_locs_in / n_emitters` — localizations pooled per emitter |
| `final_μ`, `final_shape` | learned count-distribution parameters |
| `n_partitions` | number of spatial partitions the data was split into |
| `tau_um` | applied excess uncertainty `τ̂` (µm); the `se_adjust` correction actually used |
| `se_adjust` | applied `(τx, τy)` provenance |

Sanity checks (drill into `info.info.diagnostics` for the full picture): a
plausible `compression` is a handful to tens of localizations per emitter; a
surprisingly large `final_μ` signals under-splitting (too few emitters, each
absorbing too many localizations). The diagnostics' `posterior_k` should be
*concentrated* on a count rather than smeared across many `K`, and the split /
merge / birth / death `acceptance_rates` should be nonzero — all-zero means the
chain is stuck at one `K`. A large `tau_um` means the reported `σ` was badly
underestimated upstream.

## Notes & caveats

- **Clean first.** Spurious and multi-emitter fits violate BaGoL's single-emitter
  assumption; filter them out before grouping or they corrupt the count.
- **Calibrated σ matters.** Underestimated precision over-splits emitters. Prefer
  calibrated, frame-connected input and keep `se_adjust = :auto` so the finder
  repairs residual error (it self-guards against double-correcting an already
  σ-corrected SMLD).
- **It is expensive.** BaGoL runs a full MCMC chain; this is the costliest step
  and is checkpointed by default. Large datasets are partitioned automatically and
  run in parallel — start Julia with `-t auto`.
- **The output is emitters, not localizations.** Downstream steps (e.g. a render)
  operate on the grouped emitters, so a BaGoL render looks sparser and sharper than
  the raw-localization render.

## References

- M. Fazel, *et al.* "High-Precision Estimation of Emitter Positions using Bayesian
  Grouping of Localizations." *Nature Communications* **13**, 7152 (2022).
  [doi:10.1038/s41467-022-34894-2](https://doi.org/10.1038/s41467-022-34894-2)

See the [SMLMBaGoL documentation](https://github.com/JuliaSMLM/SMLMBaGoL.jl) for the
RJMCMC sampler, priors, MAP-N estimators, and uncertainty correction in full.
