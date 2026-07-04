```@meta
CurrentModule = SMLMAnalysis
```

# Density Filter

Isolated localizations — single fits with no nearby support — are usually noise:
spurious detections, mis-fits, or unlabeled background. The density filter drops
them by counting each localization's neighbors and rejecting those with too few.
It is selected by a `DensityFilterConfig` and is a **native** SMLMAnalysis step
(built on a KD-tree from
[NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl)); no
upstream package backs it.

```julia
analyze(smld, DensityFilterConfig())   # → (filtered_smld, StepInfo)
```

## When to use / prerequisites

- Run on a `BasicSMLD` of localizations, typically **late** in the pipeline —
  after [drift correction](@ref "Drift Correction") and
  [frame connection](@ref "Frame Connection"). Those steps sharpen and
  consolidate localizations, so the neighbor counts reflect true structure rather
  than blur or repeated blinks.
- Most useful when real structures are genuinely denser than the noise floor, so
  the neighbor-count distribution separates into an *isolated* population and a
  *clustered* population. Sparse or uniformly distributed data has no such valley
  to exploit — set an explicit cutoff instead (see below).

## Inputs, returns & artifacts

- **Input:** the current `smld`.
- **Returns:** `(filtered_smld, StepInfo)`. The per-step
  [`DensityFilterInfo`](@ref) is on `StepInfo.info`.
- **Artifacts** (when `outdir` is set): `stats.md` (input/output/rejected counts
  and the chosen threshold) and `neighbor_histogram.png` (the neighbor-count
  distribution with the threshold marked) at `STANDARD` verbosity, plus the saved
  config and info. A filtered SMLD is checkpointed to `smld_density.jld2` at
  `Checkpoint.ALL`.

## Concept

For every localization the step counts how many other localizations lie within
`n_sigma` *combined localization uncertainties*. The neighbor test is per-pair
and precision-aware: localizations ``i`` and ``j`` are neighbors when their
center-to-center distance satisfies

```math
d_{ij} < \texttt{n\_sigma} \cdot \sqrt{\sigma_i^2 + \sigma_j^2},
```

where each ``\sigma = \sqrt{\sigma_x^2 + \sigma_y^2}`` is the radial CRLB
uncertainty (µm). Because the radius scales with each pair's own precision, dense
high-precision regions (small ``\sigma``) require closer neighbors to count,
while loosely localized points get a proportionally larger search radius. A 2-D
KD-tree over the ``(x, y)`` coordinates (µm) does a coarse range query
(`n_sigma · 2 · max σ`), then the exact per-pair test refines each candidate.

The keep/reject threshold is a minimum neighbor count. With `min_neighbors = :auto`
the step picks it by **valley detection** on the neighbor-count histogram: it
smooths the histogram, finds the rightmost significant peak (the clustered
population), and places the threshold at the local minimum (valley) between the
origin and that peak. If the distribution is not clearly bimodal it falls back to
a conservative heuristic (and warns) rather than guess a valley.

## Configuration

| field | typical/default | meaning |
|---|---|---|
| `n_sigma` | `2.0` | neighbor radius in combined-uncertainty units; ``j`` counts as a neighbor of ``i`` when ``d_{ij} < \texttt{n\_sigma}\sqrt{\sigma_i^2+\sigma_j^2}`` |
| `min_neighbors` | `:auto` | minimum neighbor count to keep a localization; `:auto` chooses it by valley detection on the histogram, or pass an `Int` for an explicit cutoff |

See the [API Reference](@ref) for the complete field list.

```julia
# Automatic threshold (valley detection)
(filtered, info) = analyze(smld, DensityFilterConfig())

# Explicit cutoff: keep only localizations with ≥ 2 neighbors within 3σ
(filtered, info) = analyze(smld, DensityFilterConfig(n_sigma = 3.0, min_neighbors = 2))

info.info.threshold   # the neighbor-count threshold actually applied
```

## Output & interpretation

`StepInfo.summary` reports the headline numbers:

| field | meaning |
|---|---|
| `n_before` | localizations entering the step |
| `n_after` | localizations kept |
| `n_rejected` | `n_before - n_after` (the isolated localizations dropped) |
| `threshold` | the minimum-neighbor threshold applied (chosen or supplied) |

`DensityFilterInfo` carries the same `n_before`, `n_after`, `threshold` plus
`elapsed_s`.

Sanity checks: open `neighbor_histogram.png` and confirm the red threshold line
sits in the valley between the low-count (isolated) and high-count (clustered)
populations. A reasonable filter rejects a modest tail of isolated points, not
the bulk of the data. Heed the warnings: `:auto` emits one when the distribution
peaks at very low neighbor counts ("most emitters appear isolated") or is
unimodal — in those cases it returns a conservative default (`3`) or keeps almost
everything (`1`), which usually means density filtering is the wrong tool for
that dataset or that an explicit `min_neighbors` should be set.

## Notes & caveats

- **Run it late, and order matters.** Neighbor counts depend on the entire
  localization set, so running before drift correction or frame connection
  inflates apparent isolation. The step is repeatable — re-run with a different
  `n_sigma`/`min_neighbors` to tune.
- **`:auto` assumes bimodality.** If the isolated and clustered populations do not
  separate cleanly, prefer an explicit integer `min_neighbors`.
- **It is a pure spatial filter.** Coordinates and uncertainties are unchanged;
  only emitters are removed.
- The coarse KD-tree radius uses the *global* maximum ``\sigma``; the exact
  per-pair uncertainty test is what determines neighbor membership.

## References

This is a native SMLMAnalysis step with no external method citation; the
automatic threshold is the histogram valley-detection heuristic described above.
See the [API Reference](@ref) for `DensityFilterConfig` and
[`DensityFilterInfo`](@ref).
