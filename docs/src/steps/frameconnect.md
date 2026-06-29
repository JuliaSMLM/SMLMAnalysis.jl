```@meta
CurrentModule = SMLMAnalysis
```

# Frame Connection

A single fluorophore typically blinks on for several consecutive frames before
going dark, so one emitter generates a short burst of repeated localizations.
Frame connection recognizes those repeats and links them into a single,
higher-precision localization — the average of a track is more precise than any
one blink. It is selected by a `FrameConnectConfig` and backed by
[SMLMFrameConnection](https://github.com/JuliaSMLM/SMLMFrameConnection.jl), which
solves a spatiotemporal **Linear Assignment Problem** (LAP) over proximity and
estimated blinking kinetics.

```julia
analyze(smld, FrameConnectConfig(max_frame_gap = 5))   # → (combined_smld, StepInfo)
```

## When to use / prerequisites

- Run on a `BasicSMLD` of `Emitter2DFit` localizations, after
  [detection/fitting](@ref "Detection & Fitting") and usually after a
  [quality filter](@ref "Quality Filter"). It is normally run before
  [drift correction](@ref "Drift Correction").
- Inputs must carry positive position uncertainties (`σ_x`, `σ_y` > 0) — these
  weight the LAP cost and the MLE combination. Frames are 1-based.
- Most useful when the same dye fluoresces over multiple frames (dSTORM,
  DNA-PAINT). It also lets you drop single-frame events as noise via
  `track_length`.

## Inputs, returns & artifacts

- **Input:** the current `smld`.
- **Returns:** `(combined_smld, StepInfo)`. `combined_smld` holds the recombined,
  higher-precision localizations (one per track). The full upstream
  `FrameConnectInfo` is on `StepInfo.info`; its `.connected` field is the
  *uncombined* input with a `track_id` assigned to each localization, which the
  pipeline also surfaces on [`AnalysisResult`](@ref)`.smld_connected` for
  track-level analysis.
- **Artifacts** (when `outdir` is set): the config and `FrameConnectInfo` are
  written; at `STANDARD` verbosity `track_histogram.png` (localizations-per-track
  distribution with mean/median and an estimated `k_off`) and `stats.md` are
  added. The combined SMLD is checkpointed to `smld_combined.jld2` at
  `Checkpoint.EXPENSIVE`. When calibration is enabled, `uncertainty_calibration.png`,
  `shift_histogram.png`, and `drift_jitter.png` are also written.

## Concept

After preclustering spatially/temporally adjacent detections, the algorithm
estimates the fluorophore's on/off/bleach rates and emitter density from the data
itself, then assigns localizations to emitters by solving a LAP whose costs trade
off spatial distance against blinking probability. Connected localizations are
then combined by an MLE weighted mean. See the SMLMFrameConnection documentation
for the full algorithm; the method is from Schodt & Lidke (2021).

### Uncertainty calibration (optional, inside this step)

Calibration is **not** a separate pipeline step — it is configured *inside*
`FrameConnectConfig` via its `calibration=` field (a `CalibrationConfig`). When
set, frame-to-frame jitter within tracks is analyzed to fit
`observed_var = A + B·CRLB_var`, recovering an additive motion variance
(`σ_motion² = A`) and a CRLB scale factor (`k² = B`). Corrected uncertainties
`Σ = σ_motion²·I + k²·Σ_CRLB` are applied *before* combination, so tracks are
weighted with calibrated errors in a single pass. The diagnostics land on
`FrameConnectInfo.calibration` (a `CalibrationResult`).

## Configuration

`FrameConnectConfig` is re-exported from SMLMFrameConnection (the "upstream owns
the config" idiom). The fields you set most often:

| field | typical/default | meaning |
|---|---|---|
| `max_frame_gap` | `5` (dSTORM `10`–`20`) | max dark-frame gap bridged within a track; raise for long dark states |
| `max_sigma_dist` | `5.0` | distance threshold as a multiple of localization error; higher links over larger distances |
| `track_length` | `nothing` | inclusive `(min, max)` locs-per-track filter; `(2.0, Inf)` drops single-frame blinks |
| `max_neighbors` | `2` | nearest neighbors inspected during preclustering |
| `n_density_neighbors` | `2` | preclusters used for the local density estimate |
| `calibration` | `nothing` | a `CalibrationConfig` to enable uncertainty calibration (see above) |

`CalibrationConfig` fields: `clamp_k_to_one` (default `true`, since CRLB is a
lower bound), `filter_high_chi2` (default `false`), and `chi2_filter_threshold`
(default `6.0`). See the SMLMFrameConnection documentation for the complete list.

```julia
# Basic frame connection, dropping single-frame blinks
(combined, info) = analyze(smld,
    FrameConnectConfig(max_frame_gap = 5, track_length = (2.0, Inf)))

# With uncertainty calibration
(combined, info) = analyze(smld,
    FrameConnectConfig(max_frame_gap = 5,
                       calibration = CalibrationConfig(filter_high_chi2 = true)))

connected = info.info.connected   # uncombined locs with track_id
```

## Output & interpretation

The step's `StepInfo.summary` reports the headline numbers:

| field | meaning |
|---|---|
| `n_before` | input localizations (`FrameConnectInfo.n_input`) |
| `n_after` | output tracks / combined localizations (`n_combined`) |
| `n_filtered` | tracks dropped by the `track_length` filter |
| `compression` | `n_before / n_after` — mean localizations per track |
| `k_scale`, `sigma_motion_nm`, `mean_chi2` | added only when calibration was applied |

Sanity checks: `compression` should be ≳ 1 (a few-to-many for blinky dyes); a
value of ~1 means almost nothing connected — `max_frame_gap`/`max_sigma_dist`
may be too tight, or the data is genuinely sparse. When calibration runs,
`mean_chi2` should be near `2.0` (two axes) if uncertainties are well-calibrated,
and `k_scale ≥ 1` reflects CRLB being a lower bound.

## Notes & caveats

- **Calibration lives here, not as its own step.** Enable it through
  `calibration=CalibrationConfig(...)`; its results appear on
  `FrameConnectInfo.calibration`, never as a standalone `StepInfo`.
- **Combined vs. connected.** The returned `smld` is the *combined* set (one loc
  per track). If you need per-blink localizations with their track assignment,
  use `info.info.connected` / [`AnalysisResult`](@ref)`.smld_connected`.
- **Tune `max_frame_gap` to the dye.** dSTORM dyes with long dark states want a
  larger gap (10–20); too large over-links distinct emitters, too small fragments
  one emitter into many tracks.
- **2D only.** SMLMFrameConnection operates on `Emitter2DFit` data.

## References

- D. J. Schodt and K. A. Lidke. "Spatiotemporal Clustering of Repeated
  Super-Resolution Localizations via Linear Assignment Problem." *Frontiers in
  Bioinformatics* **1**, 724325 (2021).
  [doi:10.3389/fbinf.2021.724325](https://doi.org/10.3389/fbinf.2021.724325)

See the [SMLMFrameConnection documentation](https://github.com/JuliaSMLM/SMLMFrameConnection.jl)
for the algorithm in full and all configuration options.
