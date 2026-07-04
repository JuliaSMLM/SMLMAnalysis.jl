```@meta
CurrentModule = SMLMAnalysis
```

# Rendering

A localization list is not an image — the rendering step turns the `BasicSMLD`
into a super-resolution picture by accumulating every emitter onto a finely
subdivided pixel grid. It is selected by a `RenderConfig` and backed by
[SMLMRender](https://github.com/JuliaSMLM/SMLMRender.jl), which offers several
rendering strategies, intensity- or field-based coloring, and PNG export.

```julia
analyze(smld, RenderConfig(zoom = 20, colormap = :inferno))   # → (smld, StepInfo)
```

Rendering is a **pass-through**: it writes the image to disk and returns the
**input `smld` unchanged**. That is deliberate — a render can sit anywhere in the
pipeline and be repeated at different zooms or colormaps without altering the
data flowing to later steps. The pipeline treats it as a non-SMLD-producing step.

## When to use / prerequisites

- Run on a `BasicSMLD` of localizations — typically as the **final** step, after
  [detection/fitting](@ref "Detection & Fitting"), a [quality filter](@ref
  "Quality Filter"), [frame connection](@ref "Frame Connection"), and
  [drift correction](@ref "Drift Correction").
- Because it does not consume or transform the data, it is freely
  **repeatable**: drop in several `RenderConfig` steps (different strategies,
  zoom levels, colormaps, or `color_by` fields) and each writes its own image.
- To capture the image array in code (rather than only on disk), call
  `SMLMRender.render(smld, cfg)` or the step helper `render_step(smld, cfg)`
  directly — both return `(image, RenderInfo)`.

## Inputs, returns & artifacts

- **Input:** the current `smld`.
- **Returns:** `(smld, StepInfo)` — the *same* `smld` passes through. The
  `RenderInfo` (output size, pixel size, strategy, color mode) is on
  `StepInfo.info`.
- **Artifacts** (when `outdir` is set): a PNG written into the step folder with a
  descriptive name built from the run — `<strategy>_<colormap>_<color_by>_<zoom>x.png`
  (e.g. `gaussianrender_inferno_20x.png`); `config.toml` and `info.toml`; and at
  `STANDARD` verbosity a `stats.md`. No checkpoint is written (nothing changes).
  If you set `RenderConfig(filename = ...)` yourself, that path is used instead of
  the auto-generated name.

## Concept

Each localization is placed on an output grid finer than the camera pixels. The
**strategy** decides what gets drawn at each position; the **coloring** decides
what value maps to color. See the
[SMLMRender documentation](https://github.com/JuliaSMLM/SMLMRender.jl) for the
details of each.

- **Strategies** (`strategy`): `HistogramRender()` (fast binning),
  `GaussianRender()` (smooth blobs, the default), `CircleRender()`
  (anti-aliased circles at localization precision), and `EllipseRender()`.
- **Coloring**: *intensity-based* — accumulate counts and apply `colormap`
  (`:inferno`, `:hot`, `:magma`); or *field-based* — set `color_by` to an emitter
  field (`:z`, `:photons`, `:frame`, `:σ_x`) and color by that value with a
  perceptual map like `:turbo` or `:viridis`. `CircleRender`/`EllipseRender`
  require `color_by` or a manual `color`.
- **Resolution**: `zoom` renders the exact camera FOV with subdivided pixels
  (`zoom=20` on a 128×128 camera → 2560×2560), giving predictable, reproducible
  sizes; `pixel_size` (in nm) instead uses the data bounds with a margin. `roi`
  (zoom mode) renders a subset of the FOV as camera pixel ranges.

## Configuration

`RenderConfig` is re-exported from SMLMRender (the "upstream owns the config"
idiom). Its fields map 1:1 to `render()` keyword arguments. The ones you set most:

| field | typical/default | meaning |
|---|---|---|
| `strategy` | `GaussianRender()` | rendering algorithm: `HistogramRender`, `GaussianRender`, `CircleRender`, `EllipseRender` |
| `zoom` | `20` (`nothing`) | subdivide camera pixels by this factor; renders exact FOV |
| `pixel_size` | `nothing` | output pixel size in nm (data-bounds mode); alternative to `zoom` |
| `roi` | `nothing` | render a camera-pixel subregion, `(x_range, y_range)` (zoom mode) |
| `colormap` | `:inferno` (`nothing`) | intensity colormap |
| `color_by` | `nothing` | emitter field for field-based coloring (`:z`, `:photons`, `:frame`) |
| `clip_percentile` | `0.99` | intensity clip before normalizing; `nothing` = saturate mode |
| `categorical` | `false` | categorical palette for integer fields |
| `scalebar` | `false` | overlay a scale bar |
| `scalebar_length` | `nothing` | scale-bar length in µm (`nothing` = auto) |
| `scalebar_position` | `:br` | corner: `:br`, `:bl`, `:tr`, `:tl` |
| `backend` | `:cpu` | compute backend (`:cpu`, `:cuda`, `:metal`) |

See the SMLMRender documentation for the full field list (`target`, `color`,
`field_range`, `field_clip_percentiles`, `scalebar_color`, …) and defaults.

```julia
# Two renders of the same SMLD — both pass smld through unchanged
(smld, _) = analyze(smld, RenderConfig(zoom = 20, colormap = :inferno, scalebar = true))
(smld, _) = analyze(smld, RenderConfig(zoom = 10, color_by = :photons, colormap = :viridis))

# Get the image array directly
(img, rinfo) = SMLMRender.render(smld, RenderConfig(zoom = 20))
```

## Output & interpretation

The step's `StepInfo.summary` reports:

| field | meaning |
|---|---|
| `n_locs` | emitters actually rendered (`RenderInfo.n_emitters_rendered`) |
| `strategy` | strategy used (`:gaussian`, `:histogram`, `:circle`, `:ellipse`) |
| `output_size` | `(height, width)` of the output image in pixels |

`stats.md` additionally records the output pixel size (nm), color mode, and
render time. Sanity checks: `output_size` should equal `zoom ×` the camera
dimensions (within the `roi`) in zoom mode; `n_locs` should match the input
localization count unless some fell outside the rendered FOV. If a histogram
render looks washed out, lower `clip_percentile` (or set it to `nothing` for
saturate mode); for field coloring, confirm the chosen `color_by` field exists on
the emitters.

## Notes & caveats

- **Pass-through, not a transform.** The returned `smld` is the input; later
  steps see the data exactly as before. Position the render anywhere and repeat
  it freely.
- **No checkpoint.** Because nothing changes, no SMLD checkpoint is written for
  this step.
- **Coordinates are in microns**; `zoom`/`roi` are expressed in camera pixels,
  `pixel_size`/`scalebar_length` in nm/µm respectively.
- **GPU rendering** is available via `backend = :cuda` / `:metal` for large
  outputs.

## References

No external citation — super-resolution histogram/Gaussian rendering is a
standard SMLM visualization. See the
[SMLMRender documentation](https://github.com/JuliaSMLM/SMLMRender.jl) for the
rendering strategies, coloring modes, and all configuration options.
