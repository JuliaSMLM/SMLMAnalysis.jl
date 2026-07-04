```@meta
CurrentModule = SMLMAnalysis
```

# Composite Render

A composite render overlays several single-color reconstructions into one
multi-color image: each channel is rendered to its own intensity map, tinted a
chosen color, and the tinted maps are summed into a single RGB picture. The step
is selected by a `CompositeRenderConfig` and is a **native SMLMAnalysis** step
built on [SMLMRender](https://github.com/JuliaSMLM/SMLMRender.jl) — it simply
calls `SMLMRender.render` with a per-channel color assignment.

Unlike the single-color [Rendering](@ref "Rendering") step, this is a **multi-channel** step:
it operates on a `Vector{BasicSMLD}` (one SMLD per color/target), not a single
SMLD, and is dispatched accordingly.

```julia
analyze(smlds, CompositeRenderConfig(zoom = 20))   # smlds::Vector{BasicSMLD}
#   → (smlds, StepInfo)        # SMLDs pass through unmodified; image goes to disk
```

## When to use / prerequisites

- This is a [Multi-Channel](@ref "Multi-Channel Analysis") step. It needs
  **two or more channels** — one `BasicSMLD` per color — so it only makes sense
  in a multi-target workflow, not a single-color pipeline.
- It is almost always driven by [`MultiTargetConfig`](@ref): each channel runs
  its own single-color pipeline first, then composite render (and any
  [Cross-Alignment](@ref) / [Cross-Correlation](@ref)) runs over the resulting
  vector of SMLDs. Run it after the per-channel localizations are finished, and
  after any cross-alignment when you want the registered overlay.
- It produces an image only; it does **not** modify the localizations.

## Inputs, returns & artifacts

- **Input:** a `Vector{BasicSMLD}`, one entry per channel, in the same order as
  the channel colors.
- **Returns:** `(smlds, StepInfo)`. The SMLDs are returned **unchanged**
  (pass-through). The render metadata is on `StepInfo.info`, a
  [`CompositeRenderInfo`](@ref).
- **Artifacts** (when `outdir` is set), written to
  `outdir/NN_compositerender/`:
  - `<strategy>_<zoom>x.png` — the RGB composite image (e.g.
    `gaussianrender_20x.png`).
  - `config.toml`, `info.toml` — the step config and the upstream `RenderInfo`.
  - `stats.md` (at `STANDARD` verbosity) — strategy, zoom, per-channel
    localization counts, output size, pixel size, and timing.

  When driven by `MultiTargetConfig` the step directory lives under
  `outdir/composite/`, alongside a top-level `README.md` documenting the color
  scheme and channel labels.

## Concept

Each channel is rendered independently into a grayscale intensity image (using
the chosen [Rendering](@ref "Rendering") `strategy`), then mapped onto an RGB color and the
channels are added together. The result is the familiar two- or three-color
overlay where regions of co-localization appear as the blend of the constituent
colors (e.g. cyan + magenta → white). The actual rendering, normalization, and
color mixing are done by `SMLMRender.render(smlds; colors=...)`; see the
SMLMRender documentation for the rendering strategies and how intensities are
mapped to color.

Two display choices interact with the render strategy and are resolved
automatically when left at their defaults:

- **`normalize_each`** — whether each channel is normalized to its own dynamic
  range before mixing. Auto: `false` for `HistogramRender` (saturate mode),
  `true` for the others (clip + normalize).
- **`clip_percentile`** — intensity clipping for the normalize/clip path. For
  `HistogramRender` the default `0.99` is treated as "no clip" (saturate).

## Configuration

`CompositeRenderConfig` mirrors the rendering knobs of `SMLMRender.render`. The
fields that matter most:

| field | typical / default | meaning |
|---|---|---|
| `strategy` | `GaussianRender()` | rendering strategy: `GaussianRender`, `HistogramRender`, `CircleRender`, … |
| `zoom` | `20.0` | output magnification (super-resolution pixels per camera pixel) |
| `colors` | `nothing` | per-channel colors as `Vector{Symbol}`; `nothing` = inherit from `MultiTargetConfig.colors` |
| `clip_percentile` | `0.99` | intensity clip for normalize mode; `nothing` = saturate |
| `normalize_each` | `nothing` | per-channel normalization; `nothing` = auto (false for histogram, true otherwise) |
| `scalebar` | `true` | draw a scale bar on the composite |
| `scalebar_length` | `nothing` | scale-bar length in µm; `nothing` = auto |
| `scalebar_position` | `:br` | corner: `:br`, `:bl`, `:tr`, `:tl` |
| `scalebar_color` | `:white` | scale-bar color |

```julia
# Two-color overlay, colors inherited from the MultiTargetConfig
mt = MultiTargetConfig(
    labels = [:IgG, :C1q],          # colors default to cyan / magenta
    steps  = [CompositeRenderConfig(zoom = 20.0, strategy = GaussianRender())],
    outdir = "output/cell1/",
)
(result, info) = analyze(channels, mt)

# Or call the step directly on a vector of SMLDs
(smlds, sinfo) = analyze([smld_red, smld_green],
    CompositeRenderConfig(zoom = 20.0, colors = [:red, :green]);
    outdir = "output/")
```

If you set `colors` on the config it overrides the `MultiTargetConfig` defaults;
otherwise the channel colors flow down from the parent config.

## Output & interpretation

`StepInfo.summary` reports:

| field | meaning |
|---|---|
| `strategy` | rendering strategy used (Symbol) |
| `zoom` | magnification factor |
| `n_channels` | number of channels overlaid |
| `output_size` | composite image size in pixels |

Sanity checks: `n_channels` should equal the number of SMLDs you passed (and the
number of colors); `output_size` should scale with `zoom`. Inspect the PNG — the
channels should land on top of each other where the targets co-localize. If the
overlay looks misregistered, run a [Cross-Alignment](@ref) step before the
composite render.

## Notes & caveats

- **Pass-through.** The SMLDs are returned unchanged; the only product is the
  image. Re-render at a different `zoom`/`strategy` freely without re-running the
  pipeline.
- **Colors must match channels.** When driven by `MultiTargetConfig`, the number
  of colors must equal the number of channels (the driver enforces this).
- **Histogram vs. Gaussian.** `HistogramRender` defaults to saturate mode (no
  per-channel normalization, no clip); the smooth strategies normalize and clip.
  Override `normalize_each` / `clip_percentile` only if the auto behavior is not
  what you want.

## References

No external citation — this is standard per-channel RGB compositing. For the
rendering strategies, intensity normalization, and color mapping, see the
[SMLMRender documentation](https://github.com/JuliaSMLM/SMLMRender.jl), and the
[Multi-Channel](@ref "Multi-Channel Analysis") workflow page for the surrounding
multi-target pipeline.
