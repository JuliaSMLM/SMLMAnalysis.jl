```@meta
CurrentModule = SMLMAnalysis
```

# [Edge Classification](@id edgeclassify-step)

Cellular SMLM data has a boundary: localizations belong to a cell interior, sit on
its membrane, or fall outside it in the surrounding background. The edge
classification step carves one or more **cell masks** directly from the
localization cloud and labels every localization `:interior` / `:membrane` /
`:outside`. It is backed by
[SMLMClustering](https://github.com/JuliaSMLM/SMLMClustering.jl) and, like
clustering, is selected entirely by the **type** of its config:

- `KdeValleyConfig` — the dSTORM-tuned default. Finds the cell boundary from a KDE
  valley between the dense interior and the sparse background.
- `OuterPolygonConfig` — carves the outer ring of the dominant cell as a polygon.

Both are `<: AbstractEdgeClassifyConfig`, re-exported from SMLMClustering as `const`
aliases, so they are available unqualified after `using SMLMAnalysis`.

```julia
(smld, info) = analyze(smld, KdeValleyConfig())     # → (smld, StepInfo)
(smld, info) = analyze(smld, OuterPolygonConfig())
```

## Non-destructive by design

The step is **state-modifying but emitter-preserving**: the returned SMLD wraps the
*same* emitters — nothing is filtered out. The authoritative per-emitter class is
**not** written back onto the emitters (a per-emitter side-list would desync the
moment a downstream step subsets emitters); it lives in the step's
`EdgeClassifyInfo`, reached at the classify point:

```julia
(smld, step_info) = analyze(smld, KdeValleyConfig())
info = step_info.info                 # ::EdgeClassifyInfo
info.class                            # per-emitter Vector{Symbol}
interior_mask(info)                   # BitVector of the :interior emitters
in_cell(info)                         # :interior ∪ :membrane
interior_fraction(info)               # scalar summary
```

Only the cell-mask **geometry** travels downstream, mirrored into `metadata`:

- `metadata["edge_cells"]` — the `MultiCellMask` (`Vector{CellPolygon}`, largest-first)
- `metadata["edge_outer_polygon"]` — the dominant cell's outer ring (back-compat)

A later spatial-statistics step reads that geometry as its observation region — e.g.
`HopkinsConfig` uses `edge_cells` / `edge_outer_polygon` to restrict its clustering
tendency estimate to inside the cell.

## When to use

- Run late in the pipeline, on a cleaned `BasicSMLD`, when downstream analysis must
  distinguish cell interior from membrane or background — e.g. restricting a Hopkins
  clustering-tendency estimate to the cell, or reporting an interior fraction.
- Because it removes nothing, it is safe to run before clustering or rendering; the
  emitters and their ordering are unchanged.

## Outputs

At `verbose >= STANDARD` with an `outdir`, the step writes `config.toml`, `info.toml`,
and an edge report (diagnostics plus figures) under `{step}_edge_classify/`. The
StepInfo summary records the gate `:method`, emitter counts, and the number of cells.

See the [SMLMClustering documentation](https://github.com/JuliaSMLM/SMLMClustering.jl)
for the algorithm details and the full per-config parameter list.
