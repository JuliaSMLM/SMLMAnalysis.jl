```@meta
CurrentModule = SMLMAnalysis
```

# Multi-Dataset Acquisitions

Most real acquisitions are not one continuous movie. A long experiment is split
into several blocks — multiple cells, several ROIs on one coverslip, or a single
long movie broken into segments. SMLMAnalysis treats these as **datasets** that
live inside one `BasicSMLD`, processed together but tracked separately. This page
covers what a dataset is, why the split matters, and the three ways to feed
multi-dataset data into [Detection & Fitting](@ref "Detection & Fitting").

## What a "dataset" is

A **dataset** is an independent acquisition block within a single localization
set — a cell, an ROI, or a registered segment of one long movie. Every emitter
records which block it came from in its `dataset` field (see
[Data Model & Provenance](@ref)), and the `BasicSMLD` knows how many blocks it
holds. Datasets are not user-numbered: their boundaries come from the *shape of
the data you supply* (one image stack, one file, or one frame range per dataset).

## Why split into datasets

- **Memory efficiency.** [Detection & Fitting](@ref "Detection & Fitting") loops
  over datasets one at a time — detect, fit, append, free. With the file-based
  modes below, only one dataset's images are resident at once, so arbitrarily
  large acquisitions fit in memory.
- **Per-dataset frame numbering.** Each emitter's `frame` is relative to its own
  dataset (`1` to `n_frames_per_dataset`), **not** a global frame index. This is
  required by [Drift Correction](@ref): the Legendre-polynomial drift basis
  normalizes each dataset's time axis to `[-1, 1]`, so frame counting must restart
  at each dataset boundary.
- **Dataset tracking.** Because `emitter.dataset` is carried through every step,
  downstream operations (notably drift correction's per-dataset trajectories and
  inter-dataset alignment) know exactly which block each localization belongs to.

## SMLD structure

After detection and fitting, the combined `BasicSMLD` encodes the split in two
fields:

| Field | Meaning |
|-------|---------|
| `smld.n_frames` | frames **per dataset** (not the global total) |
| `smld.n_datasets` | number of datasets in the set |
| `emitter.dataset` | which dataset (1-based) an emitter belongs to |
| `emitter.frame` | frame index **within** that dataset (1-based) |

So an acquisition of 4 datasets × 2000 frames has `n_frames == 2000` and
`n_datasets == 4`; the total of 8000 frames is `n_frames * n_datasets`.

!!! note "n_frames is per-dataset"
    If you expect `n_frames` to be the grand total and see a smaller number,
    that is correct — it is frames *per dataset*. Multiply by `n_datasets` for the
    total.

## Supplying multi-dataset data

`DetectFitConfig` resolves datasets from whichever data source you provide. There
are three modes:

| Mode | How datasets are defined |
|------|--------------------------|
| In-memory | one 3D image stack per element of a `Vector` |
| Single file | MIC blocks (auto), or explicit `dataset_frames` ranges |
| Multiple files | one file per element of `paths` |

### In-memory: a vector of image stacks

Pass a `Vector` of `H × W × F` arrays — each element is one dataset:

```julia
image_stacks = [stack1, stack2, stack3, stack4]   # 4 datasets
(result, info) = analyze(image_stacks, config)

# Step-by-step is the same input:
(smld, info) = analyze(image_stacks,
    DetectFitConfig(camera = cam, boxer = BoxerConfig(boxsize = 9)))
```

A single 3D array is wrapped automatically into a one-dataset set, so
single-acquisition code needs no special handling.

### Single file split by frame range or block

Point `DetectFitConfig` at one H5 file and let it carve the datasets. MIC-format
files expose acquisition **blocks** that auto-detect as datasets; for other
formats, give explicit `dataset_frames` ranges:

```julia
# MIC blocks become datasets automatically
config = AnalysisConfig(camera = cam, steps = [
    DetectFitConfig(path = "data.h5", h5_format = :mic,
                    boxer = BoxerConfig(boxsize = 9)),
    RenderConfig(zoom = 20),
])
(result, info) = analyze(config)          # no data argument — loaded from file

# Or split one stack into explicit per-dataset frame ranges
DetectFitConfig(path = "movie.h5",
                dataset_frames = [1:2000, 2001:4000, 4001:6000])
```

`h5_format` accepts `:auto` (default), `:smart`, or `:mic`. Each source is loaded,
processed, and freed before the next — the memory-efficient path for large files.

### Multiple files, one per dataset

Give a `paths` vector to treat each file as its own dataset:

```julia
DetectFitConfig(paths = ["d1.h5", "d2.h5", "d3.h5", "d4.h5"],
                boxer = BoxerConfig(boxsize = 9))
```

### Selecting a subset of datasets

`datasets` picks a subset of the resolved sources — useful for re-running a few
blocks. It applies uniformly across all three modes, and the selected blocks are
reindexed to a contiguous `1:length(datasets)` in the output:

```julia
DetectFitConfig(path = "data.h5", h5_format = :mic, datasets = [1, 2, 5])
```

The original source indices are preserved in `DetectFitInfo.selected_source_indices`.

## Connecting to drift correction

Multi-dataset structure exists largely to serve [Drift Correction](@ref), whose
two modes map directly onto how the datasets relate:

- **`dataset_mode = :registered`** (default) — for spatially registered blocks
  (the stage returns to roughly the same position between datasets). Each dataset
  gets its own intra-dataset drift polynomial, then inter-dataset alignment via
  entropy optimization over the spatial overlap.
- **`dataset_mode = :continuous`** — for one long acquisition split into segments,
  where drift accumulates continuously across the dataset boundaries.

See [Drift Correction](@ref) for the full mode reference, chunking guidance, and
quality settings. For how datasets thread through the rest of the pipeline, see
[The Pipeline Model](@ref) and [Running a Pipeline](@ref).
