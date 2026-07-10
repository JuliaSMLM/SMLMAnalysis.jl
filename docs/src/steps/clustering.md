```@meta
CurrentModule = SMLMAnalysis
```

# [Clustering](@id clustering-step)

A super-resolution point cloud is rarely uniform: receptors aggregate, complexes
assemble, and the background scatters at random. The clustering step asks two
distinct questions of that cloud — *which points belong to the same structure?*
and *is there any structure at all?* — and is backed by
[SMLMClustering](https://github.com/JuliaSMLM/SMLMClustering.jl). It exposes two
families of config, selected entirely by **type**:

- **Cluster-labeling** (`DBSCANConfig`, `HDBSCANConfig`, `VoronoiConfig`,
  `HierarchicalConfig`) — assigns a cluster id to every localization and returns
  a **labeled SMLD**.
- **Spatial statistics** (`HopkinsConfig`, `VoronoiDensityConfig`) — computes a
  read-only summary and returns the **SMLD unchanged**.

```julia
analyze(smld, DBSCANConfig(eps_nm = 50.0, min_points = 5))  # → (labeled_smld, StepInfo)
analyze(smld, HopkinsConfig())                              # → (smld, StepInfo)
```

## When to use / prerequisites

- Run on a `BasicSMLD` of localizations, normally late in the pipeline — after
  [detection/fitting](@ref "Detection & Fitting"), a [quality filter](@ref
  "Quality Filter"), and usually [frame connection](@ref "Frame Connection") and
  [drift correction](@ref "Drift Correction"), so clusters reflect real structure
  rather than blinking or drift.
- Use a **labeling** backend when you want to segment and count structures, then
  [render](@ref "Rendering") or export them colored by id. Use a **statistics**
  backend (`HopkinsConfig`) to check clustering tendency *before* committing to a
  labeling pass, or to attach a per-emitter density feature for your own
  thresholding.
- Most length parameters are given in **nm** for convenience even though
  `BasicSMLD` positions are in µm; each backend documents its own units upstream.

## Inputs, returns & artifacts

- **Input:** the current `smld`.
- **Returns:** `(smld_out, StepInfo)`.
  - *Labeling configs:* `smld_out` is a deep-copied, labeled SMLD — the input is
    **not mutated**. Each `emitter.id` carries the cluster id: `0` = noise /
    unclustered, `1..K` = clusters. `StepInfo.info` is a `ClusterInfo`.
  - *Statistics configs:* `smld_out` is the **same** SMLD, untouched. The result
    scalar and any per-emitter/per-dataset vectors live on
    `StepInfo.info::ClusterStatisticsInfo` (`info.statistic`, `info.extras`).
- **Artifacts** (when `outdir` is set, written into `{NN}_{backend}/` — e.g.
  `03_dbscan/`, `04_hopkins/`): at `STANDARD` verbosity, `config.toml` and
  `info.toml`. Labeling backends additionally checkpoint the labeled SMLD to
  `smld_clustered.jld2`, but only at `Checkpoint.ALL` (clustering is not treated
  as an "expensive" step, so this is *not* written by the default
  `Checkpoint.EXPENSIVE`).

The step directory name is the backend name itself (`dbscan`, `hdbscan`,
`voronoi`, `hierarchical`, `hopkins`, `voronoidensity`), derived from the config
type — so a pipeline that clusters then renders the labels reads cleanly on disk.

## Concept

The two `analyze()` methods are thin wrappers over SMLMClustering's two verbs.
Labeling dispatches to `SMLMClustering.cluster`, which discovers `K` clusters and
writes an integer label onto each `emitter.id`; because the label lives on `id`, a
subsequent [render](@ref "Rendering") (or any downstream step) sees the clustered
SMLD with no extra plumbing. Statistics dispatches to
`SMLMClustering.cluster_statistics`, which is read-only — `HopkinsConfig` returns
a clustering-tendency scalar `H`, Voronoi density returns per-emitter densities.

The backends span density methods (DBSCAN, HDBSCAN), tessellation density (Voronoi
/ SR-Tesseler), and agglomerative linkage (Hierarchical). The algorithm math, the
full backend catalog, and per-backend tuning live in the
[SMLMClustering documentation](https://github.com/JuliaSMLM/SMLMClustering.jl) —
this page is the pipeline-integration view, not a re-documentation of it.

## Configuration

All cluster configs are owned by SMLMClustering (the "upstream owns the config"
idiom); SMLMAnalysis re-exports the common ones as `const` aliases. Dispatch is
on the **abstract** types `AbstractClusterConfig` / `AbstractStatisticsConfig`,
so *any* SMLMClustering backend of those families composes as a step — including
ones SMLMAnalysis does not alias (reach them as `SMLMClustering.XConfig`).

Backend selector and the field you tune first:

| Config | Family | Key field(s) |
|---|---|---|
| `DBSCANConfig` | density `ε`-neighborhood | `eps_nm` (neighborhood radius), `min_points` |
| `HDBSCANConfig` | hierarchical density | `min_points` (core-distance *k*), `min_cluster_size` |
| `VoronoiConfig` | tessellation density (SR-Tesseler) | `density_factor`, `min_points` |
| `HierarchicalConfig` | agglomerative linkage | `cut_threshold` *or* `n_clusters`, `linkage` |
| `HopkinsConfig` | statistic (tendency `H`) | `n_samples`, `random_repeats`, `seed` |
| `VoronoiDensityConfig` | statistic (per-emitter density) | — (uses `use_3d` / `per_dataset` only) |

Fields shared by the labeling backends: `min_points` (minimum cluster size; for
HDBSCAN this is instead the core-distance neighbor count), `per_dataset = true`
(cluster each dataset independently — ids are then local to a dataset), `use_3d`
(2D vs 3D; Voronoi is 2D only), and `remove_unclustered = false` (when `true`,
noise points are dropped from the returned SMLD rather than labeled `0`).

```julia
# Density-based clustering, then render colored by cluster id
(labeled, info) = analyze(smld, DBSCANConfig(eps_nm = 50.0, min_points = 5))
info.info.n_clusters                       # number of clusters found
(img, _) = analyze(labeled, RenderConfig(zoom = 20))   # labels flow through emitter.id

# Read-only: is the data clustered at all?
(_, h) = analyze(smld, HopkinsConfig())
h.info.statistic                           # Hopkins H (≈0.5 random, →1 clustered)
```

See the SMLMClustering documentation for the complete field list, defaults, and
the additional backends.

## Output & interpretation

`StepInfo.summary` reports the headline numbers. For **labeling** backends
(`ClusterInfo`):

| field | meaning |
|---|---|
| `algorithm` | backend used (`:dbscan`, `:voronoi`, …) |
| `n_locs_in` | localizations entering the step |
| `n_clustered` | localizations assigned to a cluster (`id ≥ 1`) |
| `n_noise` | localizations left as noise (`id = 0`) |
| `n_clusters` | number of clusters `K` discovered |

For **statistics** backends (`ClusterStatisticsInfo`): `algorithm`, `n_locs_in`,
`statistic_name`, and the scalar `statistic`.

Sanity checks: `n_clustered + n_noise` should equal `n_locs_in` (unless
`remove_unclustered = true`). A Hopkins `H` near `0.5` means the cloud is
consistent with spatial randomness — no clusters to find — while `H` approaching
`1` indicates genuine aggregation. If DBSCAN returns one giant cluster, `eps_nm`
is too large; if everything is noise, it is too small (or `min_points` too high).

## Notes & caveats

- **Labels are not stable across runs.** `1..K` are discovered per run and only
  meaningful within it; with `per_dataset = true` the pair `(dataset, id)` uniquely
  identifies a cluster. A later labeling step reusing `emitter.id` overwrites an
  earlier one.
- **All clustering backends are re-exported.** `DBSCANConfig`, `HDBSCANConfig`,
  `HierarchicalConfig`, and `VoronoiConfig` each have a `const` alias in
  SMLMAnalysis, so they are available unqualified after `using SMLMAnalysis`. A
  backend without an alias would need to be constructed qualified as
  `SMLMClustering.<Name>(...)`.
- **Statistics steps change nothing downstream.** They are observational; place
  them anywhere without affecting the threaded SMLD.

## References

- M. Ester, H.-P. Kriegel, J. Sander, X. Xu. "A Density-Based Algorithm for
  Discovering Clusters in Large Spatial Databases with Noise." *Proc. 2nd Int.
  Conf. on Knowledge Discovery and Data Mining (KDD-96)*, 226–231 (1996).
- F. Levet, E. Hosy, A. Kechkar, C. Butler, A. Beghin, D. Choquet, J.-B. Sibarita.
  "SR-Tesseler: a method to segment and quantify localization-based
  super-resolution microscopy data." *Nature Methods* **12**, 1065–1071 (2015).
  [doi:10.1038/nmeth.3579](https://doi.org/10.1038/nmeth.3579)

See the [SMLMClustering documentation](https://github.com/JuliaSMLM/SMLMClustering.jl)
for the algorithms, the full backend catalog, and all configuration options.
