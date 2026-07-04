```@meta
CurrentModule = SMLMAnalysis
```

# [References](@id references-page)

The primary literature for the methods SMLMAnalysis orchestrates, grouped by
step. Each [step page](@ref "Pipeline Steps: Overview") repeats its own reference
and links to the backing package's documentation for the algorithm details.

## Detection & fitting

- **Detection (SMLMBoxer).** F. Huang, T. M. P. Hartwich, F. E. Rivera-Molina,
  *et al.* "Video-rate nanoscopy using sCMOS camera-specific single-molecule
  localization algorithms." *Nature Methods* **10**, 653–658 (2013).
  [doi:10.1038/nmeth.2488](https://doi.org/10.1038/nmeth.2488)
- **MLE fitting (GaussMLE).** C. S. Smith, N. Joseph, B. Rieger, K. A. Lidke.
  "Fast, single-molecule localization that achieves theoretically minimum
  uncertainty." *Nature Methods* **7**, 373–375 (2010).
  [doi:10.1038/nmeth.1449](https://doi.org/10.1038/nmeth.1449)
- **Localization precision / CRLB.** The reported `σ` are the *exact* CRLB from
  the MLE Fisher information — Smith *et al.* 2010 (above), extended to per-pixel
  sCMOS noise by Huang *et al.* 2013 (above) — not an analytical approximation.
- **sCMOS noise model.** Huang *et al.* 2013 (above).
- **Multi-emitter fitting & goodness-of-fit (`pvalue`).** F. Huang, S. L.
  Schwartz, J. M. Byars, K. A. Lidke. "Simultaneous multiple-emitter fitting for
  single molecule super-resolution imaging." *Biomedical Optics Express* **2**,
  1377–1393 (2011).
  [doi:10.1364/BOE.2.001377](https://doi.org/10.1364/BOE.2.001377)

## Frame connection

- D. J. Schodt, K. A. Lidke. "Spatiotemporal Clustering of Repeated
  Super-Resolution Localizations via Linear Assignment Problem." *Frontiers in
  Bioinformatics* **1**, 724325 (2021).
  [doi:10.3389/fbinf.2021.724325](https://doi.org/10.3389/fbinf.2021.724325)

## Drift correction

- J. Cnossen, T. J. Cui, C. Joo, C. Smith. "Drift correction in localization
  microscopy using entropy minimization." *Optics Express* **29**, 27961–27974
  (2021). [doi:10.1364/OE.426620](https://doi.org/10.1364/OE.426620)
- M. J. Wester, D. J. Schodt, H. Mazloom-Farsibaf, M. Fazel, S. Pallikkuth,
  K. A. Lidke. "Robust, fiducial-free drift correction for super-resolution
  imaging." *Scientific Reports* **11**, 23672 (2021).
  [doi:10.1038/s41598-021-02850-7](https://doi.org/10.1038/s41598-021-02850-7)

## Bayesian grouping

- M. Fazel, *et al.* "High-Precision Estimation of Emitter Positions using
  Bayesian Grouping of Localizations." *Nature Communications* **13**, 7152
  (2022).
  [doi:10.1038/s41467-022-34894-2](https://doi.org/10.1038/s41467-022-34894-2)

## Clustering

- **DBSCAN.** M. Ester, H.-P. Kriegel, J. Sander, X. Xu. "A Density-Based
  Algorithm for Discovering Clusters in Large Spatial Databases with Noise."
  *Proc. 2nd Int. Conf. on Knowledge Discovery and Data Mining (KDD-96)*,
  226–231 (1996).
- **Voronoi / SR-Tesseler.** F. Levet, *et al.* "SR-Tesseler: a method to
  segment and quantify localization-based super-resolution microscopy data."
  *Nature Methods* **12**, 1065–1071 (2015).
  [doi:10.1038/nmeth.3579](https://doi.org/10.1038/nmeth.3579)

See the [SMLMClustering documentation](https://github.com/JuliaSMLM/SMLMClustering.jl)
for the full set of backends and their individual references.

## Cross-correlation (multi-channel)

- **Pair correlation, theory.** B. D. Ripley. "Modelling spatial patterns."
  *Journal of the Royal Statistical Society B* **39**, 172–212 (1977).
- **Pair correlation in SMLM.** P. Sengupta, T. Jovanovic-Talisman, D. Skoko,
  *et al.* "Probing protein heterogeneity in the plasma membrane using PALM and
  pair correlation analysis." *Nature Methods* **8**, 969–975 (2011).
  [doi:10.1038/nmeth.1704](https://doi.org/10.1038/nmeth.1704)

## Methods native to SMLMAnalysis

The [Quality Filter](@ref), [Intensity Filter](@ref), and [Density Filter](@ref)
steps are implemented in SMLMAnalysis. Their methods (Poisson upper-tail
multi-emitter rejection against an estimated excitation field; neighbor-count
density filtering with automatic threshold selection) are described in full on
their step pages.

## How to cite

If you use SMLMAnalysis in your research, please cite the package together with
the primary references for the specific methods your pipeline used (the steps
above). A package citation entry (`CITATION.bib`) will accompany the registered
release.
