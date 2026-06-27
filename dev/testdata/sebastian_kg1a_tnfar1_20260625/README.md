# Sebastian Kg1a TNFaR1 (1:250, ImStrF1) dSTORM (acquired 2026-06-25)

Full-dataset run of Sebastian's newest TNFaR1 acquisition using **papers-geometry's
exact SR/dSTORM config** (`SR_STEPS_BAGOL`, verbatim — see `run.jl` header for provenance).
Same AF647 dSTORM setup as the geometry cohort (camera + optics unchanged), so the config
transfers directly (confirmed by Keith via @papers-geometry-kitt, 2026-06-25).

## Data (NAS, read-only)
`/mnt/nas/gillette/Sebastian/20260625_Kg1a_TNFaR1-1-250_ImStrF1`
- MIC H5, 256×256, multi-block (one movie/cell), per-pixel sCMOS calibration embedded.
- Conditions: `1nM_F1_30ms`, `1nM_F1_40ms`, `1nM_F1_50ms`, `500pM_F1_40ms`
  (Cell_01..04 each; 500pM has Cell_01..02). 14 cells, ~0.57 GB each.

## Output
`dev/output/sebastian_kg1a_tnfar1_20260625/<condition>/Cell_NN/<NN_step>/`
(under the repo so the 5nm render is nav-previewable; shared FS so descent writes / kitt reads).
Per-cell markers in `dev/output/.../_markers/` (`*.done` holds the repo-relative 5nm-render path).

## Run (on descent, GPU 0)
```
CUDA_VISIBLE_DEVICES=0 julia -t auto \
  --project=/home/kalidke/julia_shared_dev/SMLMAnalysis \
  dev/testdata/sebastian_kg1a_tnfar1_20260625/run.jl
```
Optional cell filter: append `<condition>/Cell_NN` substrings, e.g. `... run.jl 1nM_F1_30ms 500pM`.
Re-running skips cells whose `_markers/<tag>.done` exists.

## Config = geometry SR_STEPS_BAGOL (01..09)
detectfit(box9, σ0.130, min200, XYNBS×20) → filter(phot≥200, prec≤10nm, p≥1e-6, σ∈[100,150]nm)
→ frameconnect(gap5, 5σ, clamp_k) → drift(deg3, **:registered**, iterative, shift_scale .050)
→ 05 gaussian 20x inferno (**5nm render**) → 06 hist 10x turbo/time → 07 circle 50x turbo/time
→ 08 BaGoL(se_adjust=:auto, 6000/2000, μ10 shape1 learn, maxpart500, overlap .025, post_px 0)
→ 09 gaussian 50x inferno (post-BaGoL MAP-N).

## Notes / TODO
- **drift `:registered`** is geometry's value (exact-match). Prior 20260610 TNFaR1 harness used
  `:continuous` (the MIC blocks are one movie). Flip `SR_DRIFT_DATASET_MODE` in `run.jl` if drift
  looks wrong on this block structure.
- Camera uses each cell's OWN embedded sCMOS calibration (same camera as geometry, fresher cal).
