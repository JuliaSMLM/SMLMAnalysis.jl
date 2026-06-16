# Sebastian TNFaR1 dSTORM (acquired 2026-06-10)

Full-dataset SMLMAnalysis run on Sebastian's TNFaR1 acquisition (Gillette share),
using the papers-geometry / paper-genmab layout: symlinked data in, per-cell results
written **back to the NAS**.

## Layout (NAS symlinks — local setup, not tracked; recreate below)
```
data/raw/tnfar1        -> /mnt/nas/gillette/Sebastian/20260610_TNFaR1            (acquisition)
data/results/juliasmlm -> /mnt/nas/gillette/Sebastian/20260610_TNFaR1/juliasmlm  (results)
```
`run.jl` reads cells under `data/raw/tnfar1/<condition>/Cell_NN/Label_01/Data_*.h5`
and writes one subfolder per cell to `data/results/juliasmlm/<condition>/Cell_NN/<step>/`.

Recreate the symlinks (run from this folder):
```
ln -sfn /mnt/nas/gillette/Sebastian/20260610_TNFaR1 data/raw/tnfar1
mkdir -p /mnt/nas/gillette/Sebastian/20260610_TNFaR1/juliasmlm
ln -sfn /mnt/nas/gillette/Sebastian/20260610_TNFaR1/juliasmlm data/results/juliasmlm
```
(Results location matches paper-genmab's `Genmab/Data/juliasmlm` — co-located with the
data on the NAS. Repoint the `data/results/juliasmlm` symlink to move it.)

## Run
```
julia -t auto --project=/home/kalidke/julia_shared_dev/SMLMAnalysis \
    dev/testdata/sebastian_tnfar1_20260610/run.jl
```

## Data
MIC H5, 256×256, 4000 frames / 4 blocks, per-pixel sCMOS calibration. Conditions:
MeOH, MeOH_{100pM,500pM,2nM}, PFA, PFA_250nM (Cell_0N each).

## STATUS / TODO
- **Pixel size 97.8 nm is a placeholder** — confirm the Gillette acquisition scope's value.
- Detection params seeded from the hexabody run — **per-cell counts are low until tuned**.
- `dataset_mode=:continuous` (the 4 blocks are one movie, not registered stage positions).
- `MAX_GB=3.0` skips `PFA_250nM/Cell_02` (3.5 GB) + `Cell_03` (16.6 GB); raise to include them.
