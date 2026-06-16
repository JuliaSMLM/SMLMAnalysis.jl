# dev/testdata — full pipeline runs on real data (papers-style layout)

One subfolder per real dataset, mirroring the papers-geometry / paper-genmab convention:
**symlinked data in, per-cell results written back to the NAS**.

```
dev/testdata/<dataset>/
  run.jl                     # full pipeline over every cell in the dataset
  README.md                  # this dataset's symlink targets + run command + status
  data/raw/<name>        ->  symlink to the NAS acquisition dir
  data/results/<name>    ->  symlink to a NAS results dir (outputs go BACK to the NAS)
                             one subfolder per cell: <condition>/Cell_NN/<step>/
```

## Conventions
- **Data is never copied or committed.** `data/raw/<name>` and `data/results/<name>` are
  NAS symlinks (see each dataset's README to recreate them); `dev/testdata/*/data/` is
  gitignored. Outputs land on the NAS via the `data/results` symlink, not in the repo.
- **Per-cell output**, like the papers: `data/results/<name>/<condition>/Cell_NN/<step>/`.
- **Run** from the repo's main project (the `dev/` env can lag; main is kept current):
  ```
  julia -t auto --project=/home/kalidke/julia_shared_dev/SMLMAnalysis \
      dev/testdata/<dataset>/run.jl
  ```
- New dataset → copy a subfolder, recreate its `data/raw` + `data/results` symlinks
  (per its README), repoint `run.jl`, tune params.

## Datasets
- **`sebastian_tnfar1_20260610/`** — Sebastian (gillette share), TNFaR1 dSTORM, MIC H5,
  256×256 / 4000-frame / 4-block cells. MeOH & PFA fixation titration. See its README.
