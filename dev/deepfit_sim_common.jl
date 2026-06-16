# Shared constants for the deepfit simulated end-to-end chain.
# Stage 1 (dev/learn_psf_sim.jl, CUDA-hidden) and stages 2-5 (dev/deepfit_sim.jl, GPU)
# both include this so the spatial scale + paths line up. psf.h5 is the bridge.
const OUTDIR     = joinpath(@__DIR__, "deepfit_sim_out")
const PIXEL_SIZE = 0.1      # µm/px — matches the deepfit camera + PSF sampling
const ROI        = 16       # bead ROI / psflearning roi_size (px)
const FOV        = 64       # structure-sim FOV (px); divisible by 4 for the depth-2 U-Net
const PSF_PATH   = joinpath(OUTDIR, "01_psflearning", "psf.h5")
