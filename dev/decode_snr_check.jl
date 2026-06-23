#
# dev/decode_snr_check.jl — model-free branch test (dumps 1+2 of the standard-DECODE A/B).
#
# Renders a single emitter with the learned VECTOR PSF (+ a plain Parseval SCALAR baseline)
# at the TRAINING conditions (photons=1000, bg=5) and z = 0, ±0.5µm — the Decode z-range that
# DOMINATES training: a vector PSF defocuses hard, so most training spots are the broad ones,
# not the sharp focus spot. Reports ABSOLUTE peak-pixel shot-noise SNR (the detectability
# decider — not a ratio to a maybe-still-overbright scalar) + peak/sum energy fraction + FWHM.
#
# Routes the over-detection debug WITHOUT a 40-min train:
#   flat / low-abs-SNR ............... branch A   (rendering / PSF-sampling / photon budget)
#   sharp clean peak across z ........ branch B   (bump loss_count weight)
#   clean but BROAD at defocus ....... branch B'  (concentrate p̂ + lean on Δx/Δy offset channels)
#
# No Reactant/cuDNN — pure integrate_pixels (CPU). Run: julia --project=dev dev/decode_snr_check.jl
#
include(joinpath(@__DIR__, "deepfit_sim_common.jl"))   # lightweight: constants only
using MicroscopePSFs, SMLMData, CairoMakie, Printf

const _BG = 5.0
const _PH = 1000.0
const _SR = 25          # render ROI (px) — wide enough to capture the defocused vector wings

function _spot(psf, z)
    cam  = IdealCamera(_SR, _SR, PIXEL_SIZE)
    c    = _SR * PIXEL_SIZE / 2
    unit = Float64.(integrate_pixels(psf, cam, Emitter3D(c, c, Float64(z), 1.0)))  # unit-photon PSF
    spot = _PH .* unit
    peak = maximum(spot); tot = sum(spot)
    (spot = spot, peak = peak, tot = tot, frac = peak / max(tot, 1e-9),
     snr = peak / sqrt(peak + _BG), fwhm = 2 * sqrt(count(>=(peak / 2), spot) / π))
end

function main()
    vpsf = MicroscopePSFs.load_psf(get(ENV, "DECODE_PSF", PSF_PATH))   # DECODE_PSF env to compare an alternate PSF (else the stage-1 psf.h5)
    spsf = ScalarPSF(1.4, 0.68, 1.516)              # plain Parseval scalar baseline (NA, λ, n)
    zs   = [0.0, 0.5, -0.5]
    @printf "photons=%.0f  bg=%.0f  roi=%dx%d  pixel=%.2fum\n" _PH _BG _SR _SR PIXEL_SIZE
    @printf "%-7s %6s %9s %9s %9s %8s %8s\n" "PSF" "z(um)" "peak_ph" "sum_ph" "peak/sum" "FWHM_px" "absSNR"
    res = Dict()
    for (nm, psf) in (("VECTOR", vpsf), ("SCALAR", spsf)), z in zs
        s = _spot(psf, z); res[(nm, z)] = s
        @printf "%-7s %+5.1f %9.1f %9.1f %9.3f %8.2f %8.1f\n" nm z s.peak s.tot s.frac s.fwhm s.snr
    end
    fig = Figure(size = (3 * 250, 2 * 250 + 46))
    Label(fig[0, 1:3], "single emitter: VECTOR (learned) vs SCALAR @ z=0,±0.5µm — photons=$(Int(_PH)) bg=$(Int(_BG)); absSNR=peak/√(peak+bg)", fontsize = 12)
    for (ri, nm) in enumerate(("VECTOR", "SCALAR")), (ci, z) in enumerate(zs)
        s = res[(nm, z)]
        ax = Axis(fig[ri, ci], aspect = DataAspect(), yreversed = true,
                  title = "$nm  z=$(z)µm   peak=$(round(Int, s.peak))ph  SNR=$(round(s.snr, digits = 1))  pk/sum=$(round(s.frac, digits = 3))", titlesize = 9)
        heatmap!(ax, s.spot', colormap = :inferno); hidedecorations!(ax)
    end
    CairoMakie.save(joinpath(OUTDIR, "decode_snr_check.png"), fig)
    println("\nwrote ", joinpath(OUTDIR, "decode_snr_check.png"))
end

main()
