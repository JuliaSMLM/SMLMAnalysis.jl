using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Random

const SMLM_TEST_FULL = lowercase(get(ENV, "SMLM_TEST_FULL", "false")) in ("true", "1", "yes")

@testset "fast" begin
    @testset "Types" begin
        # Test AnalysisInfo constructor
        info = AnalysisInfo()
        @test info.elapsed_s == 0.0
        @test isempty(info.steps)
        @test isempty(info.step_infos)

        # Test AnalysisInfo with data
        steps = Dict{Symbol, Any}(:test => (a=1, b=2))
        info = AnalysisInfo(1.5, steps, StepInfo[])
        @test info.elapsed_s == 1.5
        @test info.steps[:test].a == 1

        # Test StepInfo with typed info
        cfg = FilterConfig()
        filter_info = FilterInfo(100, 80, 0.5)
        step_info = StepInfo(1, cfg, 0.5, Dict{Symbol,Any}(:n_before => 100); info=filter_info)
        @test step_info.info !== nothing
        @test step_info.info isa FilterInfo
        @test step_info.info.n_before == 100
        @test step_info.info.n_after == 80
        @test step_info.elapsed_s == 0.5

        # Test StepInfo without info
        step_info2 = StepInfo(2, cfg, 0.3, Dict{Symbol,Any}())
        @test step_info2.info === nothing

        # Test native info structs
        # Back-compat 7-arg constructor (defaults selected_source_indices to nothing)
        di = DetectFitInfo([], [], 2, 1000, 950, 5000, 1.5)
        @test di.n_datasets == 2
        @test di.n_rois == 1000
        @test di.n_fits == 950
        @test di.selected_source_indices === nothing

        # Full 8-arg constructor with provenance
        di_sel = DetectFitInfo([], [], 3, 500, 450, 1000, 0.5, [1, 3, 5])
        @test di_sel.selected_source_indices == [1, 3, 5]
        @test di_sel.n_datasets == 3

        dfi = DensityFilterInfo(1000, 800, 5, 0.3)
        @test dfi.n_before == 1000
        @test dfi.threshold == 5
    end

    @testset "analyze dispatch" begin
        # Verify analyze() dispatch methods exist for each step config type
        @test hasmethod(analyze, Tuple{Vector{<:AbstractArray{<:Real,3}}, DetectFitConfig})
        @test hasmethod(analyze, Tuple{AbstractArray{<:Real,3}, DetectFitConfig})
        @test hasmethod(analyze, Tuple{DetectFitConfig})
        @test hasmethod(analyze, Tuple{BasicSMLD, FilterConfig})
        @test hasmethod(analyze, Tuple{BasicSMLD, FrameConnectConfig})
        @test hasmethod(analyze, Tuple{BasicSMLD, DriftConfig})
        @test hasmethod(analyze, Tuple{BasicSMLD, DensityFilterConfig})
        @test hasmethod(analyze, Tuple{BasicSMLD, RenderConfig})

        # Verify old step function names are not exported
        @test !isdefined(Main, :detectfit)
        @test !isdefined(Main, :filter_step)
        @test !isdefined(Main, :frameconnect_step)
        @test !isdefined(Main, :driftcorrect_step)
        @test !isdefined(Main, :densityfilter_step)
        @test !isdefined(Main, :render_step)

        # DetectFitConfig camera field
        cfg = DetectFitConfig()
        @test cfg.camera === nothing
        cam = IdealCamera(64, 64, 0.1)
        cfg2 = DetectFitConfig(camera=cam, boxer=BoxerConfig(boxsize=7))
        @test cfg2.camera === cam
        @test cfg2.boxer.boxsize == 7
        @test cfg2.fitter.psf_model isa GaussianXYNBS

        # DetectFitConfig.datasets selection field
        @test cfg.datasets === nothing                          # default is no selection
        cfg_range = DetectFitConfig(datasets=1:19)
        @test cfg_range.datasets == 1:19
        @test cfg_range.datasets isa UnitRange{Int}
        cfg_sparse = DetectFitConfig(datasets=[1, 2, 3, 5, 7])
        @test cfg_sparse.datasets == [1, 2, 3, 5, 7]
        @test cfg_sparse.datasets isa Vector{Int}

        # _select_sources: pass-through when nothing, bounds-checked otherwise
        src = [(i=j,) for j in 1:5]
        @test SMLMAnalysis._select_sources(src, nothing) === src
        @test SMLMAnalysis._select_sources(src, [1, 3, 5]) == [src[1], src[3], src[5]]
        @test SMLMAnalysis._select_sources(src, 2:4) == src[2:4]
        @test_throws ErrorException SMLMAnalysis._select_sources(src, [1, 6])
        @test_throws ErrorException SMLMAnalysis._select_sources(src, [0, 1])
    end

    @testset "CalibrationConfig re-export" begin
        # CalibrationConfig is re-exported from SMLMFrameConnection
        @test CalibrationConfig === SMLMFrameConnection.CalibrationConfig
        @test CalibrationResult === SMLMFrameConnection.CalibrationResult

        # CalibrationConfig can be nested in FrameConnectConfig
        cal_cfg = CalibrationConfig(clamp_k_to_one=true)
        fc_cfg = FrameConnectConfig(max_frame_gap=5, calibration=cal_cfg)
        @test fc_cfg.calibration !== nothing
        @test fc_cfg.calibration.clamp_k_to_one == true

        # Default is nothing (no calibration)
        fc_default = FrameConnectConfig()
        @test fc_default.calibration === nothing
    end

    @testset "Info struct subtypes" begin
        # All info structs should be AbstractSMLMInfo subtypes
        @test StepInfo <: AbstractSMLMInfo
        @test DetectFitInfo <: AbstractSMLMInfo
        @test FilterInfo <: AbstractSMLMInfo
        @test DensityFilterInfo <: AbstractSMLMInfo
        @test CompositeRenderInfo <: AbstractSMLMInfo
        @test CrossAlignInfo <: AbstractSMLMInfo
        @test AnalysisInfo <: AbstractSMLMInfo
    end

    @testset "Multi-target step types" begin
        # Type hierarchy
        @test AbstractMultiTargetStep <: AbstractSMLMConfig
        @test CompositeRenderConfig <: AbstractMultiTargetStep
        @test CrossAlignConfig <: AbstractMultiTargetStep

        # CompositeRenderConfig defaults
        cr = CompositeRenderConfig()
        @test cr.strategy isa GaussianRender
        @test cr.zoom == 20.0
        @test cr.colors === nothing
        @test cr.clip_percentile === :auto
        @test cr.normalize_each === nothing
        @test cr.scalebar == true
        @test cr.scalebar_position == :br

        # CompositeRenderConfig with custom fields
        cr2 = CompositeRenderConfig(strategy=HistogramRender(), zoom=10.0, colors=[:red, :blue])
        @test cr2.strategy isa HistogramRender
        @test cr2.zoom == 10.0
        @test cr2.colors == [:red, :blue]

        # CrossAlignConfig defaults
        ca = CrossAlignConfig()
        @test ca.method == :entropy
        @test ca.maxn == 100
        @test ca.histbinsize == 0.05

        # CrossAlignConfig custom
        ca2 = CrossAlignConfig(method=:fft, maxn=50)
        @test ca2.method == :fft
        @test ca2.maxn == 50

        # step_name dispatch
        @test SMLMAnalysis.step_name(cr) == "compositerender"
        @test SMLMAnalysis.step_name(ca) == "crossalign"

        # analyze dispatch methods exist for multi-target steps
        @test hasmethod(analyze, Tuple{Vector{<:BasicSMLD}, CompositeRenderConfig})
        @test hasmethod(analyze, Tuple{Vector{<:BasicSMLD}, CrossAlignConfig})

        # MultiTargetConfig with steps vector
        mt = MultiTargetConfig(
            labels=[:A, :B],
            steps=[
                CompositeRenderConfig(zoom=20.0),
                CrossAlignConfig(),
                CompositeRenderConfig(zoom=10.0, strategy=HistogramRender()),
            ],
            outdir="/tmp/test_mt",
        )
        @test length(mt.steps) == 3
        @test mt.steps[1] isa CompositeRenderConfig
        @test mt.steps[2] isa CrossAlignConfig
        @test mt.steps[3] isa CompositeRenderConfig
        @test mt.colors == [:cyan, :magenta]

        # AlignConfig/AlignInfo re-exports
        @test AlignConfig === SMLMDriftCorrection.AlignConfig
        @test AlignInfo === SMLMDriftCorrection.AlignInfo
    end

    @testset "crosscorr g(r)" begin
        # Exercises the internal _compute_crosscorr for the three fixed defects:
        # (a) bin-width consistency when r_max is NOT an integer multiple of dr,
        # (b) exactly-coincident cross-channel pairs (dist==0) are counted, not
        #     dropped, and (c) CSR normalization. CPU-only; no GPU needed.
        T = Float64

        # Minimal emitter builder — only x/y matter to _compute_crosscorr; every
        # other field is a constant dummy. Layout mirrors the Emitter2DFitSigma
        # constructor used in the HDF5 round-trip testset above.
        mkem(x, y, i) = GaussMLE.Emitter2DFitSigma{T}(
            x, y, 1000.0, 5.0, 0.13,
            0.01, 0.012, 0.003, 20.0, 0.5, 0.002,
            0.4, 1, 1, 0, i)
        mksmld(xs, ys, cam) = BasicSMLD(
            [mkem(xs[i], ys[i], i) for i in eachindex(xs)],
            cam, 1, 1, Dict{String,Any}())

        cam = IdealCamera(64, 64, 0.1)   # 6.4 × 6.4 μm FOV
        xlo, xhi = first(cam.pixel_edges_x), last(cam.pixel_edges_x)
        ylo, yhi = first(cam.pixel_edges_y), last(cam.pixel_edges_y)
        nx() = xlo + rand() * (xhi - xlo)
        ny() = ylo + rand() * (yhi - ylo)

        # (a) Non-divisor consistency: r_max=1.0 is not a multiple of dr=0.03.
        Random.seed!(7)
        cfg_nd = CrossCorrConfig(r_max=1.0, dr=0.03)
        sa = mksmld([nx() for _ in 1:50], [ny() for _ in 1:50], cam)
        sb = mksmld([nx() for _ in 1:50], [ny() for _ in 1:50], cam)
        (r_nd, g_nd, area_nd) = SMLMAnalysis._compute_crosscorr(sa, sb, cfg_nd)
        # r_centers are spaced by exactly dr — one shared width for counting,
        # annulus areas, and the r-axis.
        @test all(isapprox.(diff(r_nd), cfg_nd.dr; atol=1e-9))
        # Last bin's outer edge (= last center + dr/2) reaches at least r_max.
        @test (last(r_nd) + cfg_nd.dr / 2) >= cfg_nd.r_max - 1e-9
        @test length(r_nd) == ceil(Int, cfg_nd.r_max / cfg_nd.dr)

        # (b) Zero-distance kept: A and B share exactly-coincident points (plus
        # random background). The coincident pairs must elevate the first bin.
        Random.seed!(11)
        ncoinc = 30
        # Coincident locations on a diagonal, spacing ≈ 0.27 μm ≫ dr, so only the
        # exact same-location A/B pairs (dist==0) fall into bin 1.
        cxs = [0.4 + 0.19k for k in 0:(ncoinc - 1)]
        cys = [0.4 + 0.19k for k in 0:(ncoinc - 1)]
        bg_ax = [nx() for _ in 1:100]; bg_ay = [ny() for _ in 1:100]
        bg_bx = [nx() for _ in 1:100]; bg_by = [ny() for _ in 1:100]
        sa_z = mksmld(vcat(cxs, bg_ax), vcat(cys, bg_ay), cam)
        sb_z = mksmld(vcat(cxs, bg_bx), vcat(cys, bg_by), cam)
        (r_z, g_z, _) = SMLMAnalysis._compute_crosscorr(sa_z, sb_z, CrossCorrConfig())
        @test g_z[1] > 1   # coincident cross-channel pairs are counted, not dropped
        @test g_z[1] > 2   # and clearly elevated above the CSR background

        # (c) CSR sanity: two independent uniform-random channels give g(r) ≈ 1
        # across mid-range bins (away from r→0 noise and the far edge).
        Random.seed!(1234)
        npts = 3000
        ax = [nx() for _ in 1:npts]; ay = [ny() for _ in 1:npts]
        bx = [nx() for _ in 1:npts]; by = [ny() for _ in 1:npts]
        sa_c = mksmld(ax, ay, cam)
        sb_c = mksmld(bx, by, cam)
        (r_c, g_c, _) = SMLMAnalysis._compute_crosscorr(sa_c, sb_c, CrossCorrConfig())
        mid = 15:35        # r ≈ 0.145–0.345 μm
        @test isapprox(sum(g_c[mid]) / length(mid), 1.0; atol=0.1)
        @test all(0.7 .< g_c[mid] .< 1.3)
    end

    @testset "Bleaching fit degeneracy guard" begin
        Random.seed!(42)

        # Flat data: constant ~30 locs/frame with Poisson-like noise, 19k frames.
        # Reproduces the @hstirf case that produced a=-42873, k≈0 with unbounded NelderMead.
        # Expected: reject as degenerate, return nothing.
        flat = [max(0, round(Int, 30 + randn() * sqrt(30))) for _ in 1:19_000]
        @test SMLMAnalysis._estimate_bleaching_rate(flat) === nothing

        # Pure exponential decay: should recover parameters accurately (non-regression).
        t_exp = 1:5000
        exp_data = [max(0, round(Int, 10 + 50 * exp(-0.0005 * i) + randn() * 2)) for i in t_exp]
        res_exp = SMLMAnalysis._estimate_bleaching_rate(exp_data)
        @test res_exp !== nothing
        @test isapprox(res_exp.k_bleach, 5e-4, rtol=0.1)
        @test isapprox(res_exp.offset, 10.0, atol=2.0)
        @test isapprox(res_exp.N_0, 50.0, atol=3.0)
        @test res_exp.r_squared > 0.9

        # Bleach-then-flat (realistic DNA-PAINT): should still fit with physical params.
        t_bf = 1:10_000
        bf_data = [max(0, round(Int, 20 + 30 * exp(-0.001 * i) + randn() * 1.5)) for i in t_bf]
        res_bf = SMLMAnalysis._estimate_bleaching_rate(bf_data)
        @test res_bf !== nothing
        @test res_bf.offset >= 0     # physical bound held
        @test res_bf.N_0 >= 0        # physical bound held
        @test res_bf.k_bleach > 0
    end

    @testset "crop axis conventions" begin
        # crop_images(imgs, roi_x, roi_y) == imgs[roi_y, roi_x, :]: roi_x indexes
        # columns (x, dim 2), roi_y indexes rows (y, dim 1). Encode (row, col) into
        # each pixel so an axis swap is caught. Locks a convention documented but
        # otherwise untested (SMART transposes on load, MIC does not, overlays
        # transpose before drawing — all easy to get backwards).
        nrow, ncol, nfr = 6, 8, 3
        img = [1000r + 10c + f for r in 1:nrow, c in 1:ncol, f in 1:nfr]
        roi_x, roi_y = 3:6, 2:4        # columns, rows
        cropped = crop_images(img, roi_x, roi_y)
        @test size(cropped) == (length(roi_y), length(roi_x), nfr)
        @test cropped == img[roi_y, roi_x, :]
        @test cropped[1, 1, 1] == 1000 * first(roi_y) + 10 * first(roi_x) + 1

        # crop_camera uses the same convention: roi_x → x-edges, roi_y → y-edges.
        cam = IdealCamera(ncol, nrow, 0.1)   # IdealCamera(nx=cols, ny=rows, px)
        cc = crop_camera(cam, roi_x, roi_y)
        @test cc.pixel_edges_x == cam.pixel_edges_x[first(roi_x):last(roi_x)+1]
        @test cc.pixel_edges_y == cam.pixel_edges_y[first(roi_y):last(roi_y)+1]
        @test length(cc.pixel_edges_x) - 1 == length(roi_x)   # x pixel count = #cols
        @test length(cc.pixel_edges_y) - 1 == length(roi_y)   # y pixel count = #rows
    end

    @testset "SMLD HDF5 round-trip" begin
        # Locks the σ_xy regression: save_smld/load_smld must preserve every
        # emitter field, including the position covariance σ_xy that the
        # GaussianXYNBS → Emitter2DFitSigma path carries. This bug survived
        # because the prior tests only constructed types, never round-tripped.
        cam = IdealCamera(8, 8, 0.1)
        T = Float64

        mktempdir() do dir
            # Emitter2DFitSigma (16 fields) — the primary GaussianXYNBS output.
            es = [GaussMLE.Emitter2DFitSigma{T}(
                    0.1i, 0.2i, 1000.0 + i, 5.0, 0.13,     # x, y, photons, bg, σ
                    0.01, 0.012, 0.003, 20.0, 0.5, 0.002,  # σ_x, σ_y, σ_xy, σ_photons, σ_bg, σ_σ
                    0.4, i, 1, 0, i)                        # pvalue, frame, dataset, track_id, id
                  for i in 1:5]
            smld = BasicSMLD(es, cam, 10, 1, Dict{String,Any}())
            path = joinpath(dir, "sigma.h5")
            save_smld(path, smld)
            loaded = load_smld(path)

            @test length(loaded.emitters) == 5
            @test eltype(loaded.emitters) <: GaussMLE.Emitter2DFitSigma
            for (a, b) in zip(smld.emitters, loaded.emitters)
                for f in fieldnames(GaussMLE.Emitter2DFitSigma)
                    @test getfield(a, f) ≈ getfield(b, f)
                end
            end
            @test loaded.emitters[3].σ_xy ≈ 0.003   # the field that used to vanish

            # Emitter2DFitSigmaXY (18 fields) — GaussianXYNBSXSY output.
            exy = [GaussMLE.Emitter2DFitSigmaXY{T}(
                    0.1i, 0.2i, 1000.0 + i, 5.0, 0.13, 0.14, # x, y, photons, bg, σx, σy
                    0.01, 0.012, 0.003, 20.0, 0.5,           # σ_x, σ_y, σ_xy, σ_photons, σ_bg
                    0.002, 0.0021, 0.4, i, 1, 0, i)          # σ_σx, σ_σy, pvalue, frame, dataset, track_id, id
                  for i in 1:4]
            smld_xy = BasicSMLD(exy, cam, 10, 1, Dict{String,Any}())
            pxy = joinpath(dir, "sigmaxy.h5")
            save_smld(pxy, smld_xy)
            loaded_xy = load_smld(pxy)

            @test length(loaded_xy.emitters) == 4
            @test eltype(loaded_xy.emitters) <: GaussMLE.Emitter2DFitSigmaXY
            for (a, b) in zip(smld_xy.emitters, loaded_xy.emitters)
                for f in fieldnames(GaussMLE.Emitter2DFitSigmaXY)
                    @test getfield(a, f) ≈ getfield(b, f)
                end
            end
            @test loaded_xy.emitters[2].σ_xy ≈ 0.003

            # Schema validation: a valid HDF5 file that isn't an SMLD fails with a
            # friendly ErrorException, not a raw KeyError deep in the read.
            bogus = joinpath(dir, "bogus.h5")
            SMLMAnalysis.HDF5.h5open(bogus, "w") do f
                f["junk"] = 1
            end
            @test_throws ErrorException load_smld(bogus)
        end
    end

    @testset "GaussMLE emitter types + abstract-eltype round-trip" begin
        # Regression coverage for the AbstractEmitter type-erasure + serialization fix:
        #  - _with_dataset must handle EVERY advertised GaussMLE emitter type
        #    (GaussianXYNB → Emitter2DFitGaussMLE, AstigmaticXYZNB → Emitter3DFitGaussMLE),
        #    or detectfit throws a MethodError the moment those fitters are used.
        #  - save_smld/load_smld must round-trip those types and the standard-3D
        #    off-diagonal covariances, AND must key off the concrete emitter even when
        #    the SMLD is typed BasicSMLD{T,AbstractEmitter} (as the pre-narrowing
        #    pipeline produced) rather than degrading to Emitter2DFit and dropping σ.
        cam = IdealCamera(8, 8, 0.1)
        T = Float64

        # _with_dataset must not MethodError for the fixed-width / astigmatic types.
        e2g = GaussMLE.Emitter2DFitGaussMLE{T}(0.1, 0.2, 1000.0, 5.0,
                0.01, 0.012, 0.003, 20.0, 0.5, 0.4, 1, 1, 0, 1)
        e3g = GaussMLE.Emitter3DFitGaussMLE{T}(0.1, 0.2, 0.3, 1000.0, 5.0,
                0.01, 0.012, 0.02, 0.003, 0.001, 0.002, 20.0, 0.5, 0.4, 1, 1, 0, 1)
        @test SMLMAnalysis._with_dataset(e2g, 7).dataset == 7
        @test SMLMAnalysis._with_dataset(e2g, 7).σ_xy ≈ 0.003
        @test SMLMAnalysis._with_dataset(e3g, 7).dataset == 7
        @test SMLMAnalysis._with_dataset(e3g, 7).σ_yz ≈ 0.002

        mktempdir() do dir
            # Emitter2DFitGaussMLE (GaussianXYNB) round-trip.
            g2 = [GaussMLE.Emitter2DFitGaussMLE{T}(0.1i, 0.2i, 1000.0 + i, 5.0,
                    0.01, 0.012, 0.003, 20.0, 0.5, 0.4, i, 1, 0, i) for i in 1:4]
            s2 = BasicSMLD(g2, cam, 10, 1, Dict{String,Any}())
            p2 = joinpath(dir, "g2.h5"); save_smld(p2, s2); l2 = load_smld(p2)
            @test eltype(l2.emitters) <: GaussMLE.Emitter2DFitGaussMLE
            for (a, b) in zip(s2.emitters, l2.emitters), f in fieldnames(GaussMLE.Emitter2DFitGaussMLE)
                @test getfield(a, f) ≈ getfield(b, f)
            end

            # Emitter3DFitGaussMLE (AstigmaticXYZNB) round-trip: z + full covariance.
            g3 = [GaussMLE.Emitter3DFitGaussMLE{T}(0.1i, 0.2i, 0.3i, 1000.0 + i, 5.0,
                    0.01, 0.012, 0.02, 0.003, 0.001, 0.002, 20.0, 0.5, 0.4, i, 1, 0, i) for i in 1:4]
            s3 = BasicSMLD(g3, cam, 10, 1, Dict{String,Any}())
            p3 = joinpath(dir, "g3.h5"); save_smld(p3, s3); l3 = load_smld(p3)
            @test eltype(l3.emitters) <: GaussMLE.Emitter3DFitGaussMLE
            @test l3.emitters[2].z ≈ 0.6
            @test l3.emitters[2].σ_xz ≈ 0.001
            @test l3.emitters[2].σ_yz ≈ 0.002

            # Standard Emitter3DFit: off-diagonal covariances σ_xz/σ_yz survive
            # (they were never written before this fix).
            e3 = [Emitter3DFit{T}(0.1i, 0.2i, 0.3i, 1000.0 + i, 5.0,
                    0.01, 0.012, 0.02, 20.0, 0.5;
                    σ_xy=0.003, σ_xz=0.001, σ_yz=0.002, frame=i, dataset=1, id=i) for i in 1:4]
            s3s = BasicSMLD(e3, cam, 10, 1, Dict{String,Any}())
            p3s = joinpath(dir, "e3.h5"); save_smld(p3s, s3s); l3s = load_smld(p3s)
            @test eltype(l3s.emitters) <: Emitter3DFit
            @test l3s.emitters[2].σ_xz ≈ 0.001
            @test l3s.emitters[2].σ_yz ≈ 0.002

            # Abstract-eltype SMLD (as the pipeline produced before narrowing):
            # save_smld must reload it as concrete Emitter2DFitSigma, NOT degrade to
            # Emitter2DFit and drop the PSF-width σ.
            abs_v = AbstractEmitter[GaussMLE.Emitter2DFitSigma{T}(
                        0.1i, 0.2i, 1000.0 + i, 5.0, 0.13,
                        0.01, 0.012, 0.003, 20.0, 0.5, 0.002,
                        0.4, i, 1, 0, i) for i in 1:4]
            @test eltype(abs_v) == AbstractEmitter
            s_abs = BasicSMLD(abs_v, cam, 10, 1, Dict{String,Any}())
            pabs = joinpath(dir, "abs.h5"); save_smld(pabs, s_abs); labs = load_smld(pabs)
            @test eltype(labs.emitters) <: GaussMLE.Emitter2DFitSigma   # NOT Emitter2DFit
            @test labs.emitters[2].σ ≈ 0.13                             # PSF-width σ preserved
            @test labs.emitters[2].σ_xy ≈ 0.003
        end
    end
end

if SMLM_TEST_FULL
    @testset "thorough" begin
        @testset "end-to-end analyze() tuple contract" begin
            # Runs the full detect/fit → filter → render pipeline on synthetic
            # CPU-only data and checks the (result, info) contract the whole
            # package rests on. Mirrors the precompile workload (known-good),
            # but as an assertable test rather than a build-time smoke run.
            Random.seed!(1)
            cam = IdealCamera(32, 32, 0.1)
            sim = StaticSMLMConfig(density = 5.0, σ_psf = 0.13, nframes = 50, ndatasets = 1)
            (_, si) = simulate(sim;
                pattern  = Nmer2D(n = 8, d = 0.05),
                molecule = GenericFluor(photons = 5.0e4, k_off = 20.0, k_on = 0.04),
                camera   = cam)
            (imgs, _) = gen_images(si.smld_model, SMLMAnalysis.MicroscopePSFs.GaussianPSF(0.13);
                dataset = 1, bg = 20.0, poisson_noise = true)

            cfg = AnalysisConfig(
                DetectFitConfig(boxer  = BoxerConfig(boxsize = 7, psf_sigma = 0.13),
                                fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu)),
                FilterConfig(photons = (100.0, Inf)),
                RenderConfig(zoom = 10);
                camera  = cam,
                verbose = Verbosity.SILENT,
            )
            (result, info) = analyze([imgs], cfg)

            @test result isa AnalysisResult
            @test info isa AnalysisInfo
            @test result.smld isa BasicSMLD
            @test length(result.smld.emitters) > 0
            @test length(info.step_infos) == 3
            @test info.elapsed_s >= 0

            # detectfit narrows the accumulated AbstractEmitter[] to a concrete type;
            # the default GaussianXYNBS fitter yields Emitter2DFitSigma.
            @test isconcretetype(eltype(result.smld.emitters))
            @test eltype(result.smld.emitters) <: GaussMLE.Emitter2DFitSigma

            # The fitted output round-trips through HDF5 losslessly — concrete type and
            # PSF-width σ preserved (σ was silently lost when the SMLD was AbstractEmitter-typed).
            mktempdir() do dir
                p = joinpath(dir, "pipeline.h5")
                save_smld(p, result.smld)
                reloaded = load_smld(p)
                @test length(reloaded.emitters) == length(result.smld.emitters)
                @test eltype(reloaded.emitters) <: GaussMLE.Emitter2DFitSigma
                @test reloaded.emitters[1].σ ≈ result.smld.emitters[1].σ
                @test reloaded.emitters[1].σ_xy ≈ result.smld.emitters[1].σ_xy
            end

            # Multi-dataset detectfit: exercises the per-dataset loop (dataset field,
            # per-dataset frame numbering, n_datasets tracking) that the in-memory and
            # file-based paths share. Reuse the same stack as two datasets.
            (smld2, si2) = analyze([imgs, imgs],
                DetectFitConfig(camera = cam,
                                boxer  = BoxerConfig(boxsize = 7, psf_sigma = 0.13),
                                fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu));
                verbose = Verbosity.SILENT)
            @test si2.info.n_datasets == 2
            @test smld2.n_datasets == 2
            @test smld2.n_frames == size(imgs, 3)          # equal-length → per-dataset count
            @test Set(e.dataset for e in smld2.emitters) == Set([1, 2])
            @test all(1 <= e.frame <= size(imgs, 3) for e in smld2.emitters)  # frames are per-dataset
        end
    end
else
    @info "Skipping thorough tests; set SMLM_TEST_FULL=1 to enable"
end
