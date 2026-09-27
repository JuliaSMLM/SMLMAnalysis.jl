using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Aqua
using Random
using TOML
using Statistics

const SMLM_TEST_FULL = lowercase(get(ENV, "SMLM_TEST_FULL", "false")) in ("true", "1", "yes")

@testset "fast" begin
    @testset "exports resolve" begin
        # Every exported name must have a defined binding. This used to catch a
        # dangling `export AbstractCamera`: SMLMData and Makie (via CairoMakie)
        # both export an AbstractCamera, and the ambiguity left the module
        # binding undeclared even though `export` listed it. AbstractCamera is
        # no longer exported (it's an extension hook, reached as
        # `SMLMAnalysis.SMLMData.AbstractCamera`), but this sweep still catches
        # any future collision a new dependency introduces.
        dangling = [n for n in names(SMLMAnalysis) if !isdefined(SMLMAnalysis, n)]
        @test isempty(dangling)
        @test IdealCamera <: SMLMAnalysis.SMLMData.AbstractCamera
    end

    @testset "exact export set" begin
        # Every export must justify itself (Keith's ruling, 2026-09):
        # build-a-pipeline / run-it / pick-a-result-out / load-or-save / manage-the-guide.
        # This is the literal list, so any future export is a deliberate edit here.
        keep = [
            :analyze, :stepinfo, :stepinfos, :Verbosity, :Checkpoint,
            :AnalysisConfig, :MultiTargetConfig,
            :DetectFitConfig, :FilterConfig, :DensityFilterConfig, :IntensityFilterConfig,
            :CompositeRenderConfig, :CrossAlignConfig, :CrossCorrConfig,
            :FrameConnectConfig, :CalibrationConfig, :DriftConfig, :AlignConfig,
            :RenderConfig, :BaGoLConfig, :BoxerConfig, :GaussMLEConfig,
            :GaussianXYNB, :GaussianXYNBS, :GaussianXYNBSXSY, :AstigmaticXYZNB,
            :HistogramRender, :GaussianRender, :CircleRender, :EllipseRender,
            :DBSCANConfig, :HDBSCANConfig, :HierarchicalConfig, :VoronoiConfig,
            :HopkinsConfig, :VoronoiDensityConfig, :OuterPolygonConfig, :KdeValleyConfig,
            :IdealCamera, :SCMOSCamera,
            :save_smld, :load_smld, :smld_info, :load_smart_h5, :load_smart_h5_info,
            :smart_h5_to_array, :load_mic_h5, :load_mic_h5_info, :load_mic_h5_block,
            :build_camera_from_mic_h5,
            :install_agent_guide, :uninstall_agent_guide, :agent_guide_status,
        ]
        @test Set(names(SMLMAnalysis)) == Set([:SMLMAnalysis; keep])
    end

    @testset "non-exported but public names resolve" begin
        # api_overview.md's "Non-exported but public" section promises every one of
        # these resolves as `SMLMAnalysis.Name`, exported or not. Hard-coded here (not
        # parsed from the doc) so a naming collision with a new `using`d dependency —
        # like the CairoMakie/Makie `AbstractCamera` one this guards against — is
        # caught immediately instead of silently leaving a dangling binding.
        public_names = [
            # Owned result/info structs and extension hooks.
            :AnalysisResult, :AnalysisInfo, :StepInfo, :AbstractSMLMConfig, :AbstractSMLMInfo,
            :AbstractMultiTargetStep, :MultiTargetResult, :MultiTargetInfo, :DetectFitInfo,
            :FilterInfo, :DensityFilterInfo, :IntensityFilterInfo, :BaGoLInfo,
            :CompositeRenderInfo, :CrossAlignInfo, :CrossCorrInfo, :step_name, :step_outdir,
            # Upstream re-exports, qualified with their owning package in the docs but
            # also reachable one level down as SMLMAnalysis.Name.
            :AbstractCamera, :BasicSMLD, :Emitter2DFit, :Emitter3DFit, :ROIBatch,
            :StaticSMLMConfig, :simulate, :gen_images, :Nmer2D, :Line2D, :GenericFluor,
            :fit, :frameconnect, :CalibrationResult, :driftcorrect, :align_smld, :AlignInfo,
            :run_bagol, :BaGoLDiagnostics, :render, :cluster, :cluster_statistics,
            :AbstractClusterConfig, :AbstractStatisticsConfig, :ClusterInfo, :ClusterStatisticsInfo,
            :in_cell, :interior_mask, :interior_fraction, :AbstractEdgeClassifyConfig,
            :EdgeClassifyInfo, :CellPolygon, :MultiCellMask,
        ]
        for n in public_names
            @test (try; getglobal(SMLMAnalysis, n); true; catch; false; end)
        end
        @test SMLMAnalysis.AbstractCamera === SMLMAnalysis.SMLMData.AbstractCamera
    end

    @testset "Aqua" begin
        # `ambiguities=(recursive=false,)`: recursive ambiguity checking also flags
        # method ambiguities defined entirely inside our upstream dependencies
        # (SMLMData/SMLMRender/etc.), which are not ours to fix here.
        # Skipped in the downgrade-compat CI job (SMLM_DOWNGRADE_CI=true): there the
        # oldest-allowed upstream versions carry their own ambiguities, the action
        # merges test extras into [deps] (so Aqua flags itself as stale), and the
        # persistent-task probe cannot precompile. Aqua runs on every other job.
        if get(ENV, "SMLM_DOWNGRADE_CI", "false") != "true"
            Aqua.test_all(SMLMAnalysis; ambiguities=(recursive=false,))
        end
    end

    @testset "docs cover every SMLMAnalysis-owned export" begin
        # Every name DEFINED in SMLMAnalysis (not merely re-exported from an upstream
        # package) must appear in a ```@docs block under docs/src, or the hosted manual
        # silently omits it. `checkdocs=:none` in docs/make.jl cannot catch this because
        # the upstream modules are in `modules=`; this is the scoped check.
        docsdir = joinpath(dirname(@__DIR__), "docs", "src")
        documented = Set{String}()
        for (root, _, files) in walkdir(docsdir), f in files
            endswith(f, ".md") || continue
            inblock = false
            for ln in eachline(joinpath(root, f))
                s = strip(ln)
                if startswith(s, "```@docs")
                    inblock = true
                elseif inblock && startswith(s, "```")
                    inblock = false
                elseif inblock && !isempty(s)
                    push!(documented, s)
                end
            end
        end
        owned = String[]
        for n in names(SMLMAnalysis)
            n === :SMLMAnalysis && continue
            obj = getfield(SMLMAnalysis, n)
            obj isa Module && continue            # Verbosity / Checkpoint: documented as tables
            parentmodule(obj) === SMLMAnalysis || continue
            push!(owned, string(n))
        end
        @test !isempty(owned)
        undocumented = sort!(setdiff(owned, documented))
        @test isempty(undocumented)
    end

    @testset "Types" begin
        # Test SMLMAnalysis.AnalysisInfo constructor
        info = SMLMAnalysis.AnalysisInfo()
        @test info.elapsed_s == 0.0
        @test isempty(info.step_infos)

        # Test SMLMAnalysis.AnalysisInfo with data
        cfg0 = FilterConfig()
        si0 = SMLMAnalysis.StepInfo(1, cfg0, 0.2, Dict{Symbol,Any}(); info=SMLMAnalysis.FilterInfo(10, 8, 0.2))
        info = SMLMAnalysis.AnalysisInfo(1.5, SMLMAnalysis.StepInfo[si0])
        @test info.elapsed_s == 1.5
        @test length(info.step_infos) == 1
        @test info.step_infos[1].info isa SMLMAnalysis.FilterInfo

        # Test SMLMAnalysis.StepInfo with typed info
        cfg = FilterConfig()
        filter_info = SMLMAnalysis.FilterInfo(100, 80, 0.5)
        step_info = SMLMAnalysis.StepInfo(1, cfg, 0.5, Dict{Symbol,Any}(:n_before => 100); info=filter_info)
        @test step_info.info !== nothing
        @test step_info.info isa SMLMAnalysis.FilterInfo
        @test step_info.info.n_before == 100
        @test step_info.info.n_after == 80
        @test step_info.elapsed_s == 0.5

        # Test SMLMAnalysis.StepInfo without info
        step_info2 = SMLMAnalysis.StepInfo(2, cfg, 0.3, Dict{Symbol,Any}())
        @test step_info2.info === nothing

        # Test native info structs
        # Back-compat 7-arg constructor (defaults selected_source_indices to nothing)
        di = SMLMAnalysis.DetectFitInfo([], [], 2, 1000, 950, 5000, 1.5)
        @test di.n_datasets == 2
        @test di.n_rois == 1000
        @test di.n_fits == 950
        @test di.selected_source_indices === nothing

        @test di.n_nonfinite == 0

        # Full 8-arg constructor with provenance
        di_sel = SMLMAnalysis.DetectFitInfo([], [], 3, 500, 450, 1000, 0.5, [1, 3, 5])
        @test di_sel.selected_source_indices == [1, 3, 5]
        @test di_sel.n_datasets == 3
        @test di_sel.n_nonfinite == 0

        # _drop_nonfinite_emitters: drops emitters with a NaN/Inf in any AbstractFloat field.
        e_ok = SMLMAnalysis.Emitter2DFit(1.0, 1.0, 100.0, 5.0, 0.01, 0.01, 0.0, 1.0, 1.0, 1, 1, 0, 1)
        e_nan_y = SMLMAnalysis.Emitter2DFit(1.0, NaN, 100.0, 5.0, 0.01, 0.01, 0.0, 1.0, 1.0, 1, 1, 0, 2)
        e_inf_photons = SMLMAnalysis.Emitter2DFit(1.0, 1.0, Inf, 5.0, 0.01, 0.01, 0.0, 1.0, 1.0, 1, 1, 0, 3)
        kept, n_dropped = SMLMAnalysis._drop_nonfinite_emitters([e_ok, e_nan_y, e_inf_photons])
        @test kept == [e_ok]
        @test n_dropped == 2

        dfi = SMLMAnalysis.DensityFilterInfo(1000, 800, 5, 0.3)
        @test dfi.n_before == 1000
        @test dfi.threshold == 5

        # densityfilter's per-point query radius n_sigma*sqrt(σ_i^2+max_σ^2) is an
        # exact upper bound for the pair test dist < n_sigma*sqrt(σ_i^2+σ_j^2) it
        # gates on, so switching the KD-tree query from a single global radius to a
        # per-point one must not change which emitters survive. Cross-check against
        # an independent brute-force count of that same pair test.
        cam = IdealCamera(64, 64, 0.1)
        xs = [0.0, 0.05, 0.0, 0.05, 5.0]
        ys = [0.0, 0.0, 0.05, 0.05, 5.0]
        σs = [0.01, 0.01, 0.01, 0.05, 0.2]
        emitters = [SMLMAnalysis.Emitter2DFit(xs[i], ys[i], 1000.0, 10.0, σs[i], σs[i], 0.0, 1.0, 1.0, i, 1, 0, i) for i in 1:5]
        smld_df = SMLMAnalysis.BasicSMLD(emitters, cam, 1, 1, Dict{String,Any}())
        cfg_df = DensityFilterConfig(n_sigma=3.0, min_neighbors=1)
        filtered_df, _ = SMLMAnalysis.densityfilter_step(smld_df, cfg_df)

        σ = [sqrt(e.σ_x^2 + e.σ_y^2) for e in emitters]
        expected_counts = zeros(Int, 5)
        for i in 1:5, j in 1:5
            i == j && continue
            dist = sqrt((emitters[i].x - emitters[j].x)^2 + (emitters[i].y - emitters[j].y)^2)
            dist < cfg_df.n_sigma * sqrt(σ[i]^2 + σ[j]^2) && (expected_counts[i] += 1)
        end
        @test length(filtered_df.emitters) == count(>=(cfg_df.min_neighbors), expected_counts)

        # Non-finite σ is rejected up front rather than silently propagating.
        bad_emitters = [SMLMAnalysis.Emitter2DFit(0.0, 0.0, 1000.0, 10.0, NaN, 0.01, 0.0, 1.0, 1.0, 1, 1, 0, 1)]
        smld_bad = SMLMAnalysis.BasicSMLD(bad_emitters, cam, 1, 1, Dict{String,Any}())
        @test_throws ArgumentError SMLMAnalysis.densityfilter_step(smld_bad, cfg_df)

        # IntensityFilter: non-finite photons rejected up front. Needs >=100
        # emitters, otherwise the step's own "too few emitters" early return
        # would skip the check before it is reached.
        many = [SMLMAnalysis.Emitter2DFit(0.01i, 0.01i, 1000.0, 10.0, 0.01, 0.01, 0.0, 1.0, 1.0, i, 1, 0, i) for i in 2:100]
        bad_photon = SMLMAnalysis.Emitter2DFit(0.5, 0.5, NaN, 10.0, 0.01, 0.01, 0.0, 1.0, 1.0, 1, 1, 0, 1)
        smld_ifbad = SMLMAnalysis.BasicSMLD(vcat([bad_photon], many), cam, 1, 1, Dict{String,Any}())
        @test_throws ArgumentError SMLMAnalysis.intensityfilter_step(smld_ifbad, IntensityFilterConfig())

        # psf_sigma: an explicit (lo, hi) — including a 0.0 lower bound — is
        # always applied, unlike :auto's "skip when degenerate" behavior.
        mkem_sigma(σ, i) = GaussMLE.Emitter2DFitSigma{Float64}(
            0.0, 0.0, 1000.0, 5.0, σ, 0.01, 0.01, 0.0, 20.0, 0.5, 0.002, 1.0, 1, 1, 0, i)
        smld_ps = SMLMAnalysis.BasicSMLD([mkem_sigma(0.05, 1), mkem_sigma(0.3, 2)], cam, 1, 1, Dict{String,Any}())
        filtered_ps, _ = SMLMAnalysis.filter_step(smld_ps, FilterConfig(psf_sigma=(0.0, 0.1)))
        @test length(filtered_ps.emitters) == 1
        @test filtered_ps.emitters[1].σ == 0.05
    end

    @testset "analyze dispatch" begin
        # Verify analyze() dispatch methods exist for each step config type
        @test hasmethod(analyze, Tuple{Vector{<:AbstractArray{<:Real,3}}, DetectFitConfig})
        @test hasmethod(analyze, Tuple{AbstractArray{<:Real,3}, DetectFitConfig})
        @test hasmethod(analyze, Tuple{DetectFitConfig})
        @test hasmethod(analyze, Tuple{SMLMAnalysis.BasicSMLD, FilterConfig})
        @test hasmethod(analyze, Tuple{SMLMAnalysis.BasicSMLD, FrameConnectConfig})
        @test hasmethod(analyze, Tuple{SMLMAnalysis.BasicSMLD, DriftConfig})
        @test hasmethod(analyze, Tuple{SMLMAnalysis.BasicSMLD, DensityFilterConfig})
        @test hasmethod(analyze, Tuple{SMLMAnalysis.BasicSMLD, RenderConfig})

        # analyze(data, config::AnalysisConfig) with a data type _normalize_data
        # doesn't recognize (e.g. a String — data is never itself a file path;
        # use `nothing` with a file-based DetectFitConfig instead).
        @test_throws ArgumentError analyze("some/path.h5", AnalysisConfig(camera=IdealCamera(8, 8, 0.1)))

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

        # _inject_camera (AnalysisConfig pipeline path): pixel_size/qe on the
        # DetectFitConfig would be silently ignored (the pipeline camera always
        # wins), so it's rejected instead of accepted-and-dropped.
        @test SMLMAnalysis._inject_camera(cfg, cam).camera === cam    # no pixel_size/qe: fine
        @test_throws ArgumentError SMLMAnalysis._inject_camera(DetectFitConfig(pixel_size=0.1), cam)
        @test_throws ArgumentError SMLMAnalysis._inject_camera(DetectFitConfig(qe=0.9), cam)
        @test SMLMAnalysis._inject_camera(DetectFitConfig(camera=cam, pixel_size=0.1), cam).camera === cam  # camera already set: unchanged

        # DetectFitConfig.datasets selection field
        @test cfg.datasets === nothing                          # default is no selection
        cfg_range = DetectFitConfig(datasets=1:19)
        @test cfg_range.datasets == 1:19
        @test cfg_range.datasets isa Vector{Int}   # concrete field type: any AbstractVector{Int} is accepted but stored as Vector{Int}
        cfg_sparse = DetectFitConfig(datasets=[1, 2, 3, 5, 7])
        @test cfg_sparse.datasets == [1, 2, 3, 5, 7]
        @test cfg_sparse.datasets isa Vector{Int}

        # _select_sources: pass-through when nothing, bounds-checked otherwise
        src = [(i=j,) for j in 1:5]
        @test SMLMAnalysis._select_sources(src, nothing) === src
        @test SMLMAnalysis._select_sources(src, [1, 3, 5]) == [src[1], src[3], src[5]]
        @test SMLMAnalysis._select_sources(src, 2:4) == src[2:4]
        @test_throws ArgumentError SMLMAnalysis._select_sources(src, [1, 6])
        @test_throws ArgumentError SMLMAnalysis._select_sources(src, [0, 1])

        # No kwargs... catch-all on step analyze() methods: a misspelled keyword
        # must raise MethodError, not silently vanish into a kwargs sink.
        smld_empty = SMLMAnalysis.BasicSMLD(SMLMAnalysis.Emitter2DFit{Float64}[], cam, 1, 1, Dict{String,Any}())
        @test_throws MethodError analyze(smld_empty, FilterConfig(); bogus_kwarg=1)
    end

    @testset "Verbosity/Checkpoint validation" begin
        cam = IdealCamera(64, 64, 0.1)

        # Valid levels round-trip through AnalysisConfig / MultiTargetConfig construction.
        @test AnalysisConfig(camera=cam, verbose=Verbosity.SILENT, checkpoint=Checkpoint.NONE).verbose == Verbosity.SILENT
        @test AnalysisConfig(camera=cam, verbose=Verbosity.DEBUG, checkpoint=Checkpoint.ALL).checkpoint == Checkpoint.ALL
        @test MultiTargetConfig(labels=[:A], outdir="x", verbose=Verbosity.DEBUG).verbose == Verbosity.DEBUG

        # Out-of-range verbose/checkpoint must raise ArgumentError at construction.
        @test_throws ArgumentError AnalysisConfig(camera=cam, verbose=-1)
        @test_throws ArgumentError AnalysisConfig(camera=cam, verbose=Verbosity.DEBUG + 1)
        @test_throws ArgumentError AnalysisConfig(camera=cam, checkpoint=-1)
        @test_throws ArgumentError AnalysisConfig(camera=cam, checkpoint=Checkpoint.ALL + 1)
        @test_throws ArgumentError MultiTargetConfig(labels=[:A], outdir="x", verbose=-1)
        @test_throws ArgumentError MultiTargetConfig(labels=[:A], outdir="x", verbose=Verbosity.DEBUG + 1)

        # Same validation at the _run_pipeline entry point, for direct calls that
        # bypass AnalysisConfig entirely.
        steps = SMLMAnalysis.AbstractSMLMConfig[]
        @test_throws ArgumentError SMLMAnalysis._run_pipeline(nothing, steps, cam, nothing, -1)
        @test_throws ArgumentError SMLMAnalysis._run_pipeline(nothing, steps, cam, nothing, Verbosity.STANDARD, -1)
        @test_throws ArgumentError SMLMAnalysis._run_pipeline(nothing, steps, cam, nothing, Verbosity.STANDARD, Checkpoint.ALL + 1)
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
        # All info structs should be SMLMAnalysis.AbstractSMLMInfo subtypes
        @test SMLMAnalysis.StepInfo <: SMLMAnalysis.AbstractSMLMInfo
        @test SMLMAnalysis.DetectFitInfo <: SMLMAnalysis.AbstractSMLMInfo
        @test SMLMAnalysis.FilterInfo <: SMLMAnalysis.AbstractSMLMInfo
        @test SMLMAnalysis.DensityFilterInfo <: SMLMAnalysis.AbstractSMLMInfo
        @test SMLMAnalysis.CompositeRenderInfo <: SMLMAnalysis.AbstractSMLMInfo
        @test SMLMAnalysis.CrossAlignInfo <: SMLMAnalysis.AbstractSMLMInfo
        @test SMLMAnalysis.AnalysisInfo <: SMLMAnalysis.AbstractSMLMInfo
    end

    @testset "Multi-target step types" begin
        # Type hierarchy
        @test SMLMAnalysis.AbstractMultiTargetStep <: SMLMAnalysis.AbstractSMLMConfig
        @test CompositeRenderConfig <: SMLMAnalysis.AbstractMultiTargetStep
        @test CrossAlignConfig <: SMLMAnalysis.AbstractMultiTargetStep

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
        @test ca.align isa AlignConfig
        @test ca.align.method == :entropy
        @test ca.align.maxn == 100
        @test ca.align.histbinsize == 0.05

        # CrossAlignConfig custom: upstream AlignConfig passed through as-is
        ca2 = CrossAlignConfig(align=AlignConfig(method=:fft, maxn=50, verbose=1))
        @test ca2.align.method == :fft
        @test ca2.align.maxn == 50
        @test ca2.align.verbose == 1

        # step_name dispatch
        @test SMLMAnalysis.step_name(cr) == "compositerender"
        @test SMLMAnalysis.step_name(ca) == "crossalign"

        # analyze dispatch methods exist for multi-target steps
        @test hasmethod(analyze, Tuple{Vector{<:SMLMAnalysis.BasicSMLD}, CompositeRenderConfig})
        @test hasmethod(analyze, Tuple{Vector{<:SMLMAnalysis.BasicSMLD}, CrossAlignConfig})

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

    @testset "multi-target config.toml is valid with two nested-config steps" begin
        # CompositeRenderConfig (strategy) and CrossAlignConfig (align) each carry a
        # nested config field. Written into the same multi_target_config.toml under
        # separate [[steps]] entries, a bare `[strategy]` / `[align]` table floats to
        # the document root and collides across steps -- TOML.parsefile then fails
        # with "key already defined". table_prefix="steps." scopes each nested table
        # to its own [[steps]] array element instead.
        mktempdir() do dir
            mt = MultiTargetConfig(
                labels=[:A, :B],
                steps=[
                    CompositeRenderConfig(),
                    CrossAlignConfig(),
                    CompositeRenderConfig(),
                ],
                outdir=dir,
            )
            SMLMAnalysis._save_multitarget_config!(mt)
            parsed = TOML.parsefile(joinpath(dir, "multi_target_config.toml"))
            @test length(parsed["steps"]) == 3
            @test haskey(parsed["steps"][1], "strategy")
            @test haskey(parsed["steps"][2], "align")
            @test haskey(parsed["steps"][3], "strategy")
        end
    end

    @testset "multi-target saves and returns aligned channels" begin
        # Phase-2 dispatch + result assembly on hand-made SMLDs (no detectfit): channel
        # B is channel A shifted by a known offset; after CrossAlign the saved file,
        # result[:B].smld and result.smlds[2] must all be the aligned data.
        rng = MersenneTwister(3)
        cam = IdealCamera(64, 64, 0.1)
        N = 1500
        xs = 1.0 .+ 4.4 .* rand(rng, N); ys = 1.0 .+ 4.4 .* rand(rng, N)
        mk(dx, dy) = SMLMAnalysis.BasicSMLD([SMLMAnalysis.Emitter2DFit{Float64}(xs[i] + dx, ys[i] + dy, 1000.0, 10.0,
                                    0.01, 0.01, 50.0, 2.0; frame=1 + (i % 10)) for i in 1:N],
                               cam, 10, 1, Dict{String,Any}())
        a, b = mk(0.0, 0.0), mk(0.08, -0.05)
        labels = [:A, :B]
        (state, si) = analyze([a, b], CrossAlignConfig();
            outdir=nothing, step_number=1, verbose=0)
        meanxy(s) = (sum(e.x for e in s.emitters) / N, sum(e.y for e in s.emitters) / N)
        off(s1, s2) = hypot((meanxy(s2) .- meanxy(s1))...)
        @test off(state[1], state[2]) < 0.010          # known ~94 nm offset removed to < 10 nm

        channels = Dict{Symbol,SMLMAnalysis.AnalysisResult}(:A => SMLMAnalysis.AnalysisResult(a, a, nothing),
                                               :B => SMLMAnalysis.AnalysisResult(b, b, nothing))
        mktempdir() do dir
            SMLMAnalysis._finalize_channels!(channels, state, labels, dir; verbose=0)
            mtr = SMLMAnalysis.MultiTargetResult(labels, state, channels, [si], dir)
            @test mtr[:B].smld === mtr.smlds[2] === state[2]
            @test mtr[:B].smld_connected === b                # pre-alignment data kept
            saved = load_smld(joinpath(dir, "smld_B.h5"))
            @test [e.x for e in saved.emitters] == [e.x for e in state[2].emitters]
            @test off(saved, b) > 0.05                        # file holds aligned, not raw, B
        end

        # A state that is not one SMLD per label is a contract violation: throw
        # rather than silently leave the channel results alone (a caller that
        # continued on to _write_composite_readme! would otherwise BoundsError).
        ch2 = Dict{Symbol,SMLMAnalysis.AnalysisResult}(:A => SMLMAnalysis.AnalysisResult(a, nothing, nothing),
                                          :B => SMLMAnalysis.AnalysisResult(b, nothing, nothing))
        mktempdir() do dir
            @test_throws ArgumentError SMLMAnalysis._finalize_channels!(ch2, state[1:1], labels, dir; verbose=0)
            @test ch2[:B].smld === b
            @test isempty(readdir(dir))

            # Right length, but an untyped Vector{Any} container: rejected even
            # though its actual elements are BasicSMLDs -- _write_composite_readme!
            # requires Vector{<:SMLMAnalysis.BasicSMLD}, so this would otherwise MethodError
            # there instead of failing loudly here.
            @test_throws ArgumentError SMLMAnalysis._finalize_channels!(ch2, Any[state[1], state[2]], labels, dir; verbose=0)

            # A view is an AbstractVector{<:SMLMAnalysis.BasicSMLD} but not a Vector -- same
            # rejection, for the same reason (_write_composite_readme! requires
            # exactly Vector{<:SMLMAnalysis.BasicSMLD}).
            @test_throws ArgumentError SMLMAnalysis._finalize_channels!(ch2, view(state, 1:2), labels, dir; verbose=0)
        end
    end

    @testset "cross-align carries edge geometry with the emitters" begin
        # align_smld only moves emitters; edge-classification metadata
        # (edge_outer_polygon / edge_cells) must move with them or a saved
        # post-alignment SMLD pairs aligned emitters with a stale mask (Codex #51).
        CP = SMLMAnalysis.SMLMClustering.CellPolygon
        rng = MersenneTwister(11)
        cam = IdealCamera(64, 64, 0.1)
        N = 1500
        xs = 1.0 .+ 4.4 .* rand(rng, N); ys = 1.0 .+ 4.4 .* rand(rng, N)
        sq(x0, s) = NTuple{2,Float64}[(x0, x0), (x0 + s, x0), (x0 + s, x0 + s), (x0, x0 + s)]
        shiftpts(pts, dx, dy) = NTuple{2,Float64}[(p[1] + dx, p[2] + dy) for p in pts]

        outer_a = sq(1.0, 4.0)
        cells_a = [CP(sq(1.0, 4.0), [sq(1.5, 1.0)])]
        dx, dy = 0.094, 0.0
        outer_b = shiftpts(outer_a, dx, dy)
        cells_b = [CP(shiftpts(cells_a[1].outer, dx, dy),
                      [shiftpts(h, dx, dy) for h in cells_a[1].holes])]

        mkgeom(ddx, ddy, md) = SMLMAnalysis.BasicSMLD([SMLMAnalysis.Emitter2DFit{Float64}(xs[i] + ddx, ys[i] + ddy, 1000.0,
                                    10.0, 0.01, 0.01, 50.0, 2.0; frame=1 + (i % 10)) for i in 1:N],
                               cam, 10, 1, md)
        a = mkgeom(0.0, 0.0, Dict{String,Any}("edge_outer_polygon" => outer_a, "edge_cells" => cells_a))
        b = mkgeom(dx, dy, Dict{String,Any}("edge_outer_polygon" => outer_b, "edge_cells" => cells_b))

        # Default CrossAlignConfig() uses AlignConfig's default transform=:shift, so
        # every vertex and every emitter is offset by the exact same [dx, dy] —
        # emitter shift and polygon-vertex shift must agree exactly (not just close).
        (aligned, si) = analyze([a, b], CrossAlignConfig(); verbose=0)
        @test si.info isa SMLMAnalysis.CrossAlignInfo

        mean_x(s) = sum(e.x for e in s.emitters) / length(s.emitters)
        mean_y(s) = sum(e.y for e in s.emitters) / length(s.emitters)
        emitter_dx = mean_x(b) - mean_x(aligned[2])
        emitter_dy = mean_y(b) - mean_y(aligned[2])

        # Every vertex (outer ring, cell outer ring, and the cell's hole) must have
        # moved by the same [dx, dy] the emitters did, not merely a nearby amount.
        vertex_shift_matches(before, after) = all(
            isapprox(pb[1] - pa[1], emitter_dx; atol=1e-9) &&
            isapprox(pb[2] - pa[2], emitter_dy; atol=1e-9)
            for (pb, pa) in zip(before, after))

        @test vertex_shift_matches(outer_b, aligned[2].metadata["edge_outer_polygon"])
        aligned_cell = aligned[2].metadata["edge_cells"][1]
        @test vertex_shift_matches(cells_b[1].outer, aligned_cell.outer)
        @test vertex_shift_matches(cells_b[1].holes[1], aligned_cell.holes[1])

        # The reference channel (aligned[1] === a) is never touched.
        @test aligned[1] === a
        @test aligned[1].metadata["edge_outer_polygon"] == outer_a
        @test aligned[1].metadata["edge_cells"][1].outer == cells_a[1].outer
    end

    @testset "cross-align drops edge geometry under :affine alignment" begin
        # align_smld's :affine transform composes two sequential affine passes,
        # but info.diagnostic only records each pass's own coefficients summed
        # together, which is NOT the composed map (Codex measured a 50nm
        # mismatch at x=50μm with a=0.1, a2=0.01). Recovering the exact composed
        # map from the aligned emitter positions is unsound in general
        # (degenerate/collinear layouts give a spurious exact fit) and Float32
        # data can legitimately fail a tight residual tolerance, so the ruling is
        # to not carry geometry through :affine at all: both metadata keys are
        # dropped, with one @warn per channel.
        CP = SMLMAnalysis.SMLMClustering.CellPolygon
        cam = IdealCamera(16, 16, 0.1)
        em = [SMLMAnalysis.Emitter2DFit{Float64}(0.1i, 0.1i, 1000.0, 5.0, 0.01, 0.01, 20.0, 0.5; frame=i) for i in 1:3]
        outer = NTuple{2,Float64}[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)]
        hole = NTuple{2,Float64}[(1.5, 1.5), (2.0, 1.5), (2.0, 2.0)]
        md = Dict{String,Any}("edge_outer_polygon" => outer, "edge_cells" => [CP(outer, [hole])])
        ref = SMLMAnalysis.BasicSMLD(em, cam, 3, 1, Dict{String,Any}())
        chan = SMLMAnalysis.BasicSMLD(em, cam, 3, 1, deepcopy(md))
        aligned = [ref, chan]

        info = SMLMDriftCorrection.AlignInfo(
            [zeros(2), zeros(2)],
            SMLMDriftCorrection.AbstractAlignTransform[SMLMDriftCorrection.AffineTransform2D(0.0, 1.0, 0.0, 0.0),
                                                        SMLMDriftCorrection.AffineTransform2D(0.0, 1.0, 0.0, 0.0)],
            0.01, :entropy, :affine, :cpu, nothing)

        @test_logs (:warn, r"not carried through :affine") SMLMAnalysis._align_edge_geometry!(aligned, info)
        @test !haskey(aligned[2].metadata, "edge_outer_polygon")
        @test !haskey(aligned[2].metadata, "edge_cells")
        @test aligned[1].metadata == Dict{String,Any}()   # reference untouched
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
        mksmld(xs, ys, cam) = SMLMAnalysis.BasicSMLD(
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
        # Reproduces a real-data case that produced a=-42873, k≈0 with unbounded NelderMead.
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
        cropped = SMLMAnalysis.crop_images(img, roi_x, roi_y)
        @test size(cropped) == (length(roi_y), length(roi_x), nfr)
        @test cropped == img[roi_y, roi_x, :]
        @test cropped[1, 1, 1] == 1000 * first(roi_y) + 10 * first(roi_x) + 1

        # crop_camera uses the same convention: roi_x → x-edges, roi_y → y-edges.
        cam = IdealCamera(ncol, nrow, 0.1)   # IdealCamera(nx=cols, ny=rows, px)
        cc = SMLMAnalysis.crop_camera(cam, roi_x, roi_y)
        @test cc.pixel_edges_x == cam.pixel_edges_x[first(roi_x):last(roi_x)+1]
        @test cc.pixel_edges_y == cam.pixel_edges_y[first(roi_y):last(roi_y)+1]
        @test length(cc.pixel_edges_x) - 1 == length(roi_x)   # x pixel count = #cols
        @test length(cc.pixel_edges_y) - 1 == length(roi_y)   # y pixel count = #rows
    end

    @testset "roi refuses a step-level DetectFit camera" begin
        # roi crops only the pipeline camera; a camera on DetectFitConfig would stay
        # full-frame and offset every localization by the crop origin.
        cam = IdealCamera(16, 16, 0.1)
        cfg = AnalysisConfig(camera=cam, roi=(x=3:10, y=3:10),
                             steps=[DetectFitConfig(camera=cam)], outdir=nothing)
        # Match the message: all-zero frames also throw an (unrelated) ArgumentError downstream.
        @test_throws r"DetectFitConfig has its own camera" analyze(zeros(Float32, 16, 16, 2), cfg)
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
            smld = SMLMAnalysis.BasicSMLD(es, cam, 10, 1, Dict{String,Any}())
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
            smld_xy = SMLMAnalysis.BasicSMLD(exy, cam, 10, 1, Dict{String,Any}())
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
        #  - detectfit sets each emitter's dataset field by mutating it in place
        #    (all advertised emitter types are mutable structs, so this needs no
        #    per-type dispatch — unlike the old positional-reconstruction
        #    `_with_dataset`, removed once every type was confirmed mutable), for
        #    EVERY advertised GaussMLE emitter type (GaussianXYNB → Emitter2DFitGaussMLE,
        #    AstigmaticXYZNB → Emitter3DFitGaussMLE) as well as the standard ones.
        #  - save_smld/load_smld must round-trip those types and the standard-3D
        #    off-diagonal covariances, AND must key off the concrete emitter even when
        #    the SMLD is typed SMLMAnalysis.BasicSMLD{T,AbstractEmitter} (as the pre-narrowing
        #    pipeline produced) rather than degrading to SMLMAnalysis.Emitter2DFit and dropping σ.
        cam = IdealCamera(8, 8, 0.1)
        T = Float64

        # Setting .dataset in place must work for the fixed-width / astigmatic types
        # and leave every other field untouched.
        e2g = GaussMLE.Emitter2DFitGaussMLE{T}(0.1, 0.2, 1000.0, 5.0,
                0.01, 0.012, 0.003, 20.0, 0.5, 0.4, 1, 1, 0, 1)
        e3g = GaussMLE.Emitter3DFitGaussMLE{T}(0.1, 0.2, 0.3, 1000.0, 5.0,
                0.01, 0.012, 0.02, 0.003, 0.001, 0.002, 20.0, 0.5, 0.4, 1, 1, 0, 1)
        e2g.dataset = 7
        @test e2g.dataset == 7
        @test e2g.σ_xy ≈ 0.003
        e3g.dataset = 7
        @test e3g.dataset == 7
        @test e3g.σ_yz ≈ 0.002

        mktempdir() do dir
            # Emitter2DFitGaussMLE (GaussianXYNB) round-trip.
            g2 = [GaussMLE.Emitter2DFitGaussMLE{T}(0.1i, 0.2i, 1000.0 + i, 5.0,
                    0.01, 0.012, 0.003, 20.0, 0.5, 0.4, i, 1, 0, i) for i in 1:4]
            s2 = SMLMAnalysis.BasicSMLD(g2, cam, 10, 1, Dict{String,Any}())
            p2 = joinpath(dir, "g2.h5"); save_smld(p2, s2); l2 = load_smld(p2)
            @test eltype(l2.emitters) <: GaussMLE.Emitter2DFitGaussMLE
            for (a, b) in zip(s2.emitters, l2.emitters), f in fieldnames(GaussMLE.Emitter2DFitGaussMLE)
                @test getfield(a, f) ≈ getfield(b, f)
            end

            # Emitter3DFitGaussMLE (AstigmaticXYZNB) round-trip: z + full covariance.
            g3 = [GaussMLE.Emitter3DFitGaussMLE{T}(0.1i, 0.2i, 0.3i, 1000.0 + i, 5.0,
                    0.01, 0.012, 0.02, 0.003, 0.001, 0.002, 20.0, 0.5, 0.4, i, 1, 0, i) for i in 1:4]
            s3 = SMLMAnalysis.BasicSMLD(g3, cam, 10, 1, Dict{String,Any}())
            p3 = joinpath(dir, "g3.h5"); save_smld(p3, s3); l3 = load_smld(p3)
            @test eltype(l3.emitters) <: GaussMLE.Emitter3DFitGaussMLE
            @test l3.emitters[2].z ≈ 0.6
            @test l3.emitters[2].σ_xz ≈ 0.001
            @test l3.emitters[2].σ_yz ≈ 0.002

            # Standard SMLMAnalysis.Emitter3DFit: off-diagonal covariances σ_xz/σ_yz survive
            # (they were never written before this fix).
            e3 = [SMLMAnalysis.Emitter3DFit{T}(0.1i, 0.2i, 0.3i, 1000.0 + i, 5.0,
                    0.01, 0.012, 0.02, 20.0, 0.5;
                    σ_xy=0.003, σ_xz=0.001, σ_yz=0.002, frame=i, dataset=1, id=i) for i in 1:4]
            s3s = SMLMAnalysis.BasicSMLD(e3, cam, 10, 1, Dict{String,Any}())
            p3s = joinpath(dir, "e3.h5"); save_smld(p3s, s3s); l3s = load_smld(p3s)
            @test eltype(l3s.emitters) <: SMLMAnalysis.Emitter3DFit
            @test l3s.emitters[2].σ_xz ≈ 0.001
            @test l3s.emitters[2].σ_yz ≈ 0.002

            # Abstract-eltype SMLD (as the pipeline produced before narrowing):
            # save_smld must reload it as concrete Emitter2DFitSigma, NOT degrade to
            # SMLMAnalysis.Emitter2DFit and drop the PSF-width σ.
            abs_v = SMLMAnalysis.AbstractEmitter[GaussMLE.Emitter2DFitSigma{T}(
                        0.1i, 0.2i, 1000.0 + i, 5.0, 0.13,
                        0.01, 0.012, 0.003, 20.0, 0.5, 0.002,
                        0.4, i, 1, 0, i) for i in 1:4]
            @test eltype(abs_v) == SMLMAnalysis.AbstractEmitter
            s_abs = SMLMAnalysis.BasicSMLD(abs_v, cam, 10, 1, Dict{String,Any}())
            pabs = joinpath(dir, "abs.h5"); save_smld(pabs, s_abs); labs = load_smld(pabs)
            @test eltype(labs.emitters) <: GaussMLE.Emitter2DFitSigma   # NOT SMLMAnalysis.Emitter2DFit
            @test labs.emitters[2].σ ≈ 0.13                             # PSF-width σ preserved
            @test labs.emitters[2].σ_xy ≈ 0.003
        end
    end

    @testset "step checkpoint is versioned HDF5" begin
        # _save_step_smld writes through save_smld; load_smld must read it back unchanged.
        cam = IdealCamera(16, 16, 0.1)
        em = [SMLMAnalysis.Emitter2DFit{Float64}(0.1i, 0.2i, 1000.0 + i, 5.0, 0.01, 0.012, 20.0, 0.5;
                                    frame=i, dataset=1 + (i % 2)) for i in 1:6]
        smld = SMLMAnalysis.BasicSMLD(em, cam, 6, 2, Dict{String,Any}())
        dm = SMLMDriftCorrection.LegendrePolynomial(smld; degree=2)
        mktempdir() do dir
            p = SMLMAnalysis._save_step_smld(joinpath(dir, "03_driftcorrect"), smld;
                                             filename="smld_corrected.h5", drift_model=dm)
            @test p == joinpath(dir, "03_driftcorrect", "smld_corrected.h5") && isfile(p)
            s2 = load_smld(p)
            @test s2.emitters isa Vector{SMLMAnalysis.Emitter2DFit{Float64}}
            @test all(getfield(a, f) == getfield(b, f) for (a, b) in zip(em, s2.emitters)
                      for f in fieldnames(SMLMAnalysis.Emitter2DFit{Float64}))
            @test (s2.n_frames, s2.n_datasets) == (6, 2)
            @test s2.camera.pixel_edges_x == cam.pixel_edges_x
            @test s2.metadata["drift_correction"]["model_type"] == "LegendrePolynomial"
            @test SMLMAnalysis._save_step_smld(nothing, smld; filename="x.h5") === nothing
        end
    end

    @testset "edge geometry metadata round-trip" begin
        # Edge classification mirrors its cell mask into metadata; save_smld must keep it.
        CP = SMLMAnalysis.SMLMClustering.CellPolygon
        sq(x0, s) = NTuple{2,Float64}[(x0, x0), (x0 + s, x0), (x0 + s, x0 + s), (x0, x0 + s)]
        cells = [CP(sq(0.0, 4.0), [sq(0.5, 1.0), sq(2.0, 0.5)]),   # two holes
                 CP(sq(5.0, 1.0)),                                # no holes
                 CP(sq(7.0, 2.0), [sq(7.5, 0.2)])]
        outer = sq(0.0, 4.0)
        cam = IdealCamera(16, 16, 0.1)
        em = [SMLMAnalysis.Emitter2DFit{Float64}(0.1i, 0.1i, 1000.0, 5.0, 0.01, 0.01, 20.0, 0.5; frame=i) for i in 1:3]
        cellkey(cs) = [(c.outer, c.holes) for c in cs]   # CellPolygon has no ==; compare its fields
        mktempdir() do dir
            md = Dict{String,Any}("edge_cells" => cells, "edge_outer_polygon" => outer,
                                  "empty_cells" => CP[], "empty_polygon" => NTuple{2,Float64}[])
            p = joinpath(dir, "geom.h5")
            save_smld(p, SMLMAnalysis.BasicSMLD(em, cam, 3, 1, md))
            m2 = load_smld(p).metadata
            @test m2["edge_outer_polygon"] == outer
            @test m2["edge_outer_polygon"] isa Vector{NTuple{2,Float64}}
            @test m2["edge_cells"] isa Vector{CP}
            @test cellkey(m2["edge_cells"]) == cellkey(cells)
            @test m2["empty_cells"] isa Vector{CP} && isempty(m2["empty_cells"])
            @test m2["empty_polygon"] isa Vector{NTuple{2,Float64}} && isempty(m2["empty_polygon"])
        end
    end

    @testset "TOML provenance is valid" begin
        # Upstream field names outside TOML's bare-key set (DriftConfig.σ_loc) must be
        # written as quoted keys, or config.toml fails to parse.
        let io = IOBuffer()
            SMLMAnalysis._write_config_fields!(io, DriftConfig())
            parsed = TOML.parse(String(take!(io)))
            @test haskey(parsed, "σ_loc")
        end
        # Provenance files are named .toml and must parse back. The hand-rolled
        # serializer used to emit invalid TOML for tuples ((500.0, Inf)), ranges
        # (1:19), symbol vectors ([:red, :blue]), and unescaped strings. _toml_value
        # fixes every scalar-emit point; assert the round-trip actually works.
        mktempdir() do dir
            # FilterConfig: photons/precision are Tuple{Float64,Float64}; Inf must
            # serialize as TOML `inf`. FilterConfig has no nested-config fields, so
            # every value sits at the root table.
            cfg = FilterConfig(photons=(500.0, Inf), precision=(0.0, 0.02))
            SMLMAnalysis._save_config!(dir, cfg)
            parsed = TOML.parsefile(joinpath(dir, "config.toml"))   # throws if invalid TOML
            @test parsed["type"] == "FilterConfig"
            @test parsed["photons"][1] == 500.0
            @test isinf(parsed["photons"][2])        # `inf` parses back to Inf::Float64
            @test parsed["precision"] == [0.0, 0.02]

            # DetectFitConfig.datasets is an AbstractVector{Int}; a UnitRange (1:19)
            # used to emit the bare, invalid `datasets = 1:19`. It must now parse.
            cfg2 = DetectFitConfig(datasets=1:19)
            SMLMAnalysis._save_config!(dir, cfg2)
            @test TOML.parsefile(joinpath(dir, "config.toml")) isa AbstractDict  # no parse error

            # _toml_value unit behavior for the hard scalar types.
            @test SMLMAnalysis._toml_value(1:19) == "\"1:19\""       # range -> quoted, not materialized
            @test SMLMAnalysis._toml_value((500.0, Inf)) == "[500.0, inf]"
            @test SMLMAnalysis._toml_value([:red, :blue]) == "[\"red\", \"blue\"]"
            @test SMLMAnalysis._toml_value(true) == "true"

            # String escaping: embedded quote, backslash, and newline survive a
            # write -> TOML.parse round-trip (the old String branch escaped nothing).
            s = "he said \"hi\"\\\n done"
            open(joinpath(dir, "esc.toml"), "w") do io
                println(io, "s = $(SMLMAnalysis._toml_value(s))")
            end
            @test TOML.parsefile(joinpath(dir, "esc.toml"))["s"] == s

            # info.toml: a `nothing` field is OMITTED (as in config.toml), never written
            # as the string "nothing" — otherwise a key's TOML type flips between runs.
            info = (; n_after = 3, p2_estimate = nothing, field_mode = :gaussian)
            SMLMAnalysis._save_info!(dir, info)
            parsed_info = TOML.parsefile(joinpath(dir, "info.toml"))
            @test parsed_info["n_after"] == 3
            @test parsed_info["field_mode"] == "gaussian"
            @test !haskey(parsed_info, "p2_estimate")
        end
    end

    @testset "intensity-filter p₂ mixture estimator" begin
        # Regression for the unbiased double-emitter fraction estimator.
        # The mixture model must (a) recover a small KNOWN p₂ and, critically,
        # (b) return ≈0 when there are NO doubles — NOT the ~(1 - rate_percentile)
        # ≈ 5% floor that the legacy tail-ratio estimate is pinned at.
        rng = MersenneTwister(20260717)
        # Peaked, identifiable single-emitter model: Gamma(k=3, θ=1) drawn as the sum
        # of three exponentials. Field-normalized by the single-emitter p95 (→ ~1),
        # mirroring the pipeline's photons/λ(x,y).
        gamma3(r) = -log(rand(r)) - log(rand(r)) - log(rand(r))
        function make_normalized(p_true, N, r)
            n_dbl = round(Int, N * p_true)
            n_sgl = N - n_dbl
            raw = [gamma3(r) for _ in 1:(n_sgl + 2 * n_dbl)]
            scale = quantile(raw, 0.95)
            scale <= 0 && (scale = 1.0)
            singles = raw[1:n_sgl] ./ scale
            da = raw[n_sgl+1 : n_sgl+n_dbl] ./ scale
            db = raw[n_sgl+n_dbl+1 : n_sgl+2*n_dbl] ./ scale
            vcat(singles, da .+ db)   # f_double = f_single ⊛ f_single (sum of two singles)
        end
        cfg = IntensityFilterConfig(p2_n_bins = 200)
        @test cfg.p2_method === :mixture   # new behavior is the default

        # (a) p_true = 1% → recovered within a few × and clearly off the 5% floor.
        p_hat = SMLMAnalysis._estimate_p2_mixture(make_normalized(0.01, 60_000, rng), cfg)
        @test p_hat !== nothing
        @test 0.003 <= p_hat <= 0.03

        # (b) p_true = 0% → near zero, NOT pinned near 1 - rate_percentile (~0.05).
        p_hat0 = SMLMAnalysis._estimate_p2_mixture(make_normalized(0.0, 60_000, rng), cfg)
        @test p_hat0 !== nothing
        @test p_hat0 < 0.02
    end

    @testset "intensity-filter spatial binning matches the reference scan" begin
        # Reference implementation: origin/main's _spatial_bin_rates, an O(n_bins²·n)
        # per-bin scan against the exact half-open edges (see git show
        # origin/main:src/steps/intensityfilter.jl). SMLMAnalysis._spatial_bin_rates
        # computes each emitter's bin directly via _bin_index instead of scanning, so
        # this pins the single-pass version to the scan's exact edge comparisons —
        # including the float-rounding cases a floor-division guess gets wrong.
        function ref_spatial_bin_rates(xs, ys, photons, cfg)
            n_bins = cfg.n_bins
            min_count = cfg.min_bin_count
            pct = cfg.rate_percentile

            x_min, x_max = extrema(xs)
            y_min, y_max = extrema(ys)

            dx = (x_max - x_min) / n_bins
            dy = (y_max - y_min) / n_bins
            dx == 0 && (dx = 1.0)
            dy == 0 && (dy = 1.0)

            centers_x = Float64[]
            centers_y = Float64[]
            rates = Float64[]
            counts = Int[]
            rate_grid = fill(NaN, n_bins, n_bins)

            for ix in 1:n_bins, iy in 1:n_bins
                bx_lo = x_min + (ix - 1) * dx
                bx_hi = ix == n_bins ? x_max + eps(x_max) : x_min + ix * dx
                by_lo = y_min + (iy - 1) * dy
                by_hi = iy == n_bins ? y_max + eps(y_max) : y_min + iy * dy

                bin_photons = Float64[]
                for i in eachindex(xs)
                    if bx_lo <= xs[i] < bx_hi && by_lo <= ys[i] < by_hi
                        push!(bin_photons, photons[i])
                    end
                end

                if length(bin_photons) >= min_count
                    rate = quantile(bin_photons, pct)
                    push!(centers_x, (bx_lo + min(bx_hi, x_max)) / 2)
                    push!(centers_y, (by_lo + min(by_hi, y_max)) / 2)
                    push!(rates, rate)
                    push!(counts, length(bin_photons))
                    rate_grid[ix, iy] = rate
                end
            end

            x_edges = range(x_min, x_max, length=n_bins+1)
            y_edges = range(y_min, y_max, length=n_bins+1)

            return (centers_x=centers_x, centers_y=centers_y, rates=rates, counts=counts,
                    rate_grid=rate_grid, x_edges=x_edges, y_edges=y_edges)
        end

        # rate_grid carries NaN in unfilled bins; isequal (not ==) is the correct
        # "identical field-by-field" check since NaN == NaN is false but
        # isequal(NaN, NaN) is true — a bare == would fail even when the two
        # implementations agree exactly.
        function check_matches_reference(xs, ys, photons, cfg)
            ref = ref_spatial_bin_rates(xs, ys, photons, cfg)
            got = SMLMAnalysis._spatial_bin_rates(xs, ys, photons, cfg)
            @test isequal(ref, got)
        end

        # (a) Interior-edge float rounding: extrema (1.0, 2.0), 10 bins → dx=0.1;
        # x=1.2 lands in bin 3 under exact edge comparison, but a naive
        # floor((x-x_min)/dx)+1 guess can drift to bin 2 on this input.
        cfg_a = IntensityFilterConfig(n_bins=10, min_bin_count=1)
        xs_a = [1.0, 1.2, 2.0, 1.5, 1.3, 1.7, 1.9, 1.05, 1.85, 1.45]
        ys_a = collect(range(0.0, 1.0, length=10))
        photons_a = collect(100.0:100.0:1000.0)
        check_matches_reference(xs_a, ys_a, photons_a, cfg_a)
        @test SMLMAnalysis._bin_index(1.2, 1.0, 2.0, 0.1, 10) == 3

        # (b) Constant x: x_min == x_max forces dx = 1.0 (the /0 guard), which must
        # still clip bin centers to min(bx_hi, x_max).
        cfg_b = IntensityFilterConfig(n_bins=4, min_bin_count=1)
        xs_b = fill(5.0, 12)
        ys_b = collect(range(0.0, 3.0, length=12))
        photons_b = collect(10.0:10.0:120.0)
        check_matches_reference(xs_b, ys_b, photons_b, cfg_b)

        # (c) Constant y: same, on the other axis.
        xs_c = collect(range(0.0, 3.0, length=12))
        ys_c = fill(-2.0, 12)
        check_matches_reference(xs_c, ys_c, photons_b, cfg_b)

        # (d) Points exactly on interior edges x_min + k*dx (k = 0..n_bins), plus a
        # repeated x_max point ((e) the max point) to confirm the widened last-bin
        # edge still claims it rather than dropping it as out-of-range.
        cfg_d = IntensityFilterConfig(n_bins=10, min_bin_count=1)
        xs_d = vcat(collect(0.0:1.0:10.0), 10.0)   # 0,1,...,10, and a second 10.0
        ys_d = collect(range(0.0, 5.0, length=length(xs_d)))
        photons_d = collect(1.0:length(xs_d))
        check_matches_reference(xs_d, ys_d, photons_d, cfg_d)

        # (f) ~20 random seeded clouds.
        rng_bins = MersenneTwister(20260926)
        for _ in 1:20
            n = rand(rng_bins, 40:200)
            xs_r = rand(rng_bins, n) .* 10 .- 3
            ys_r = rand(rng_bins, n) .* 6 .+ 1
            photons_r = rand(rng_bins, n) .* 900 .+ 100
            cfg_r = IntensityFilterConfig(n_bins=rand(rng_bins, (4, 5, 8, 10)), min_bin_count=1)
            check_matches_reference(xs_r, ys_r, photons_r, cfg_r)
        end
    end

    @testset "stepinfo/stepinfos accessors" begin
        # Build StepInfos directly; FilterConfig → name "filter" (repeated),
        # DensityFilterConfig → name "densityfilter" (unique).
        si_a = SMLMAnalysis.StepInfo(1, FilterConfig(), 0.1, Dict{Symbol,Any}(); info=SMLMAnalysis.FilterInfo(100, 90, 0.1))
        si_b = SMLMAnalysis.StepInfo(2, DensityFilterConfig(), 0.2, Dict{Symbol,Any}(); info=SMLMAnalysis.DensityFilterInfo(90, 80, 5, 0.2))
        si_c = SMLMAnalysis.StepInfo(3, FilterConfig(), 0.3, Dict{Symbol,Any}(); info=SMLMAnalysis.FilterInfo(200, 150, 0.3))
        steps = SMLMAnalysis.StepInfo[si_a, si_b, si_c]

        for info in (SMLMAnalysis.AnalysisInfo(1.0, steps),
                     SMLMAnalysis.MultiTargetInfo(1.0, Dict{Symbol,SMLMAnalysis.AnalysisInfo}(), steps))
            # Symbol and String name lookup return the same SMLMAnalysis.StepInfo
            @test stepinfo(info, :densityfilter) === si_b
            @test stepinfo(info, "densityfilter") === si_b
            @test stepinfo(info, :densityfilter).info isa SMLMAnalysis.DensityFilterInfo

            # Repeated name: stepinfo returns the FIRST, stepinfos returns ALL in order
            @test stepinfo(info, :filter) === si_a
            @test stepinfo(info, "filter") === si_a
            @test stepinfos(info, :filter) == [si_a, si_c]
            @test [si.number for si in stepinfos(info, "filter")] == [1, 3]

            # Missing name: stepinfo throws ArgumentError listing available names,
            # stepinfos returns empty
            @test_throws ArgumentError stepinfo(info, :nope)
            @test isempty(stepinfos(info, :nope))
        end
    end

    @testset "install_agent_guide" begin
        # Arg validation
        @test_throws ArgumentError install_agent_guide(tool = :bogus)
        @test_throws ArgumentError install_agent_guide(scope = :bogus)
        @test_throws ArgumentError uninstall_agent_guide(tool = :bogus)
        @test_throws ArgumentError agent_guide_status(scope = :bogus)

        # Claude, project scope, default track=false → gitignored, namespaced, stamped.
        mktempdir() do dir
            skill = install_agent_guide(dir = dir)
            @test skill == joinpath(dir, ".claude", "skills", "smlma-ecosystem")
            @test isfile(joinpath(skill, "SKILL.md"))

            refs = readdir(joinpath(skill, "reference"))
            # SMLMAnalysis + the 10 ecosystem packages, each with a reference file.
            @test length(refs) == 11
            for name in ("SMLMAnalysis", "SMLMData", "GaussMLE", "SMLMBaGoL", "SMLMRender")
                @test "$name.md" in refs
            end

            skilltext = read(joinpath(skill, "SKILL.md"), String)
            @test occursin("name: smlma-ecosystem", skilltext)
            @test occursin("description:", skilltext)
            @test occursin("x-installer: SMLMAnalysis", skilltext)   # provenance stamp
            @test occursin("x-source-version:", skilltext)
            @test occursin("Dependency hierarchy", skilltext)

            # track=false (default) gitignores the namespaced skill dir.
            @test occursin(".claude/skills/smlma-ecosystem/", read(joinpath(dir, ".gitignore"), String))

            # The rewritten .gitignore follows the umask, like a direct `write`
            # would — not the fixed 0600 a mktemp-created temp carries over
            # (the bug _replace_atomically's private-temp-dir + write() path avoids).
            control = joinpath(dir, "control.txt")
            write(control, "control")
            @test filemode(joinpath(dir, ".gitignore")) & 0o777 == filemode(control) & 0o777

            # Doctor: freshly installed, not stale.
            st = agent_guide_status(dir = dir)
            @test st.installed
            @test !st.stale
            @test st.source_version == st.current_version

            # Own-install idempotency: re-running our OWN stamped install refreshes
            # WITHOUT overwrite (no error), returning the same path.
            @test install_agent_guide(dir = dir) == skill

            # A foreign/unstamped skill in the same dir IS refused unless overwrite.
            write(joinpath(skill, "SKILL.md"), "---\nname: smlma-ecosystem\n---\nhand-made\n")
            @test_throws ErrorException install_agent_guide(dir = dir)
            @test install_agent_guide(dir = dir, overwrite = true) == skill
            @test occursin("x-installer: SMLMAnalysis", read(joinpath(skill, "SKILL.md"), String))

            # Uninstall removes only our stamped install.
            @test uninstall_agent_guide(dir = dir) == [skill]
            @test !isdir(skill)
            @test !agent_guide_status(dir = dir).installed
            @test isempty(uninstall_agent_guide(dir = dir))   # nothing left → no-op
        end

        # Guard-gap regression: a hand-made target dir with a reference/ but NO
        # SKILL.md must be REFUSED (not silently wiped) unless overwrite=true.
        mktempdir() do dir
            skill = joinpath(dir, ".claude", "skills", "smlma-ecosystem")
            mkpath(joinpath(skill, "reference"))
            write(joinpath(skill, "reference", "keep.md"), "hand-made")
            @test_throws ErrorException install_agent_guide(dir = dir)
            @test isfile(joinpath(skill, "reference", "keep.md"))            # survived
            @test install_agent_guide(dir = dir, overwrite = true) == skill  # overwrite proceeds
            @test isfile(joinpath(skill, "SKILL.md"))
        end

        # Uninstall leaves a foreign skill untouched.
        mktempdir() do dir
            skill = joinpath(dir, ".claude", "skills", "smlma-ecosystem")
            mkpath(skill)
            write(joinpath(skill, "SKILL.md"), "---\nname: smlma-ecosystem\n---\nnot ours\n")
            @test isempty(uninstall_agent_guide(dir = dir))
            @test isdir(skill)
        end

        # Regression (data loss): overwrite=true onto a hand-made dir stamps it as ours
        # while preserving the user's own files — a later uninstall must remove ONLY
        # what we wrote (SKILL.md + reference/) and leave the still-populated directory
        # in place, never `rm -r` the user's notes.md / scripts/ along with it.
        mktempdir() do dir
            skill = joinpath(dir, ".claude", "skills", "smlma-ecosystem")
            mkpath(joinpath(skill, "scripts"))
            write(joinpath(skill, "SKILL.md"), "---\nname: smlma-ecosystem\n---\nhand-made\n")
            write(joinpath(skill, "notes.md"), "my notes")
            write(joinpath(skill, "scripts", "tool.py"), "print(1)")
            @test install_agent_guide(dir = dir, overwrite = true) == skill
            @test isfile(joinpath(skill, "notes.md"))                       # install preserved it
            @test uninstall_agent_guide(dir = dir) == [skill]
            @test isdir(skill)                                              # dir kept (non-empty)
            @test isfile(joinpath(skill, "notes.md"))                       # user's files survive
            @test isfile(joinpath(skill, "scripts", "tool.py"))
            @test !isfile(joinpath(skill, "SKILL.md"))                      # ours removed
            @test !isdir(joinpath(skill, "reference"))
            @test !agent_guide_status(dir = dir).installed
            @test isempty(uninstall_agent_guide(dir = dir))                 # nothing of ours left
        end

        # Same contract for the Codex bundle: user files inside smlm-agent-guide/ survive.
        mktempdir() do dir
            bundle = joinpath(dir, "smlm-agent-guide")
            mkpath(bundle)
            write(joinpath(bundle, "GUIDE.md"), "hand-made\n")
            write(joinpath(bundle, "mine.txt"), "keep")
            @test install_agent_guide(dir = dir, tool = :codex, overwrite = true) == bundle
            removed = uninstall_agent_guide(dir = dir, tool = :codex)
            @test bundle in removed
            @test isdir(bundle) && isfile(joinpath(bundle, "mine.txt"))
            @test !isfile(joinpath(bundle, "GUIDE.md")) && !isdir(joinpath(bundle, "reference"))
        end

        # reference/ survival across a refresh AND an uninstall (Claude): a file the
        # user adds inside reference/ is never touched, but files we generated
        # (identified by their per-file provenance stamp) ARE regenerated by a
        # refresh / removed by uninstall.
        mktempdir() do dir
            skill  = install_agent_guide(dir = dir)
            refdir = joinpath(skill, "reference")
            notes  = joinpath(refdir, "personal-notes.md")
            write(notes, "mine")

            datamd = joinpath(refdir, "SMLMData.md")
            rm(datamd)   # force a visible change so we can tell a refresh regenerates it

            @test install_agent_guide(dir = dir) == skill      # plain refresh, no overwrite
            @test read(notes, String) == "mine"                # user file survived untouched
            @test isfile(datamd)                               # ours was regenerated

            removed = uninstall_agent_guide(dir = dir)
            @test removed == [skill]
            @test isfile(notes)                                 # user file still there
            @test isdir(refdir)                                 # dir kept (notes still inside)
            @test !isfile(datamd)                               # generated files removed
            @test !isfile(joinpath(skill, "SKILL.md"))
            @test !agent_guide_status(dir = dir).installed
        end

        # Same contract for Codex.
        mktempdir() do dir
            bundle = install_agent_guide(dir = dir, tool = :codex)
            refdir = joinpath(bundle, "reference")
            notes  = joinpath(refdir, "personal-notes.md")
            write(notes, "mine")

            datamd = joinpath(refdir, "SMLMData.md")
            rm(datamd)

            @test install_agent_guide(dir = dir, tool = :codex) == bundle
            @test read(notes, String) == "mine"
            @test isfile(datamd)

            removed = uninstall_agent_guide(dir = dir, tool = :codex)
            @test bundle in removed
            @test isfile(notes)
            @test isdir(refdir)
            @test !isfile(datamd)
            @test !isfile(joinpath(bundle, "GUIDE.md"))
            @test !agent_guide_status(dir = dir, tool = :codex).installed
        end

        # Data-loss regression: a DIRECTORY (not just a symlink) sitting at a reference
        # file's name must be refused, never wiped by mv(...; force=true)'s implicit
        # recursive rm of the destination. Exercised via :claude; :codex shares the same
        # _preflight_reference_files helper.
        mktempdir() do dir
            skill  = install_agent_guide(dir = dir)
            refdir = joinpath(skill, "reference")
            datamd = joinpath(refdir, "SMLMData.md")
            rm(datamd)
            mkpath(datamd)
            write(joinpath(datamd, "keep.txt"), "keep")

            @test_throws ArgumentError install_agent_guide(dir = dir)
            @test isdir(datamd)                                 # not wiped
            @test isfile(joinpath(datamd, "keep.txt"))           # contents survived
        end

        # Same regression for AGENTS.md itself: a directory there must be refused
        # BEFORE any bundle content is written (nothing half-installed).
        mktempdir() do dir
            mkpath(joinpath(dir, "AGENTS.md"))
            @test_throws ArgumentError install_agent_guide(dir = dir, tool = :codex)
            @test !ispath(joinpath(dir, "smlm-agent-guide"))    # bundle never written
        end

        # agent_guide's writes now go through _replace_atomically (src/io/atomic.jl),
        # same as save_smld: a failed write (destination is a directory, refused
        # outright) must leave no stray temp file behind.
        mktempdir() do dir
            mkpath(joinpath(dir, "somedir"))
            @test_throws ArgumentError SMLMAnalysis._replace_atomically(joinpath(dir, "somedir")) do tmp
                write(tmp, "x")
            end
            @test readdir(dir) == ["somedir"]   # no stray temp file
        end

        if Sys.isunix()
            # (a) target itself is a symlink to a real, stamped install elsewhere:
            # neither install nor uninstall may follow it (Claude).
            mktempdir() do root
                actual_repo = joinpath(root, "actual_repo")
                mkpath(actual_repo)
                actual_target = install_agent_guide(dir = actual_repo)

                repo = joinpath(root, "repo")
                mkpath(joinpath(repo, ".claude", "skills"))
                target = joinpath(repo, ".claude", "skills", "smlma-ecosystem")
                symlink(actual_target, target; dir_target = true)

                @test isempty(uninstall_agent_guide(dir = repo))
                @test isfile(joinpath(actual_target, "SKILL.md"))   # real install untouched
                @test_throws ArgumentError install_agent_guide(dir = repo)
            end

            # (a) same for Codex.
            mktempdir() do root
                actual_repo = joinpath(root, "actual_repo")
                mkpath(actual_repo)
                actual_bundle = install_agent_guide(dir = actual_repo, tool = :codex)

                repo = joinpath(root, "repo")
                mkpath(repo)
                target = joinpath(repo, "smlm-agent-guide")
                symlink(actual_bundle, target; dir_target = true)

                @test isempty(uninstall_agent_guide(dir = repo, tool = :codex))
                @test isfile(joinpath(actual_bundle, "GUIDE.md"))   # real install untouched
                @test_throws ArgumentError install_agent_guide(dir = repo, tool = :codex)
            end

            # (b) the wrapper itself is a symlink pointing outside the install dir
            # (Claude): install must refuse rather than write through it.
            mktempdir() do dir
                skill = joinpath(dir, ".claude", "skills", "smlma-ecosystem")
                mkpath(skill)
                external = joinpath(dir, "external.md")
                write(external, "not ours")
                symlink(external, joinpath(skill, "SKILL.md"))

                @test_throws ArgumentError install_agent_guide(dir = dir, overwrite = true)
                @test read(external, String) == "not ours"   # external file untouched
            end

            # (c) a reference file itself is a symlink to an external file: a plain
            # refresh must refuse rather than write through it.
            mktempdir() do dir
                skill    = install_agent_guide(dir = dir)
                refdir   = joinpath(skill, "reference")
                external = joinpath(dir, "external-ref.md")
                write(external, "not ours")
                rm(joinpath(refdir, "SMLMData.md"))
                symlink(external, joinpath(refdir, "SMLMData.md"))

                @test_throws ArgumentError install_agent_guide(dir = dir)
                @test read(external, String) == "not ours"   # external file untouched
            end

            # (d) the user replaces a generated reference file with their own content
            # (no provenance header): refresh without overwrite refuses and leaves it
            # untouched; overwrite=true replaces it and restores the stamp.
            mktempdir() do dir
                skill  = install_agent_guide(dir = dir)
                refdir = joinpath(skill, "reference")
                datamd = joinpath(refdir, "SMLMData.md")
                write(datamd, "hand-edited, no stamp\n")

                @test_throws ErrorException install_agent_guide(dir = dir)
                @test read(datamd, String) == "hand-edited, no stamp\n"

                @test install_agent_guide(dir = dir, overwrite = true) == skill
                @test occursin("SMLMAnalysis.install_agent_guide()", read(datamd, String))
            end

            # (e) a hardlinked wrapper: install must replace it via rename, never
            # truncate it in place, so the OTHER hardlink (the user's) is untouched.
            mktempdir() do dir
                skill    = joinpath(dir, ".claude", "skills", "smlma-ecosystem")
                mkpath(skill)
                external = joinpath(dir, "external-skill.md")
                write(external, "not ours")
                @test ccall(:link, Cint, (Cstring, Cstring), external, joinpath(skill, "SKILL.md")) == 0

                @test install_agent_guide(dir = dir, overwrite = true) == skill
                @test occursin("x-installer: SMLMAnalysis", read(joinpath(skill, "SKILL.md"), String))
                @test read(external, String) == "not ours"   # other hardlink untouched
            end

            # (f) symlinked .gitignore / AGENTS.md: never written through.
            mktempdir() do dir
                external = joinpath(dir, "external-gitignore")
                write(external, "external\n")
                symlink(external, joinpath(dir, ".gitignore"))
                @test_throws ArgumentError install_agent_guide(dir = dir)
                @test read(external, String) == "external\n"
            end
            mktempdir() do dir
                external = joinpath(dir, "external-agents")
                write(external, "external\n")
                symlink(external, joinpath(dir, "AGENTS.md"))
                @test_throws ArgumentError install_agent_guide(dir = dir, tool = :codex)
                @test read(external, String) == "external\n"
            end

            # (g) Codex: the bundle directory is replaced by a symlink to a copy
            # elsewhere — uninstall must not follow it, and must leave the AGENTS.md
            # block in place (it must not go on to strip the block once the bundle
            # itself is judged unsafe).
            mktempdir() do dir
                bundle         = install_agent_guide(dir = dir, tool = :codex)
                copy_elsewhere = joinpath(dir, "bundle-copy")
                cp(bundle, copy_elsewhere)
                rm(bundle; recursive = true)
                symlink(copy_elsewhere, bundle; dir_target = true)

                @test isempty(uninstall_agent_guide(dir = dir, tool = :codex))
                agents = read(joinpath(dir, "AGENTS.md"), String)
                @test occursin("BEGIN SMLMAnalysis agent-guide", agents)
            end

            # (h) a symlinked reference/ dir inside an otherwise-stamped install: the
            # doctor reports not-installed, and uninstall leaves everything untouched.
            mktempdir() do dir
                skill          = install_agent_guide(dir = dir)
                refdir         = joinpath(skill, "reference")
                copy_elsewhere = joinpath(dir, "reference-copy")
                cp(refdir, copy_elsewhere)
                rm(refdir; recursive = true)
                symlink(copy_elsewhere, refdir; dir_target = true)

                @test !agent_guide_status(dir = dir).installed
                @test isempty(uninstall_agent_guide(dir = dir))
                @test isfile(joinpath(skill, "SKILL.md"))
            end

            # (i) Ordering regression: a symlinked .gitignore must be refused BEFORE
            # the skill bundle is written — never a half-installed guide.
            mktempdir() do dir
                external = joinpath(dir, "external-gitignore-order")
                write(external, "external\n")
                symlink(external, joinpath(dir, ".gitignore"))
                @test_throws ArgumentError install_agent_guide(dir = dir)
                @test !ispath(joinpath(dir, ".claude", "skills", "smlma-ecosystem", "SKILL.md"))
            end

            # (j) Same ordering regression for a symlinked AGENTS.md (Codex): GUIDE.md
            # must never appear.
            mktempdir() do dir
                external = joinpath(dir, "external-agents-order")
                write(external, "external\n")
                symlink(external, joinpath(dir, "AGENTS.md"))
                @test_throws ArgumentError install_agent_guide(dir = dir, tool = :codex)
                @test !isfile(joinpath(dir, "smlm-agent-guide", "GUIDE.md"))
            end
        end

        # _expand_home strips repeated leading separators after "~/" too.
        @test SMLMAnalysis._expand_home("~//repo") == joinpath(homedir(), "repo")
        @test SMLMAnalysis._expand_home("~") == homedir()

        # AGENTS.md verbatim restore: content appended after our block (including
        # indentation) must survive an uninstall untouched — not run through strip().
        mktempdir() do dir
            write(joinpath(dir, "AGENTS.md"), "# My project rules\n\nBe careful.\n")
            install_agent_guide(dir = dir, tool = :codex)
            agents_path = joinpath(dir, "AGENTS.md")
            write(agents_path, read(agents_path, String) * "\n    user_code()\n")

            uninstall_agent_guide(dir = dir, tool = :codex)
            final = read(agents_path, String)
            @test occursin("\n    user_code()\n", final)   # indentation intact
            @test occursin("My project rules", final)
            @test !occursin("BEGIN SMLMAnalysis agent-guide", final)
        end

        # `dir = "~/..."` is expanded; it must never create a literal "~" under the cwd,
        # and the doctor's resolved path must actually point under the real home.
        mktempdir() do tmp
            cd(tmp) do
                st = agent_guide_status(dir = "~/smlma-nonexistent-repo-for-test")
                @test !st.installed
                @test st.path == joinpath(homedir(), "smlma-nonexistent-repo-for-test",
                                           ".claude", "skills", "smlma-ecosystem")
                @test !ispath("~")
            end
        end

        # Regression: the OLD (unexpanded) code passed a no-mutation check like the one
        # above too, so assert against a resolved path under a TEMPORARY HOME instead —
        # isolated from the real one — to actually exercise `~` expansion.
        if Sys.isunix()
            mktempdir() do home
                withenv("HOME" => home) do
                    mktempdir() do cwd
                        cd(cwd) do
                            p = install_agent_guide(dir = "~/repo")
                            @test startswith(p, joinpath(home, "repo"))
                            @test isfile(joinpath(p, "SKILL.md"))
                            @test !ispath("~")
                        end
                    end
                end
            end
        end

        # Claude, track=true → committed (no .gitignore written).
        mktempdir() do dir
            install_agent_guide(dir = dir, track = true)
            @test !isfile(joinpath(dir, ".gitignore"))
        end

        # Codex: stamped bundle + a managed block appended to AGENTS.md that preserves
        # pre-existing content and is idempotent; uninstall strips it back out.
        mktempdir() do dir
            write(joinpath(dir, "AGENTS.md"), "# My project rules\n\nBe careful.\n")
            bundle = install_agent_guide(dir = dir, tool = :codex)
            @test bundle == joinpath(dir, "smlm-agent-guide")
            @test isfile(joinpath(bundle, "GUIDE.md"))
            @test occursin("x-installer: SMLMAnalysis", read(joinpath(bundle, "GUIDE.md"), String))
            @test length(readdir(joinpath(bundle, "reference"))) == 11

            agents = read(joinpath(dir, "AGENTS.md"), String)
            @test occursin("My project rules", agents)                   # user content kept
            @test occursin("BEGIN SMLMAnalysis agent-guide", agents)     # our block added
            @test occursin("smlm-agent-guide/GUIDE.md", agents)

            # Idempotent refresh of our own bundle (no overwrite needed).
            install_agent_guide(dir = dir, tool = :codex)
            agents2 = read(joinpath(dir, "AGENTS.md"), String)
            @test count("BEGIN SMLMAnalysis agent-guide", agents2) == 1  # not duplicated
            @test occursin("My project rules", agents2)

            # Uninstall removes bundle + our block, preserving the user's content.
            removed = uninstall_agent_guide(dir = dir, tool = :codex)
            @test joinpath(dir, "smlm-agent-guide") in removed
            @test !isdir(bundle)
            agents3 = read(joinpath(dir, "AGENTS.md"), String)
            @test occursin("My project rules", agents3)
            @test !occursin("BEGIN SMLMAnalysis agent-guide", agents3)
        end
    end

    @testset "lab convention conformance" begin
        # Explicit conformance check for the lab skills-installer convention
        # (independent per-package implementations; see the convention doc). Each
        # assert maps to one spec invariant so the lab-guide can cite this block.
        mktempdir() do dir
            skill = install_agent_guide(dir = dir)                    # tool=:claude default
            # 1. Namespaced install dir: <pkgprefix>-<skill>
            @test occursin(r"[/\\]smlma-ecosystem$", skill)
            fm = read(joinpath(skill, "SKILL.md"), String)
            # 2. Provenance stamp: all four x- fields present
            for k in ("x-installer:", "x-source-version:", "x-source-commit:", "x-installed-format:")
                @test occursin(k, fm)
            end
            # 3. copy-never-symlink: installed files are real files, not links
            @test !islink(joinpath(skill, "SKILL.md"))
            # 4. track=false default anchors a .gitignore entry
            @test occursin("/.claude/skills/smlma-ecosystem/", read(joinpath(dir, ".gitignore"), String))
            # 5. own-install refresh is idempotent (no flag, same path)
            @test install_agent_guide(dir = dir) == skill
            # 6. stamp-scoped uninstall removes our own install
            @test uninstall_agent_guide(dir = dir) == [skill]
        end
    end

    @testset "MIC H5 calibration units" begin
        # SMITE convention: RawData = Gain_stored*photons + Offset (ADU); CCDVar is the
        # dark variance in ADU². Regression coverage for the ADU→e⁻ readnoise conversion
        # (load_mic_h5_calibration_for_scmos must divide by the stored ADU/e⁻ gain, not
        # just take sqrt(variance)).
        mktempdir() do dir
            path = joinpath(dir, "cal.h5")
            offset_adu = fill(100.0, 4, 4)
            gain_stored = fill(2.0, 4, 4)  # ADU/e⁻
            var_adu2 = fill(4.0, 4, 4)     # ADU² → readnoise = sqrt(4)/2 = 1.0 e⁻ rms
            SMLMAnalysis.HDF5.h5open(path, "w") do f
                g = SMLMAnalysis.HDF5.create_group(f, "Calibration")
                g["CCDOffset"] = offset_adu
                g["CCDVar"] = var_adu2
                g["Gain"] = gain_stored
            end

            cal = SMLMAnalysis.load_mic_h5_calibration_for_scmos(path)
            @test all(cal.readnoise .≈ 1.0f0)
            @test all(cal.gain .≈ 0.5f0)     # e⁻/ADU = 1/gain_stored
            @test all(cal.offset .≈ 100.0f0)

            cam = build_camera_from_mic_h5(path; pixel_size=0.1)
            @test all(cam.readnoise .≈ 1.0f0)
            @test all(cam.gain .≈ 0.5f0)
            @test all(cam.offset .≈ 100.0f0)
        end
    end

    @testset "load_mic_h5 warns on an unreadable data block" begin
        # Data002 is a group without its nested dataset, so _resolve_data_path
        # throws for it — load_mic_h5 must @warn (naming the block and the error)
        # and skip it, rather than the previous bare `catch; continue`.
        mktempdir() do dir
            path = joinpath(dir, "corrupt.h5")
            SMLMAnalysis.HDF5.h5open(path, "w") do f
                g = SMLMAnalysis.HDF5.create_group(f, "Channel01/Zposition001")
                g["Data001"] = rand(Float32, 4, 4, 3)
                SMLMAnalysis.HDF5.create_group(g, "Data002")
            end
            @test_logs (:warn, r"load_mic_h5: skipping unreadable data block \"Data002\"") match_mode=:any load_mic_h5(path)
            images, dataset_indices = load_mic_h5(path)
            @test size(images, 3) == 3
            @test all(==(1), dataset_indices)
        end
    end

    @testset "atomic saves refuse a directory destination" begin
        # save_smld goes through _replace_atomically: temp file + rename(2). rename
        # fails on a directory destination rather than deleting it (and everything
        # inside it), so it must throw and leave the directory intact. The exception
        # type differs by Julia version, so only a throw is asserted.
        cam = IdealCamera(8, 8, 0.1)
        T = Float64
        es = [SMLMAnalysis.Emitter2DFit{T}(0.1i, 0.2i, 1000.0 + i, 5.0, 0.01, 0.012, 20.0, 0.5, 0.4, i, 1, 0, i)
              for i in 1:3]
        smld = SMLMAnalysis.BasicSMLD(es, cam, 10, 1, Dict{String,Any}())

        mktempdir() do dir
            # save_smld: destination is a directory, not a file.
            target = joinpath(dir, "out.h5")
            mkdir(target)
            inner = joinpath(target, "keepme.txt")
            write(inner, "do not delete me")

            @test_throws Exception save_smld(target, smld)
            @test isdir(target)
            @test isfile(inner)
            @test read(inner, String) == "do not delete me"
            @test isempty(filter(f -> f != "keepme.txt", readdir(target)))

            # A normal save (no pre-existing directory) still works and leaves no
            # stray temp files behind.
            good = joinpath(dir, "good.h5")
            save_smld(good, smld)
            @test isfile(good)
            @test isempty(filter(f -> f ∉ ("out.h5", "good.h5"), readdir(dir)))
        end

        # Race: a directory appears at the destination after the isdir pre-check.
        # rename(2) itself must refuse it (never a recursive delete), and the temp
        # must be cleaned up. This is the property the pre-check cannot prove.
        mktempdir() do dir
            p = joinpath(dir, "out.h5")
            @test_throws Exception SMLMAnalysis._replace_atomically(p) do tmp
                write(tmp, "x")
                mkdir(p)
                write(joinpath(p, "keep"), "k")
            end
            @test read(joinpath(p, "keep"), String) == "k"
            @test readdir(dir) == ["out.h5"]
        end

        # The writer's temp path lives in a directory only we can enter (POSIX), so
        # no other user sharing the destination directory can pre-place a symlink at it.
        if !Sys.iswindows()
            mktempdir() do dir
                p = joinpath(dir, "out.txt")
                SMLMAnalysis._replace_atomically(p) do tmp
                    @test dirname(tmp) != dir
                    @test dirname(dirname(tmp)) == dir
                    @test filemode(dirname(tmp)) & 0o777 == 0o700
                    @test !ispath(tmp)
                    write(tmp, "x")
                end
                @test read(p, String) == "x"
                @test readdir(dir) == ["out.txt"]
            end
        end

        # A failing writer leaves no temp behind and never touches the destination.
        mktempdir() do dir
            p = joinpath(dir, "out.h5")
            @test_throws ErrorException SMLMAnalysis._replace_atomically(tmp -> error("boom"), p)
            @test !ispath(p)
            @test isempty(readdir(dir))

            write(p, "old content")
            @test_throws ErrorException SMLMAnalysis._replace_atomically(p) do tmp
                write(tmp, "partial")
                error("boom")
            end
            @test read(p, String) == "old content"
            @test readdir(dir) == ["out.h5"]
        end

        if !Sys.iswindows()
            # A symlink destination is replaced as a directory entry, not written through.
            mktempdir() do dir
                target = joinpath(dir, "target.h5")
                write(target, "target content")
                link = joinpath(dir, "link.h5")
                symlink(target, link)
                save_smld(link, smld)
                @test !islink(link)
                @test isfile(link)
                @test length(load_smld(link).emitters) == 3
                @test read(target, String) == "target content"
            end

            # The saved file's mode follows the umask, as a direct create would
            # (not the 0600 a mktemp-created temp would carry over).
            mktempdir() do dir
                p = joinpath(dir, "mode.h5")
                mask = UInt32(0o022)
                old = ccall(:umask, UInt32, (UInt32,), mask)
                try
                    save_smld(p, smld)
                finally
                    ccall(:umask, UInt32, (UInt32,), old)
                end
                @test filemode(p) & 0o777 == 0o666 & ~mask
            end
        end
    end

    @testset "upstream API smoke (tiny data)" begin
        # Runs the full detect/fit -> filter -> frame-connect -> drift -> render
        # pipeline through analyze(), at default verbosity so the figure/stats
        # writers (CairoMakie) run too, on data tiny enough to stay in the fast
        # tier. Its purpose is compat detection: each upstream API gets exercised
        # cheaply on every CI run, not just in the local thorough tier.
        Random.seed!(1)
        cam = IdealCamera(32, 32, 0.1)
        sim = SMLMAnalysis.StaticSMLMConfig(density = 5.0, σ_psf = 0.13, nframes = 50, ndatasets = 2)
        (_, si) = SMLMAnalysis.simulate(sim;
            pattern  = SMLMAnalysis.Nmer2D(n = 8, d = 0.05),
            molecule = SMLMAnalysis.GenericFluor(photons = 5.0e4, k_off = 20.0, k_on = 0.04),
            camera   = cam)
        images = [SMLMAnalysis.gen_images(si.smld_model, SMLMAnalysis.MicroscopePSFs.GaussianPSF(0.13);
                              dataset = d, bg = 20.0, poisson_noise = true)[1] for d in 1:2]

        cfg = AnalysisConfig(
            DetectFitConfig(boxer  = BoxerConfig(boxsize = 7, psf_sigma = 0.13, backend = :cpu),
                            fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu)),
            FilterConfig(photons = (100.0, Inf)),
            FrameConnectConfig(max_frame_gap = 2),
            DriftConfig(degree = 1),
            RenderConfig(zoom = 5);
            camera = cam,
            outdir = mktempdir(),
        )
        t = @elapsed (result, info) = analyze(images, cfg)
        @info "upstream API smoke (tiny data) wall time" seconds=t

        @test result isa SMLMAnalysis.AnalysisResult
        @test length(result.smld.emitters) >= 1

        # stepinfo: by name, by config type, and the not-found ArgumentError.
        @test stepinfo(info, :driftcorrect).name == "driftcorrect"
        @test stepinfo(info, DriftConfig).name == "driftcorrect"
        @test_throws ArgumentError stepinfo(info, :nosuchstep)
        @test_throws ArgumentError stepinfo(info, SMLMAnalysis.CrossAlignConfig)
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
            sim = SMLMAnalysis.StaticSMLMConfig(density = 5.0, σ_psf = 0.13, nframes = 50, ndatasets = 1)
            (_, si) = SMLMAnalysis.simulate(sim;
                pattern  = SMLMAnalysis.Nmer2D(n = 8, d = 0.05),
                molecule = SMLMAnalysis.GenericFluor(photons = 5.0e4, k_off = 20.0, k_on = 0.04),
                camera   = cam)
            (imgs, _) = SMLMAnalysis.gen_images(si.smld_model, SMLMAnalysis.MicroscopePSFs.GaussianPSF(0.13);
                dataset = 1, bg = 20.0, poisson_noise = true)

            cfg = AnalysisConfig(
                DetectFitConfig(boxer  = BoxerConfig(boxsize = 7, psf_sigma = 0.13, backend = :cpu),
                                fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu)),
                FilterConfig(photons = (100.0, Inf)),
                RenderConfig(zoom = 10);
                camera  = cam,
                verbose = Verbosity.SILENT,
            )
            (result, info) = analyze([imgs], cfg)

            @test result isa SMLMAnalysis.AnalysisResult
            @test info isa SMLMAnalysis.AnalysisInfo
            @test result.smld isa SMLMAnalysis.BasicSMLD
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
                                boxer  = BoxerConfig(boxsize = 7, psf_sigma = 0.13, backend = :cpu),
                                fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu));
                verbose = Verbosity.SILENT)
            @test si2.info.n_datasets == 2
            @test smld2.n_datasets == 2
            @test smld2.n_frames == size(imgs, 3)          # equal-length → per-dataset count
            @test Set(e.dataset for e in smld2.emitters) == Set([1, 2])
            @test all(1 <= e.frame <= size(imgs, 3) for e in smld2.emitters)  # frames are per-dataset
        end

        @testset "end-to-end with outdir" begin
            # Moderate simulated data through the full pipeline with disk output,
            # checking the on-disk contract: numbered step directories, valid TOML
            # provenance, every saved SMLD loadable and non-empty, and a render PNG.
            Random.seed!(1)
            cam = IdealCamera(64, 64, 0.1)
            sim = SMLMAnalysis.StaticSMLMConfig(density = 5.0, σ_psf = 0.13, nframes = 500, ndatasets = 2)
            (_, si) = SMLMAnalysis.simulate(sim;
                pattern  = SMLMAnalysis.Nmer2D(n = 8, d = 0.05),
                molecule = SMLMAnalysis.GenericFluor(photons = 5.0e4, k_off = 20.0, k_on = 0.04),
                camera   = cam)
            images = [SMLMAnalysis.gen_images(si.smld_model, SMLMAnalysis.MicroscopePSFs.GaussianPSF(0.13);
                                  dataset = d, bg = 20.0, poisson_noise = true)[1] for d in 1:2]

            outdir = mktempdir()
            cfg = AnalysisConfig(
                DetectFitConfig(boxer  = BoxerConfig(boxsize = 7, psf_sigma = 0.13, backend = :cpu),
                                fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu)),
                FilterConfig(photons = (100.0, Inf)),
                FrameConnectConfig(max_frame_gap = 2),
                DriftConfig(degree = 1),
                RenderConfig(zoom = 5);
                camera = cam,
                outdir = outdir,
            )
            (result, info) = analyze(images, cfg)

            # Numbered step directories only -- outdir also holds a top-level
            # .cache/ (checkpoint cache, see common.jl cache_dir) that isn't a step.
            step_dirs = filter(f -> isdir(f) && occursin(r"^\d\d_", basename(f)),
                                readdir(outdir; join=true))
            @test length(step_dirs) == length(cfg.steps)

            tomls = String[]
            h5s = String[]
            pngs = String[]
            for (root, _, files) in walkdir(outdir)
                for f in files
                    endswith(f, ".toml") && push!(tomls, joinpath(root, f))
                    endswith(f, ".h5")   && push!(h5s, joinpath(root, f))
                    endswith(f, ".png")  && push!(pngs, joinpath(root, f))
                end
            end
            @test !isempty(tomls)
            for p in tomls
                @test TOML.parsefile(p) isa Dict
            end

            smld_h5s = filter(p -> occursin("smld_", basename(p)), h5s)
            @test !isempty(smld_h5s)
            last_n = 0
            for p in smld_h5s
                loaded = load_smld(p)
                @test length(loaded.emitters) > 0
                if occursin("smld_corrected.h5", p)
                    last_n = length(loaded.emitters)
                end
            end
            @test last_n > 0
            @test length(result.smld.emitters) == last_n

            @test !isempty(pngs)
        end

        @testset "multi-target orchestrator writes smld_<label>.h5 through _finalize_channels!" begin
            # A multi-target channel is always raw images (or a file path) that goes
            # through its own DetectFitConfig -- analyze(channels, MultiTargetConfig)
            # has no path for already-localized data (steps=[] leaves _run_pipeline's
            # state a raw image Vector, never a SMLMAnalysis.BasicSMLD, so it errors "Pipeline
            # produced no SMLD"). So this drives a tiny simulated 2-channel run
            # end-to-end -- the multi-target analogue of "upstream API smoke" above --
            # to catch a regression where the orchestrator stops calling
            # _finalize_channels! after the phase-2 multi-target steps (Codex #51).
            # Channel B is channel A's own pattern (same seed) offset by a known
            # 0.1 μm, so the saved file must show the *aligned* B, not raw B.
            cam = IdealCamera(32, 32, 0.1)
            dx_true = 0.1
            gen_channel(seed, dx, dy) = begin
                Random.seed!(seed)
                sim = SMLMAnalysis.StaticSMLMConfig(density = 5.0, σ_psf = 0.13, nframes = 50, ndatasets = 1)
                (_, si) = SMLMAnalysis.simulate(sim;
                    pattern  = SMLMAnalysis.Nmer2D(n = 8, d = 0.05),
                    molecule = SMLMAnalysis.GenericFluor(photons = 5.0e4, k_off = 20.0, k_on = 0.04),
                    camera   = cam)
                model = deepcopy(si.smld_model)
                for e in model.emitters
                    e.x += dx
                    e.y += dy
                end
                (imgs, _) = SMLMAnalysis.gen_images(model, SMLMAnalysis.MicroscopePSFs.GaussianPSF(0.13);
                    dataset = 1, bg = 20.0, poisson_noise = true)
                [imgs]
            end
            images_a = gen_channel(11, 0.0, 0.0)
            images_b = gen_channel(11, dx_true, 0.0)

            chan_cfg() = AnalysisConfig(
                DetectFitConfig(boxer  = BoxerConfig(boxsize = 7, psf_sigma = 0.13, backend = :cpu),
                                fitter = GaussMLEConfig(psf_model = GaussianXYNBS(), backend = :cpu)),
                FilterConfig(photons = (100.0, Inf)),
                FrameConnectConfig(max_frame_gap = 2);
                camera = cam,
            )

            outdir = mktempdir()
            mt = MultiTargetConfig(
                labels = [:A, :B],
                steps  = [CrossAlignConfig(align = AlignConfig(method = :fft))],
                outdir = outdir,
            )
            (result, info) = analyze([(images_a, chan_cfg()), (images_b, chan_cfg())], mt)

            @test result isa SMLMAnalysis.MultiTargetResult
            for label in (:A, :B)
                p = joinpath(outdir, "smld_$(label).h5")
                @test isfile(p)
                loaded = load_smld(p)
                @test [e.x for e in loaded.emitters] == [e.x for e in result[label].smld.emitters]
                @test [e.y for e in loaded.emitters] == [e.y for e in result[label].smld.emitters]
            end

            # smld_B.h5 == result.smlds[2] == result[:B].smld, all three agreeing on
            # the SAME (aligned) coordinates.
            loaded_b = load_smld(joinpath(outdir, "smld_B.h5"))
            @test [e.x for e in loaded_b.emitters] == [e.x for e in result.smlds[2].emitters]
            @test [e.y for e in loaded_b.emitters] == [e.y for e in result.smlds[2].emitters]
            @test result.smlds[2] === result[:B].smld

            # And the aligned B must actually differ from B's pre-alignment data
            # -- catches a regression where the "aligned" save silently falls
            # back to writing unaligned data. `result[:B].smld_connected` is
            # NOT a safe baseline for this: it's FrameConnectInfo's pre-COMBINE
            # linked data, and combining alone (independent of cross-align)
            # shifts the mean by its own ~0.2 μm, which would make this
            # assertion pass even with a broken cross-align. Instead, compare
            # against the channel's own last-saved pre-alignment SMLD -- the
            # frame-connect step's `smld_combined.h5`, i.e. B's phase-1 result
            # exactly as it entered CrossAlignConfig. Found by directory pattern
            # (not a hard-coded step number), since the channel's own step
            # numbering is an implementation detail of chan_cfg() above.
            b_dir = joinpath(outdir, "B")
            fc_dirs = filter(d -> occursin(r"^\d+_frameconnect$", d), readdir(b_dir))
            @test length(fc_dirs) == 1
            pre = load_smld(joinpath(b_dir, only(fc_dirs), "smld_combined.h5"))
            mean_x(s) = sum(e.x for e in s.emitters) / length(s.emitters)
            # dx_true = 0.1 μm was removed by alignment; require most of it back.
            @test abs(mean_x(result[:B].smld) - mean_x(pre)) >= 0.05
        end
    end
else
    @info "Skipping thorough tests; set SMLM_TEST_FULL=1 to enable"
end
