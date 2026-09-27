using SMLMAnalysis
using SMLMFrameConnection
using SMLMDriftCorrection
using GaussMLE
using Test
using Random
using TOML
using Statistics

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
