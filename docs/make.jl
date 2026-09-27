using SMLMAnalysis
using Documenter

DocMeta.setdocmeta!(SMLMAnalysis, :DocTestSetup, :(using SMLMAnalysis); recursive = true)

# Non-exported SMLMAnalysis-owned result/info structs and extension hooks that are
# documented explicitly (api/index.md, api/steps.md, api/io.md) via `SMLMAnalysis.Name`
# `@docs` blocks. Excluded from api/internals.md's `@autodocs Public=false` block below
# so Documenter doesn't see the same docstring twice ("duplicate docs" error).
const DOCUMENTED_ELSEWHERE = Set(
    Any[
        SMLMAnalysis.AnalysisResult, SMLMAnalysis.AnalysisInfo, SMLMAnalysis.StepInfo,
        SMLMAnalysis.DetectFitInfo, SMLMAnalysis.FilterInfo, SMLMAnalysis.IntensityFilterInfo,
        SMLMAnalysis.DensityFilterInfo, SMLMAnalysis.BaGoLInfo,
        SMLMAnalysis.MultiTargetResult, SMLMAnalysis.MultiTargetInfo,
        SMLMAnalysis.AbstractMultiTargetStep,
        SMLMAnalysis.CompositeRenderInfo, SMLMAnalysis.CrossAlignInfo, SMLMAnalysis.CrossCorrInfo,
        SMLMAnalysis.step_name, SMLMAnalysis.step_outdir,
    ]
)
not_documented_elsewhere(x) = !(x in DOCUMENTED_ELSEWHERE)

makedocs(;
    # Include the upstream packages so @docs of re-exported/aliased symbols
    # (FrameConnectConfig, DriftConfig, RenderConfig, BaGoLConfig, the cluster
    # configs, …) resolve their docstrings, which live in those packages.
    modules = [
        SMLMAnalysis,
        SMLMAnalysis.SMLMFrameConnection,
        SMLMAnalysis.SMLMDriftCorrection,
        SMLMAnalysis.SMLMRender,
        SMLMAnalysis.SMLMBaGoL,
        SMLMAnalysis.SMLMClustering,
    ],
    authors = "klidke@unm.edu",
    sitename = "SMLMAnalysis.jl",
    format = Documenter.HTML(;
        canonical = "https://JuliaSMLM.github.io/SMLMAnalysis.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "tutorial.md",
        "User Guide" => "guide.md",
        "Concepts" => [
            "Overview" => "concepts/index.md",
            "The Pipeline Model" => "concepts/pipeline.md",
            "The JuliaSMLM Ecosystem" => "concepts/ecosystem.md",
            "Data Model & Provenance" => "concepts/data_model.md",
            "The Multi-Target Model" => "concepts/multitarget.md",
        ],
        "Workflows" => [
            "Installation & Setup" => "workflows/install.md",
            "Running a Pipeline" => "workflows/running.md",
            "Multi-Dataset" => "workflows/multidataset.md",
            "Multi-Channel" => "workflows/multichannel.md",
            "I/O & Resume" => "workflows/io_resume.md",
            "Extending the Pipeline" => "workflows/extending.md",
            "AI Assistant Guide" => "workflows/agent_guide.md",
            "Troubleshooting" => "workflows/troubleshooting.md",
        ],
        "Pipeline Steps" => [
            "Overview" => "steps/index.md",
            "Detection & Fitting" => "steps/detectfit.md",
            "Quality Filter" => "steps/filter.md",
            "Intensity Filter" => "steps/intensityfilter.md",
            "Frame Connection" => "steps/frameconnect.md",
            "Drift Correction" => "steps/driftcorrect.md",
            "Density Filter" => "steps/densityfilter.md",
            "Rendering" => "steps/render.md",
            "Bayesian Grouping" => "steps/bagol.md",
            "Clustering" => "steps/clustering.md",
            "Edge Classification" => "steps/edgeclassify.md",
            "Multi-Channel Steps" => [
                "Composite Render" => "steps/composite_render.md",
                "Cross-Alignment" => "steps/cross_align.md",
                "Cross-Correlation" => "steps/crosscorr.md",
            ],
        ],
        "References" => "references.md",
        "API Reference" => [
            "Overview" => "api/index.md",
            "Step Configs & Info" => "api/steps.md",
            "Multi-Target & I/O" => "api/io.md",
            "Internals" => "api/internals.md",
        ],
    ],
    # Including the upstream modules above turns on missing-docstring checking
    # for their whole export surface; :none keeps that from flooding the build
    # (SMLMAnalysis is a facade that re-exports a large API). Revisit if a
    # scoped per-page check is wanted.
    # Upstream packages are in modules= (so re-exported @docs resolve), which also
    # makes Documenter run THEIR jldoctests — those assume their own DocTestSetup and
    # fail here. We have no jldoctests of our own, so disable doctests entirely.
    doctest = false,
    checkdocs = :none,
    warnonly = false,  # build is warning-clean (verified from the main checkout)
)

deploydocs(;
    repo = "github.com/JuliaSMLM/SMLMAnalysis.jl",
    devbranch = "main",
)
