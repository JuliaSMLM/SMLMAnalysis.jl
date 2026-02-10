using SMLMAnalysis
using Documenter

DocMeta.setdocmeta!(SMLMAnalysis, :DocTestSetup, :(using SMLMAnalysis); recursive=true)

makedocs(;
    modules=[SMLMAnalysis],
    authors="klidke@unm.edu",
    sitename="SMLMAnalysis.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaSMLM.github.io/SMLMAnalysis.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Guide" => "guide.md",
        "API Reference" => "api.md",
    ],
    warnonly=[:missing_docs],
)

deploydocs(;
    repo="github.com/JuliaSMLM/SMLMAnalysis.jl",
    devbranch="main",
)
