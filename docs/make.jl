using Documenter
using PileWave

# Setup for doctests in docstrings
DocMeta.setdocmeta!(PileWave, :DocTestSetup, :(using PileWave); recursive=true)

makedocs(;
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        # assets = ["assets/extra_styles.css"],
    ),
    modules = [PileWave],
    sitename = "PileWave.jl",
    pages=[
        "Manual" => "index.md",
    ],
    doctest = true, # :fix
)

deploydocs(
    repo = "github.com/KeitaNakamura/PileWave.jl.git",
    devbranch = "main",
    push_preview = true,
)
