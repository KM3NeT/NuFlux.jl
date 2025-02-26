using Documenter, NuFlux

makedocs(;
    modules = [NuFlux],
    authors = "Johannes Schumann, Santiago Pena Martinez",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages=[
        "Introduction" => "index.md",
    ],
    repo="https://github.com/8me/NuFlux.jl/blob/{commit}{path}#L{line}",
    sitename="NuFlux.jl",
)

deploydocs(;
    repo="github.com/8me/NuFlux.jl",
    devbranch="main"
)
