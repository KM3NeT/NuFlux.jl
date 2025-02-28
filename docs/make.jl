using Documenter, NuFlux

makedocs(;
    modules = [NuFlux],
    authors = "Johannes Schumann, Santiago Pena Martinez",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages=[
        "Introduction" => "index.md",
        "API" => "api.md",
    ],
    repo="https://github.com/KM3NeT/NuFlux.jl/blob/{commit}{path}#L{line}",
    sitename="NuFlux.jl",
)

deploydocs(;
    repo="github.com/KM3NeT/NuFlux.jl",
    devbranch="main"
)
