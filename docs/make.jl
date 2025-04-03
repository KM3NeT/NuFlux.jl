using Documenter, NuFlux

makedocs(;
    modules = [NuFlux],
    authors = "Johannes Schumann, Santiago Pena Martinez",
    sitename = "NuFlux.jl",
    format = Documenter.HTML(;
        assets = ["assets/custom.css"],
        sidebar_sitename = false,
        collapselevel = 4,
        warn_outdated = true,
    ),
    warnonly = [:missing_docs],
    checkdocs = :exports,
    pages = [
        "Introduction" => "index.md",
        "API" => "api.md",
    ],
    repo = Documenter.Remotes.URL(
        "https://git.km3net.de/simulation/NuFlux.jl/blob/{commit}{path}#L{line}",
        "https://git.km3net.de/simulation/NuFlux.jl"
    ),
)

deploydocs(;
  repo = "git.km3net.de/simulation/NuFlux.jl",
  devbranch = "main",
  push_preview=true
)
