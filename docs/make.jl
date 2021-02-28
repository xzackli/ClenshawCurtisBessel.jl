using ClenshawCurtisBessel
using Documenter

DocMeta.setdocmeta!(ClenshawCurtisBessel, :DocTestSetup, :(using ClenshawCurtisBessel); recursive=true)

makedocs(;
    modules=[ClenshawCurtisBessel],
    authors="Zack Li",
    repo="https://github.com/xzackli/ClenshawCurtisBessel.jl/blob/{commit}{path}#{line}",
    sitename="ClenshawCurtisBessel.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://xzackli.github.io/ClenshawCurtisBessel.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/xzackli/ClenshawCurtisBessel.jl",
)
