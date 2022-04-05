using GenomePermutations
using Documenter

push!(LOAD_PATH,"../src/")
DocMeta.setdocmeta!(GenomePermutations, :DocTestSetup, :(using GenomePermutations); recursive=true)

makedocs(;
    modules=[GenomePermutations],
    authors="Simone De Angelis <simone.deangelis@qmul.ac.uk> <side97.11+jl@gmail.com>",
    repo="https://github.com/sdangelis/GenomePermutations.jl/blob/{commit}{path}#{line}",
    sitename="GenomePermutations.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://sdangelis.github.io/GenomePermutations.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/sdangelis/GenomePermutations.jl",
    devbranch="main",
)
