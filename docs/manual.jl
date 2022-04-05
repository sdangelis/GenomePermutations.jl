using GenomePermutations
using Documenter

cd("docs")
DocMeta.setdocmeta!(
    GenomePermutations, :DocTestSetup, 
    :(using GenomePermutations);
    recursive=true)
# REVISE REGEX
makedocs(sitename="GenomePermutations.jl", modules = [GenomePermutations])