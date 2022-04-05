# GenomePermutations

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sdangelis.github.io/GenomePermutations.jl/dev)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sdangelis.github.io/GenomePermutations.jl/dev)
[![Build Status](https://github.com/sdangelis/GenomePermutations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sdangelis/GenomePermutations.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/sdangelis/GenomePermutations.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/sdangelis/GenomePermutations.jl)[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sdangelis.github.io/GenomePermutations.jl/stable)


A package to generate random genomic regions and run permutation tests on genomic ranges in the spirit of RegioneR but implemented in Julia.  

In particualr, it can

- Check if 2 intervals overlap
- Calculate distances between two intervals
- Same as above but with entire lists
- Generate random genomes from a list of allowed regions
- Randomise genomes
- Run a permutation test to compare overlaps and distances (see plot below)  

This is driven by a need for speedier Julia's fast loops and the minimal overheads provided by data structures like GenomicIntervals Intervals and IntervalCollection, compared to Genomic Ranges are enough to significantly improve speed over RegioneR, despite the naivety of the current Julia implementation.

## Disclaimer

Please note this package is at this stage an internal tool and has some known issues:

- not all randomization approaches from RegioneR are implemented.
The current guess and check approach has issues with segments that relatively
long and have only a limited range of possible locations.
- there are no functions for plotting or for processing and filtering bed files. BED.jl
provides basic I/O. I currently rely on R & JuliaCall, with  custom wrapper
R/Julia function to move beds from data frames to interval collections.
- the GenomicFeatures API/Documentation can be very minimal. One day I'll have the courage
to contribute to it or switch to the more high-level GenomicRanges.jl.
- this API is unstable. Function names/exported status can change and there might
be breaking changes.
- The `isin` function is defined where actually base Julia uses `issubset`.
- Linear searches are everywhere whereas interval collections are interval trees
and seem sorted. This would allow some more interesting approaches in searching/
iterating through them.

Any other issues let me know :)
