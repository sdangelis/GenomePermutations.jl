```@meta
CurrentModule = GenomePermutations
```

# GenomePermutations

Documentation for [GenomePermutations](https://github.com/sdangelis/GenomePermutations.jl).

## Index

```@index
```

## Introduction

Often we require to calculate how far or how many overlaps a series of genomic locations
have with a set of predetermined features. We can naively compare our sample with other samples
or a range of randomly selected locations along the genome.
However, the distribution of overlaps or distances depends on the lengths of the intervals

Permutation tests, where we shuffle random regions to a create a random null distribution
can provide a probabilistic framework. By generating a null distribution from intervals of 
the same length and distribution as what we observe, we can examine the location of 
our observation in the distribution to see how likely our value was due to chance alone.

GenomePermutations is an higher-performance - although there is still lots of space for optimization -  
opinionated Julia reimplementation of the permutation approach of [RegioneR](http://bioconductor.org/packages/release/bioc/html/regioneR.html)
that aims to improve ease of use, composability and performance without sacrificing statistical rigor.

GenomePermutations package is designed with composability in mind:
There are 2 main steps in running a permutation test

- instantiate a GenomeDistribution type, subtype of AbstractGenomeDist, which contains both the regions,
   the new distribution and defines methods to draw random samples
- define a distance function - sensible defaults are provided
- run a permutation test.

## Quick Start and Vignette

inspired by R's vignettes, this section will introduce a complete workflow, from loading bed files,
to running the permutation test and visualizing the output

At the moment this is an unregistered package, so we need to provide a GitHub
link to manually install it. 

Note this will install from the main 

```Julia
using Pkg

Pkg.add(url="TBD")
```

if you want the nightly, most unstable and up to date you'll have to install
from the `dev` branch

We then load our BED files with `BED.jl`

```Julia
using GenomePermutations 
using BED
using plots

CNAs =  open(BED.Reader, "CNAs.bed") do reader
    IntervalCollection(reader, true)
end

features =  open(BED.Reader, "features.bed") do reader
    IntervalCollection(reader, true)
end

regions =  open(BED.Reader, "regions.bed") do reader
    IntervalCollection(reader, true)
end
```

Calculating distances:

```Julia
distances = dist(CNAs, features)
```

Generating a genome distributions:

```Julia
dist = GenomePermutations.StartMixture("hgTest", regions, false)
```

Running a permutation test:

```Julia
test = permtest(CNAs,  features, dist, 1000)

# get the P value
pvalue(test.test)

histogram(test.ran, bins = :scott, color_palette=["grey75"]) # Ok, not the prettiest but it can fly for now
```

## Checking if intervals are in a collection

These functions extend functions in GenomicFeatures to check if collections are contained 
in intervals or vice versa.

```@docs
isin
anyin
```

## Overlap calculations

GenomePermutations contains functions to calculate overlaps between intervals and/or collections

```@docs
anyoverlapping
alloverlapping
countoverlapping
countoverlaps
```

## Distance calculations

GenomePermutations contains functions to calculate distances between intervals and/or collections.

```@docs
dist
featuredist
```

>**NOTE:*
> Distance functions uses the maximum integer representable in memory to intialise the distance for
> performance reasons. However, in the unlikely event the real distance is equal to this value
> the function will return `missing`.

## Genome Distributions

All Genome Distributions have shared methods.

All implement two methods, one to draw a random interval, and another to randomize a given interval.

The main difference is while `randomise` will use the length of the provided interval (and if by_chromosome = true, the chromosome) to, rand will always draw them at random from the lengths and chromosome distribution.  
Moreover, while `rand` will always return empty metadata with nothing type for metadata,
`randomise` should conserve type and content of the metadata from the original  - although this depends on the exact implementation of the specific genome distribution type

The following methods that extend rand, and `randomise` along a whole collection are shared by all Genome Distributions, unless otherwise specified.

```@docs
Base.rand(::AbstractGenomeDist, ::Int)
randomise(::GenomicFeatures.IntervalCollection, ::AbstractGenomeDist)
```

### Start Mixture

The core idea behind the Start Mixture genome distribution is that 
we can create a distribution of the all possible locations an interval can start and then check 
at sampling time if an interval of given length l can start at that location and still be contained in the allowed regions.
If so, and optionally after checking it is not overlapping the nascient new interval list if the the interval is valid

```@docs
StartMixture
```

### Work in progress

Currently, there are some Genome Distributions - including some lifted straight out of RegioneR -
that are sketched out but not yet implemented.

```@docs
CircularRandomiser
RegionRandomiser
AdaptiveRegions
lengthsMixture
```

if any sounds like the kind, I'm always looking for contributors.
See the [section](@ref Extending-Genome-Permutations) on extending Genome Permutations for more info

## Permutation tests

Before we discuss the permutation test function, 
We need to introduce the structure that holds the results

```@docs
GenomePermutations.PermTestResult
HypothesisTests.pvalue(::PermTestResult)
```

In the spirit of RegioneR, we can derive a quick P-value by simply

```@docs
GenomePermutations.simpleP
```

```@docs
HypothesisTests.pvalue(::SimplePTest)
GenomePermutations.SimplePTest
```

We can then run a permutation test with `permtest`

```@docs
GenomePermutations.permtest(::GenomicFeatures.IntervalCollection, ::GenomicFeatures.IntervalCollection, ::AbstractGenomeDist, ::Int64)
```

## Extending Genome Permutations

While sensible distance functions and genomeDistributions are provided,
It should be possible to add new custom `GenomDist` types.

These custom types must be a subtype of the abstract `AbstractGenomeDist` type and must conform to the specification below

```@docs
AbstractGenomeDist
```

For reference, the following shared methods are available for all concrete subtypes of `AbstractGenomeDistance`:

```Julia
Base.rand(::AbstractGenomeDist, ::Int)
randomise(::GenomicFeatures.IntervalCollection, ::AbstractGenomeDist)
Base.show(::AbstractGenomeDist)
```

It should also be possible to create a new distance function,
as long as it returns 1 value (an `Int`, or `missing`) per interval we can rely on looping helpers to create appropriate versions that can handle collections, as needed

```@docs
_loopdist
```

You can then alias your distance as needed

```Julia
function mydist(a::Interval, b::IntervalCollection) 
    _loopdist(a, b, mydist)
end 

function  mydist(a::IntervalCollection, b::IntervalCollection)
    _loopdist(a, b, mydist)
end
```

in fact simply calling `_loopdist` from an untyped function should suffice.

```Julia
function mydist(a, b) 
    _loopdist(a, b, mydist)
end
```

As long as you conform  to both specifications, your new distance and genome distribution should work with the permutation test.

## Legacy Functions

These functions were in use prior to the new Genome Distribution type systems are are deprecate due to their unstable and at times confusing output and clunky function calls.  
Will be officially deprecated and/or removed altogether by version 1.0.0

**DEPRECATION WARNING!: PLEASE DO NOT USE THEM**

```@docs
generatedistribution
randominterval
randomisegenome
randomiseregions
overlappermtest
GenomePermutations.permtest(::GenomicFeatures.IntervalCollection, ::GenomicFeatures.IntervalCollection, ::GenomicFeatures.IntervalCollection, ::Int64)
```

These functions are not only deprecated but are also to be considered internal and private

```@docs
GenomePermutations.iter_getcollection
GenomePermutations.vec_getcollection
GenomePermutations._randomgenome
GenomePermutations._randominterval
```
