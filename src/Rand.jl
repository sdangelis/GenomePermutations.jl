"""
AbstractGenomeRand is a base type for all GenomeRand objects.

All GenomeRand objects should have:

AT LEAST the following properties:
    - genome::str the name of the genome
    - regions::IntervalCollection the regions the genome distribution is 
    - distribution:: the distribution to draw samples from
    - overlaps::Bool whether the sample can overlap with itself
    - seqs the sequences of the regions used for the randomisation.
    - chr_distribution::Distribution the distribution to draw chromosomes from

Implement the follwing methods:
    - an appropriate constructor
    - rand(): a method that returns a random interval from the distribution
    - randomise(): a method that randomises a given interval 

AbstractGenomeRand has shared methods that can be used by all GenomeRand objects.
    - randomise(::AbstractGenomeRand): a method that wraps randomise() to randomise a given intervalcollection 
    - rand(::AbstractGenomeRand, n::Int): a method that wraps rand() to return a collection of n random intervals
    - pretty printing functions
please override amy abstract-level methods for custome genome rand objects as required
"""
abstract type GenomeRand end


"""
    wrapper around GenomeRand type to create methods that randomise chromosome by chromosome
"""
struct ByChr{T<: GenomeRand}
    # we need a concrete type here
end     


"""
    StartMixture(genome, regions, overlaps = false)

generates a genome rand object.

# note:
    it 
    it 
    it
    it 
"""
struct StartMixture <: GenomeRand 
    genome::String
    regions::GenomicFeatures.IntervalCollection
    distribution::Dict{String, Distributions.MixtureModel} 
    overlaps::Bool
    seqs
    len_distribution::Dict{String, Distributions.Binomial}
    chr_distribution::Distributions.Categorical
    
    # inner constructor - required to force distribution to be generated from regions
    function StartMixture(genome, regions, overlaps = false)
        
	    #generate start distrubutions - Soon to be spun off into a new function
	    distribution = Dict{String, Distributions.MixtureModel}()
        len_distribution = Dict{String, Distributions.Binomial}()
        lengths = Dict{String, Int}()
        # generate a mixture model for each region
        seqs = regions.trees
        for seq in keys(regions.trees)
            # generate a mixture model for each region
            d = Vector{Distributions.DiscreteUniform}(undef, 0)
            l = Vector{Int}(undef,0)
            for interval in regions.trees[seq]
                push!(d, Distributions.DiscreteUniform(interval.first, interval.last))
                push!(l, interval.last-interval.first)
            end
            lengths[seq] = sum(l)
            distribution[seq] = Distributions.MixtureModel(d, Distributions.Categorical(l/sum(l)))
            len_distribution[seq] = Distributions.Binomial(sum(l))
            # generate a mixture model for each region
        end
        # generate a categorical distribution for the chromosomes
        chr_distribution = Distributions.Categorical(collect(values(lengths))/sum(values(lengths)))
        new(genome, regions, distribution, overlaps, seqs, len_distribution, chr_distribution)
    end
end 


function Base.rand(d::GenomeRand) 
    # TO-DO       
    
end


function randomise() 
    # TO-DO
end

"""
"""
struct RandomiseRegions <: GenomeRand 
    # To Do:
end

"""

"""
struct AdaptiveRegions <: GenomeRand 
    # To Do:
end

"""
"""
struct lengthsMixture <: GenomeRand 
    # To Do:
end

"""
"""
struct CircularRegions <: GenomeRand 
    # To Do:
end

"""
"""

const CircularRand = ByChr{CircularRegions} # add a type alias

