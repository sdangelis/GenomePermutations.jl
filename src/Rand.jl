"""md
AbstractGenomeDist is a base type for all AbstractGenomeDist objects.

All AbstractGenomeDist objects should have:  

AT LEAST the following properties:  
  
    - genome::String the name of the genome  
    - _regions::IntervalCollection the regions the genome distribution is  
    - _distribution::Distribution the distribution to draw samples from  
    - overlaps::Bool whether the sampled regions can overlap  
    - _seqs the sequences of the regions used for the randomisation  
    - _chr_distribution::Distribution the distribution to draw chromosomes from  
    - _length_distribution::Distribution the distribution to draw lengths from  
    - on_fail::Symbol. The action to take if the randomisation fails to draw a valid interval
    shared methods recognise throw (throws an erorr), :continue (skips the interval) or :orig (returns the original interval)  
  
note properties starting with _ are PRIVATE as they are derived from the regions used, 
please DO NOT ALLOW USERS TO SET THEM DIRECTLY.  
if you want to allow users to edit any of them you need to allow them to genereate 
all of them in one go with something like setfield!(self, :regions, value).   

You must Implement the follwing methods:  
  
    - an appropriate constructor. should enforce on_fail to be one of the standard options, unless specialised methods to handle them are implemented  
    - rand(): a method that returns a random interval from the distribution  
    - randomise(): a method that returns a random interval from the distribution on the same chromosome
    and with the same length of the provided interval  

"""
abstract type AbstractGenomeDist end


# Abstract methods for AbstractGenomeDist

"""
    Base.rand(distribution::AbstractGenomeDist, n::Int)

Wraps the rand(<:AbstractGenomeDit) function to return a collection with nothing as metadata of n random intervals

Returns an interval collection with n intervals drawn from the distribution.

# Note
if n = 1, a collection of size 1 is returned, if n = 0 a 
"""
function Base.rand(distribution::AbstractGenomeDist, n::Int)
    collection = GenomicFeatures.IntervalCollection{Nothing}()
    for i in 1:n
        if distribution.overlaps
            push!(collection, rand(distribution))
        else
            push!(collection, rand(distribution; collection = collection))
        end
    end
    return collection
end 

"""
    Base.show(distribution::AbstractGenomeDist)

shows the general proprieties of any GenomeDistribution
"""
function Base.show(distribution::AbstractGenomeDist)
    """
    A Genome Distribution of type $(typeof(distribution)),
    Generated from $(length(distribution._regions)) intervals from,
    across $(length(distribution._seqs)) sequences.
    Overlaps are $(distribution.overlaps ? "allowed" : "disallowed")
    Fail behaviour: $(distribution.on_fail)

    Distributions:
    - sequence distribution: $(distribution._chr_distribution)
    - lenght distribution: $(distribution._length_distribution)
    - distribution: $(distribution._distribution)
    """
end


"""
    randomise(collection::IntervalCollection, distribution::AbstractGenomeDist) 

Waps randomise(<:AbstractGenomeDist) to randomise a given IntervalCollection
"""
function randomise(collection::GenomicFeatures.IntervalCollection{T}, distribution::AbstractGenomeDist) where {T}
    new_collection = GenomicFeatures.IntervalCollection{T}()
    for interval in collection
        push!(new_collection, randomise(interval, distribution; collection = new_collection))
    end
    return new_collection
end 


"""
    StartMixture(genome, regions, overlaps = false)

generates a StartsMixture

Under the hood it creates a ? 

#arguments

#Avaibale fields
    - genome::String 
"""
mutable struct StartMixture <: AbstractGenomeDist 
    genome::String
    regions::GenomicFeatures.IntervalCollection
    _distribution::Dict{String, Distributions.MixtureModel} 
    overlaps::Bool
    _seqs
    _len_distribution::Dict{String, Distributions.Binomial}
    _chr_distribution::Distributions.Categorical
    by_chromosome::Bool
    on_fail::Symbol
    max_tries::Int

    # inner constructor - required to force distribution to be generated from regions
    function StartMixture(genome, regions, overlaps = false; by_chromosome = false, on_fail = :throw, max_tries = 100)
        
        # Check on_fail is valid 
        if !(on_fail in [:throw, :continue, :orig])
            error("on_fail must be one of :throw, :continue, :orig")
        end

        if on_fail == :orig && overlaps == false
            @warn "returning the original interval on fail may return overlapping segments"
        end
	    #generate start distrubutions

        # initialise objects to push into
	    distribution = Dict{String, Distributions.MixtureModel}()
        len_distribution = Dict{String, Distributions.Binomial}()
        lengths = Dict{String, Int}()
        seqnames = Vector{String}()
        # generate a mixture model for each region
        for seqname in keys(regions.trees)
            # generate a mixture model for each region
            push!(seqnames, seqname) # add sequence to sequence vector
            d = Vector{Distributions.DiscreteUniform}(undef, 0)
            l = Vector{Int}(undef,0)
            for interval in regions.trees[seqname]
                push!(d, Distributions.DiscreteUniform(interval.first, interval.last))
                push!(l, interval.last-interval.first)
            end
            lengths[seqname] = sum(l)
            distribution[seqname] = Distributions.MixtureModel(d, Distributions.Categorical(l/sum(l)))
            len_distribution[seqname] = Distributions.Binomial(round(Int, Statistics.median(l)))
            # generate a mixture model for each region
        end
        # generate a categorical distribution for the chromosomes
        chr_distribution = Distributions.Categorical(collect(values(lengths))/sum(values(lengths)))
        new(genome, regions, distribution, overlaps, seqnames, len_distribution, chr_distribution, by_chromosome, on_fail, max_tries)
    end
end 


function Base.rand(distribution::AbstractGenomeDist;
     collection = GenomicFeatures.IntervalCollection{Nothing}()
    )
    
    for i in 1:distribution.max_tries
        # draw a chromosome
        chr = distribution._seqs[rand(distribution._chr_distribution)]
        # draw a length
        len = rand(distribution._len_distribution[chr])
        # draw a start position
        start = rand(distribution._distribution[chr])
        interval = GenomicFeatures.Interval(chr, start, start+len, GenomicFeatures.Strand('.'), nothing)
        # check if we belong to the regions, else we need to draw again
        isin(interval, distribution.regions) || continue
        # if we allow overlaps we can return straight away
        distribution.overlaps && return interval 
        # else we need to make sure we are not are overlapping the collection before we return
        anyoverlapping(interval, collection) || return interval
    end
    throw(ErrorException("failed to draw a valid interval after $(distribution.max_tries) tries"))
end


function randomise(interval::GenomicFeatures.Interval, distribution::StartMixture;
        collection = GenomicFeatures.IntervalCollection{Nothing}()
    ) 

    for i in 1:distribution.max_tries
        # if by chromosome we don't need to draw a new chromosome
        if distribution.by_chromosome 
            chr = interval.seqname 
        else
            chr = distribution._seqs[rand(distribution._chr_distribution)]
        end
        # use the same length as the interval
        len = interval.last - interval.first
        # draw a start position
        start = rand(distribution._distribution[chr])
        new_interval = GenomicFeatures.Interval(chr, start, start+len, interval.strand, interval.metadata)
        # check if we belong to the regions, else we need to draw again
        isin(new_interval, distribution.regions) || continue
        # if we allow overlaps we can return straight away
        distribution.overlaps && return new_interval 
        # else we need to make sure we are not are overlapping the collection before we return
        anyoverlapping(new_interval, collection) || return new_interval
    end
    if distribution.on_fail == :orig
        return(interval)
    else
        throw(ErrorException("failed to radomise $interval after $(distribution.max_tries) tries"))
    end
end


"""
To Do
The idea here is to rotate all segments by a x amount as if the chromosome was circular.
"""
struct CircularRandomiser <: AbstractGenomeDist 
    # To Do:
end


"""
To Do

the idea here is to IDK what the idea is.
"""
struct RegionRandomiser <: AbstractGenomeDist 
    # To Do:
end


"""
To Do

the idea here is to regenerate the possible distribution at each iteartion
"""
struct AdaptiveRegions <: AbstractGenomeDist 
    # To Do:
end


"""
To Do

the idea here is to create a series of possible start position for each length in a given collection
"""
struct lengthsMixture <: AbstractGenomeDist 
    # To Do:
end