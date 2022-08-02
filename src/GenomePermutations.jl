module GenomePermutations

import BED
import Distributions
import GenomicFeatures
import Statistics
import HypothesisTests

include("Overlaps.jl")
include("Dist.jl")

export anyoverlapping, alloverlapping, isin, anyin, countoverlapping
export dist 
export generatedistribution
export randomiseregions, randomisegenome
export simpleP, SimplePTest, pvalue
export permtest, overlappermtest, PermTestResult


"""
	vec_getcollection(collection::GenomicFeatures.IntervalCollection{T}, sequence::String) 

(inefficiently) retun a interval collection with only features in the specified sequence 

See Also [`GenomePermutations.iter_getcollection`](@ref)
"""
function vec_getcollection(collection, sequence::String) 
	GenomicFeatures.IntervalCollection([GenomicFeatures.Interval(i.seqname, 
	i.first, i.last, i.strand, i.metadata ) for i in collection.trees[sequence]])
end


"""
	iter_getcollection(collection::GenomicFeatures.IntervalCollection{T}, sequence::String) 

(inefficiently) retun a interval collection with only features in the specified sequence.

see also [`GenomePermutations.vec_getcollection`](@ref)
"""
function iter_getcollection(
	collection::GenomicFeatures.IntervalCollection{T}, 
	sequence::String
) where {T}
	tmp = GenomicFeatures.IntervalCollection{T}()
	for i in collection.trees[sequence]
		push!(tmp, i)
	end
	return(tmp)
end

"""
	generatedistribution(collection)

Generate a Distributions.MixtureModel from the collection 

The model assumes each position of the genome has equal probability of being sampled.
Practically, it is modelled as a mixture model of individual discrete uniform distributions,
one for each interval in the collection.
The component distributions are distribuited according to a categorical distribution, with 
probabilities proportional to the length of each interval.

# Note

At the moment it only supports collections with single sequence.
Generating a distribution from an overlapping collection is undocumented
(I assume it would add up the probability but this is untested)
"""
function generatedistribution(collection) 

	d = Vector{Distributions.DiscreteUniform}(undef, 0)
	l = Vector{Int}(undef,0)
	for i in collection
		push!(d, Distributions.DiscreteUniform(i.first, i.last))
		push!(l, i.last-i.first)
	end
	l = l/sum(l)
	return Distributions.MixtureModel(d, Distributions.Categorical(l))
end


"""
	_randominterval(interval::GenomicFeatures.Interval{T}, 
		distribution::Distributions.Sampleable, regions::GenomicFeatures.IntervalCollection{S};
		collection::GenomicFeatures.IntervalCollection{T} = GenomicFeatures.IntervalCollection{T}(), 
		allow_overlap::Bool = true, max_tries::Int = 1000) where {T, S}}

This is the private version of randominterval, it skips sequence assertions in the name of performance
Does it actually improve performance? Debeatable, in theory it prevents 100s of IFs. In practice,
I haven't tested it. Indeed, branch prediction on asymmetrical branches should mean 
that the ifs are not really that damaging to performance.
"""
function _randominterval(
	interval::GenomicFeatures.Interval{T}, 
	distribution::Distributions.Sampleable, 
	regions::GenomicFeatures.IntervalCollection{S};
	collection::GenomicFeatures.IntervalCollection{T} = GenomicFeatures.IntervalCollection{T}(), 
	allow_overlap::Bool = true, 
	max_tries::Int = 1000, 
	onfail = :throw
) where {T,S}

	il = interval.last - interval.first
	new = GenomicFeatures.Interval{T}
	
	# Can you turn this into branchless code for performance? 
	for i in 1:max_tries
		newfirst = rand(distribution,1)
		new = GenomicFeatures.Interval(interval.seqname, newfirst[1], newfirst[1]+il,
				interval.strand, interval.metadata)
		if allow_overlap
				isin(new, regions) && return new
		else 
				isin(new,regions) && !(anyoverlapping(new, collection)) && return new
		end
	end 
	# allow passsing of the orignal interval if a random interval cannot be generated.
	if onfail == :orig
		if !(allow_overlap)
			anyoverlapping(interval, collection) && @warn "Returning the original interval can currently cause overlapping segements"
			return interval
		end 
		# if we allow overlapping, we can just return the original interval
		return interval
	elseif onfail == :skip
		return nothing
	else
		error("""
			Error: can't randomise regions with $max_tries tries, 
			Please consider increasing max_iter,
			Please also check the intervals come from the same sequnce used to generate
			the distribution from the regions
		""")
	end	
end


"""
	randominterval(interval, distribution, regions, collection; allow_overlap, max_tries, onfail)

Randomises the position of an interval according to a given distribution. 

# Arguments 

- `interval::GenomicFeatures.IntervalCollection`: The interval to randomise
- `distribution::Distributions.Sampleable`: distribution of possible interval start location.
Should be a discrete distribution, although right now this is not enforced at the signature
level
- `regions::`: the regions the new interval should be fully contained in (checked by`isin`)
- `collecton::GenomicFeatures.IntervalCollection =  GenomicFeatures.IntervalCollection{T}()`: Collection to avoid overlapping 
to if allow_overlap is set to false, defaults to an empty collection of metadata .
- `allow_overlap::Bool = true`: whether to allow overlaps to collection.
- `max_tries::Int = 1000`: the maximum number of attempts to draw an interval that fits in
the region from the distribution
- `onfail::Symbol = :throw`: what to do if the interval fails to fit in the region.
Will return the orignal interval if :orig, returns nothing if :skip or throw an error if :throw.
	
# Note 

To prevent intervals in a loop from overlapping each other, use 
`allow_overlap = false` and add the result of each iteration to the collection
paseed to the function.

See Also (_randominterval)[@ref]
"""
function randominterval(
	interval::GenomicFeatures.Interval{T}, 
	distribution::Distributions.Sampleable, 
	regions::GenomicFeatures.IntervalCollection{S};
	collection::GenomicFeatures.IntervalCollection{T} = GenomicFeatures.IntervalCollection{T}(), 
	allow_overlap::Bool = true, 
	max_tries::Int = 1000, 
	onfail = :throw
) where {S,T}

	# check interval key is in regions
	interval.seqname in keys(regions.trees) || throw(KeyError("$interval.seqname not found in $regions"))	
	# check distribution is valid, given the regions
	if !(distribution == generatedistribution(iter_getcollection(regions, interval.seqname)))
		throw(KeyError("$distribution does not match the regions for $interval.seqname"))
	end 

	return _randominterval(
		interval, distribution, regions; 
		collection = collection, 
		allow_overlap = allow_overlap, 
		max_tries = max_tries,
		onfail = onfail
	)
	
end 

"""
	randomiseregions(collection, distribution, regions; allow_overlap, max_tries, onfail)

Generates randomised regions from an interval collection according to a start postion distribution.
It guarantees that the new intervals are contained in the given regions.

# Arguments

- `collection::GenomicFeatures.IntervalCollection{T}`: collection to randomise,
must have only 1 sequence.
- `distribution::Distributions.Sampleable`: distribution of possible interval
start locations, must match collection.
- `regions::GenomicFeatures.IntervalCollection{S}`: regions to randomise the intervals in.
- `allow_overlap::Bool = false.
- `max_tries::Int = 1000`: maximum number of attempts to draw an interval
that fits in the region from the distribution.
- `onfail::Symbol = :throw`: what to do if an interval fails to fit in the region.
Will return the orignal interval if :orig, skip the interval if :skip, otherwise throws an error.

# Notes

this functions assumes that the interval collection contains intervals belonging 
to a single sequence that is the same sequence used for the distribution. 
The distribution should be derived from the regions.
If the distribution and the regions don't match it will most likely fail safely
New intervals overlapping with other intervals in the nascient collection are disallowed by default.
To allow overlaps use `allow_overlap = true
"""
function randomiseregions(
	collection::GenomicFeatures.IntervalCollection{T}, 
	distribution::Distributions.Sampleable, 
	regions::GenomicFeatures.IntervalCollection{S};
	allow_overlap = false, max_tries::Int = 1000, onfail = :throw
) where {T, S}

	new_collection = GenomicFeatures.IntervalCollection{T}()
	for i in collection
		# we no longer handle errors here. If onfail is set to :skip, we skip 
		# before adding nothing to the collection 
		new_interval =  GenomePermutations._randominterval(i, distribution, regions;
				collection = new_collection, allow_overlap, max_tries, onfail = onfail)
		# Should be a quite asymmetrical branch so not too much performance drag
		isnothing(new_interval) || push!(new_collection, new_interval)
	end 
	return new_collection
end 

"""
	randomisegenome(col, regions; allow_overlap, max_tries)

Randomises a collection, sequence by sequence along the specified regions, returning the new collection

# rguments

- `col::GenomicFeatrues.IntervalCollection{T}`: collection to randomise
- `regions::GenomicFeatures.IntervalCollection{S}`: regions to randomise the collection in.
- `allow_overlap::Bool = false`: whether to allow overlapping intervals in the result.
- `max_tries::Int = 1000`: maximum number of attempts to draw an interval
that fits in the regions.
- `onfail::Symbol = :throw`: what to do if an interval fails to fit in the regions.
Will return the orignal interval if :orig, skip the interval if :skip, otherwise throws an error.

Note: it will skip any sequences in the collection not found in the target regions.
If no sequences are found returns an empty collection.
"""
function randomisegenome(
	col::GenomicFeatures.IntervalCollection{T}, 
	regions::GenomicFeatures.IntervalCollection{S}; 
	allow_overlap = false, max_tries::Int = 1000,
	onfail=:throw
) where {S,T}
	
	# iterate per chromosome, generating all distributions 
	ran = GenomicFeatures.IntervalCollection{T}()
	
	D = Dict{String, Distributions.MixtureModel}()
	for seq in keys(regions.trees)
		push!(D, seq => generatedistribution(regions.trees[seq]))
	end

	a = Dict{String, 
	GenomicFeatures.IntervalCollection{T}}()
	for seq in keys(col.trees)
			push!(a, seq => GenomePermutations.iter_getcollection(col, seq))
	end

	r = Dict{String, 
	GenomicFeatures.IntervalCollection{S}}()
	for seq in keys(regions.trees)
		push!(r, seq => GenomePermutations.iter_getcollection(regions, seq))
	end
		
	for seq in keys(col.trees)

		if !(seq in keys(regions.trees))
			@warn "sequence $seq is not found in the target regions, skipping it..." 
			continue
		else
			# This is ugly but it works 
			for interval in randomiseregions(a[seq],D[seq],r[seq]; allow_overlap, max_tries, onfail)
				push!(ran, interval)
			end
		end
	end
	return ran::GenomicFeatures.IntervalCollection{T}
end


"""
	simpleP(obs, ran, iterations; alternative :Auto)

Calculates the p value for the alternative hypothesis that the observed value
is either more or less than what we expect from the Vecor of random values.

If no altrantive is provided, will be automatically selected 
on the basis of the data.

See also: [regioneR](https://doi.org/doi:10.18129/B9.bioc.regioneR), [`SimplePTest`](@ref)
Please see regioneR for full details
"""
function simpleP(obs, ran, iterations; alternative = :Auto)
	
	if alternative == :Less
		p = (length(filter(x-> x <= obs, ran)) + 
		1)/(length(ran) + 1)
	elseif alternative == :More 
		p = ((length(filter(x-> x >= obs, ran)) + 
		1)/(length(ran) + 1))
	else
		# the median should be the same as mean for large numbers but is more robust for small numbers
		obs <= Statistics.median(ran) && return(simpleP(obs, ran, iterations; alternative = :Less))
		return(simpleP(obs, ran, iterations, alternative = :More))
	end
	return(SimplePTest(iterations, p, alternative))
end 


"""
	SimplePTest(iterations::Int, p_val::, alternative::)

A structure to save RegioneR - style simple P tests. Implements `pvalue`
"""
struct SimplePTest
	
	iterations::Int
	p_val::Float64
	alternative::Symbol
	min_p_val::Float64
	function SimplePTest(iterations, p_val, alternative)
		new(iterations, p_val, alternative, (1/(iterations+1)))
	end
end
 

"""
	PermTestResult{S <: Union{Real, Vector{<:Real}}, N <: Union{Real, Vector{<:Real}}, T <: Any}

A structure to store results of any permutation test
implements pvalue.
"""
struct PermTestResult{S <: Union{Real, Vector{<:Real}}, N <: Union{Real, Vector{<:Real}}, T <: Any}
	
	iterations::Int
	obs::S 
	ran::Vector{N}
	test::T
	tested_regions::Union{String, Nothing}
	randomised_regions::Union{String, Nothing}
	eval_function::Function
	random_strategy::String
end 


"""
	HypothesisTests.pvalue(test::GenomePermutations.PermTestResult)

Extends HypothesisTests.pvalue to get the P value of a GenomePermutations.PermTestResult.
"""
function HypothesisTests.pvalue(test::GenomePermutations.PermTestResult)
	return HypothesisTests.pvalue(test.test)
end


"""
	HypothesisTests.pvalue(test::GenomePermutations.SimplePTest)

Extends HypothesisTests.pvalue to get the P value of a GenomePermutations.SimplePTest.
"""
function HypothesisTests.pvalue(test::GenomePermutations.SimplePTest)
	return(test.p_val)
end


"""
	overlappermtest(col1, col2, regions, iterations; bed_names, allow_overlap = false,
	max_tries = 1000, onfail = :throws)

Legacy function to run an overlap permutation test across 2 regions.
will be replaced by `permtest` using countoverlapping as `f` argument.
Will need a `g` function to summarise overlaps across regions. (likely `sum`).

# Arguments

- `col1::GenomicFeatures.IntervalCollection{T}`: the collection to randomise.
- `col2::GenomicFeatures.IntervalCollection{S}`: the constant collection to test.
- `regions::GenomicFeatures.IntervalCollection{S}`: the regions col1 is randomised in.
- `iterations::Int`: the number of iterations to run.
- `bed_names::Union{String, Nothing}`: the names of the bed files to write to.
- `allow_overlap::Bool = false`: whether to allow overlaps in the randomised regions.
- `max_tries::Int`: the maximum number of tries to generate a randomised collection.
- `onfail::Symbol = :throw`: what to do if a random interval fails to fit in the region.
Will return the orignal interval if :orig, skip the interval if :skip, otherwise throws an error.
"""
function overlappermtest(col1, col2, regions, iterations;
	bed_names = (nothing,nothing),
	allow_overlap = false, 
	max_tries::Int = 1000, 
	onfail = :throw
)
		 
	# count base overlaps 
	obs = countoverlapping(col1,col2)
	if !(isin(col2,regions))  
		@warn "Intervals to randomise are not a subset of target regions for the randomisation"
	end
	if !(keys(col1.trees) ⊆ keys(regions.trees) && keys(col2.trees) ⊆ keys(regions.trees))
		throw(ErrorException("Sequences in either collection does not match selected regions"))
	end
	# iterate per chromosome, generating all distributions 
	D = Dict{String, Distributions.MixtureModel}()
	for seq in keys(regions.trees)
		push!(D, seq => generatedistribution(regions.trees[seq]))
	end
	
	# iterate per chromosome get all regions 
	
	r = Dict{String, GenomicFeatures.IntervalCollection{}}()
	for seq in keys(regions.trees)
		push!(r, seq => GenomePermutations.iter_getcollection(regions, seq))
	end
	
	a = Dict{String, GenomicFeatures.IntervalCollection{}}()
	for seq in keys(col1.trees)
		push!(a, seq => GenomePermutations.iter_getcollection(col1, seq))
	end
	
	b = Dict{String, GenomicFeatures.IntervalCollection{}}()
	for seq in keys(col2.trees)
		push!(b, seq => GenomePermutations.iter_getcollection(col2, seq))
	end
	

	# over the number of iterations
	ran = Array{Int}(undef, iterations)

	# Because we made an array of right size, no need for bounds checks
	Threads.@threads for i in 1:iterations
		tmprand = Array{Int}(undef,0)
		for seq in keys(col2.trees)
			# check if the sequence is in collection 1 too
			seq in keys(col1.trees) || continue 

			# run chromosome by chromosome.
			rand_a = randomiseregions(a[seq],D[seq],r[seq]; allow_overlap,  
									  max_tries, onfail)
			push!(tmprand, countoverlapping(rand_a, b[seq]))
		end 
		ran[i] = sum(tmprand)
	end 
	# get p value for overlap test

	test = simpleP(obs, ran, iterations)
		
	return PermTestResult(iterations, obs, ran, test, bed_names[1], bed_names[2], 
						  countoverlapping, "guess and check")
end


"""	
	_randomgenome(regions, n, C, L, names; allow_overlap = false)

Generates a new random genome with n segments from a named sequence distribution

# Arguments

- `regions`: the regions the genome to generate belongs to.
- `n::Int`: the number of segments to generate.
- `C::Distributions.Categorical`: a distibution with each cateogy representing
one of the named sequences. Lenght and order must match names.
- `L::Distributions.DiscreteUniform": the segment lenght distribution.
- `names::Vector{Strimg}`: the names of the sequences. must match the categories in C.
- `allow_overlap::Bool = false`: whether to allow overlapping segments in the result.

This is currently private as I have not tested this and I don't like the categorical distribution
squence names situation.
"""
function _randomgenome(
	regions,
	n::Int, C::Distributions.Categorical, 
	L::Distributions.DiscreteDistribution, 
	names::Vector{String};
	allow_overlap = false
)

	#generate start distrubutions - Soon to be spun off into a new function
	S = Dict{String, Distributions.MixtureModel}()
	for seq in keys(regions.trees)
		push!(S, seq => generatedistribution(regions.trees[seq]))
	end
	
	# Check sequence distribution
	Distributions.ncategories(C) == length(names) || throw(KeyError("Number of Sequence in distribution $C does not match names $names"))
	length(S) == length(names) ||  throw(KeyError("Number of sequences in $regions does not match names $names"))

	col = GenomicFeatures.IntervalCollection{Nothing}()
	for i in 1:n
		int = GenomicFeatures.Interval(names[rand(C)],1, rand(L)+1)
		push!(col, _randominterval(int, S[int.seqname], regions; collection = col, allow_overlap = allow_overlap, onfail = :throw))
	end 
	return(col)
end


"""
	permtest(col1, col2, regions, iterations, f; allow_overlap = false, max_tries::Int = 1000, onfail = :throw)

Runs a flexible permuation test, using a custsom evaluation function `f`.

# Arguments

- `col1::GenomicFeatures.IntervalCollection{}`: the collection to randomise.
- `col2::GenomicFeatures.IntervalCollection{}`: the constant collection to test.
- `regions::GenomicFeatures.IntervalCollection{}`: the regions col1 is randomised in.
- `iterations::Int`: the number of iterations to run.
- `f::Function = GenomePermutations.dist`: the evaluation function. 
Defaults to GenomePermutations.dist.
- `bed_names::Union{String, Nothing}`: the names of the bed files to write to.
- `allow_overlap::Bool = false`: whether to allow overlaps in the randomised regions.
- `max_tries::Int = 1000`: the maximum number of tries to generate a randomised collection.
- `onfail::Symbol = :throw`: what to do if a random interval fails to fit in the region.
Will return the orignal interval if :orig, skips it if :skip, otherwise throws an error.
"""
function permtest(col1, col2, regions, iterations, f::Function = GenomePermutations.dist;
	bed_names = (nothing,nothing), allow_overlap = false, max_tries::Int = 1000, onfail = :throw)

	# count base overlaps 
	obs = f(col1, col2)
	
	# should spin this off into a function, rather than copy paste from the other function
	if !(isin(col2, regions))  
		@warn "Intervals to randomise are not a subset of target regions for the randomisation"
	end
	if !(keys(col1.trees) ⊆ keys(regions.trees) && keys(col2.trees) ⊆ keys(regions.trees))
		throw(ErrorException("Sequences in either collection does not match selected regions"))
	end

	# iterate per chromosome, generating all distributions 
	D = Dict{String, Distributions.MixtureModel}()
	for seq in keys(regions.trees)
		push!(D, seq => generatedistribution(regions.trees[seq]))
	end
	
	# iterate per chromosome get all regions 
	r = Dict{String, GenomicFeatures.IntervalCollection{}}()
	for seq in keys(regions.trees)
		push!(r, seq => GenomePermutations.iter_getcollection(regions, seq))
	end
	
	a = Dict{String, GenomicFeatures.IntervalCollection{}}()
	for seq in keys(col1.trees)
		push!(a, seq => GenomePermutations.iter_getcollection(col1, seq))
	end
	
	b = Dict{String, GenomicFeatures.IntervalCollection{}}()
	for seq in keys(col2.trees)
		push!(b, seq => GenomePermutations.iter_getcollection(col2, seq))
	end
	
	# over the number of iterations
	
	ran = Array{typeof(obs)}(undef, iterations)
	Threads.@threads for i in 1:iterations
		
		tmprand = Array{typeof(obs)}(undef,0)
		for seq in keys(col2.trees)
			# check if the sequence is in collection 1 too
			seq in keys(col1.trees) || continue 

			# run chromosome by chromosome.
			rand_a = randomiseregions(a[seq],D[seq],r[seq]; allow_overlap,
									  max_tries, onfail)
			push!(tmprand, f(rand_a, b[seq]))
		end 
		ran[i] = reduce(vcat, tmprand)
	end 


	# get p value for overlap test, if we only have 1 observaton per iteration
	if length(obs) == 1
		test = simpleP(obs, ran, iterations)
	else
		ran = map(Statistics.mean, ran)
		obs = Statistics.mean(reduce(vcat, obs))
		#obs = map(Statistics.mean, map(Statistics.mean, obs))
		test = simpleP(obs, ran, iterations)
	end
		
	return PermTestResult(iterations, obs, ran, test, bed_names[1], bed_names[2], f,
						  "guess and check") 

end                  

end
