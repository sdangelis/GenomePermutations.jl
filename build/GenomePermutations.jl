module GenomePermutations

import BED
import Distributions
import GenomicFeatures
import Statistics
import HypothesisTests

export anyoverlapping, alloverlapping, isin, countoverlapping
export dist
export simpleP
export generatedistribution, randomiseregions, randomisegenome
export SimplePTest, PermTestResult
export permtest, overlappermtest


"""
	anyoverlapping(a::GenomicFeatures.Interval{T}, b::GenomicFeatures.IntervalCollection{S})

Extend GenomicFeatures.isoverlapping to linearly check if interval a.
overlaps collection. Return true or false

```jldoctest
using GenomicFeatures 
a = GenomicFeatures.Interval("chr1", 0, 10)
b = GenomicFeatures.IntervalCollection([
	GenomicFeatures.Interval("chr1", 5, 15),
	GenomicFeatures.Interval("chr1", 20, 25)
	])
anyoverlapping(a, b) 

# output

true
```

See Also ['alloverlapping'](alloverlapping), ['countoverlapping'](countoverlapping).
"""
function anyoverlapping(a::GenomicFeatures.Interval{T}, b::GenomicFeatures.IntervalCollection{S}) where {T, S}
	for interval in b
		GenomicFeatures.isoverlapping(a, interval) && return(true::Bool)
	end
	return(false::Bool)
end


"""
	alloverlapping(a::GenomicFeatures.IntervalCollection{T}, b::GenomicFeatures.IntervalCollection{S})

Linearly check if all intervals in collection a overlap collection b. 
Return true or false
	
	```jldoctest
	using GenomicFeatures 
	a = 	b = GenomicFeatures.IntervalCollection([
		GenomicFeatures.Interval("chr1", 14, 17),
		GenomicFeatures.Interval("chr1", 20, 25)
		])
	b = GenomicFeatures.IntervalCollection([
		GenomicFeatures.Interval("chr1", 5, 10),
		GenomicFeatures.Interval("chr1", 20, 25)
		])
	alloverlapping(a, b) 
	
	# output
	
	false
	```
"""
function alloverlapping(a::GenomicFeatures.IntervalCollection{T}, b::GenomicFeatures.IntervalCollection{S}) where {T,S}
	for interval in a
		anyoverlapping(interval, b) || return(false::Bool)
	end
	return(true::Bool)
end


"""
	vec_getcollection(collection::GenomicFeatures.IntervalCollection{T}, sequence::String) 

(inefficiently) retun a interval collection with only features in the specified sequence 

see also iter_getcollection
"""
function vec_getcollection(collection, sequence::String) 
	GenomicFeatures.IntervalCollection([GenomicFeatures.Interval(i.seqname, 
	i.first, i.last, i.strand, i.metadata ) for i in collection.trees[sequence]])
end


"""
	iter_getcollection(collection::GenomicFeatures.IntervalCollection{T}, sequence::String) 

(inefficiently) retun a interval collection with only features in the specified sequence.

see also vec_getcollection
"""
function iter_getcollection(collection::GenomicFeatures.IntervalCollection{T}, 
							sequence::String) where {T}
	tmp = GenomicFeatures.IntervalCollection{T}()
	for i in collection.trees[sequence]
		push!(tmp, i)
	end
	return(tmp)
end


"""
	countoverlapping(a::GenomicFeatures.IntervalCollection{T} , b::GenomicFeatures.IntervalCollection{T})

Linearly count how many intervals in collection a overlap with collection b.

```jldoctest

using GenomicFeatures 
a = GenomicFeatures.IntervalCollection([
	GenomicFeatures.Interval("chr1", 14, 17),
	GenomicFeatures.Interval("chr1", 20, 50)
	])
b = GenomicFeatures.IntervalCollection([
	GenomicFeatures.Interval("chr1", 5, 10),
	GenomicFeatures.Interval("chr1", 25, 35 ),
	GenomicFeatures.Interval("chr1", 40, 65)
	])

[countoverlapping(a, b), countoverlapping(b, a)]

# output

[1,2]
```

#Note 
As seen above, countoverlapping(a,b) != countoverlapping(b,a) as one interval in a 
can overlap with >1 intervals in b.
Should it be defined this way or be defined as total number of segments in a or b 
that overlap each other? 
"""
function countoverlapping(a::GenomicFeatures.IntervalCollection{T} ,
						  b::GenomicFeatures.IntervalCollection{S}) where {T,S}
	overlaps::Int = 0
	for interval in a
		anyoverlapping(interval, b) && (overlaps += 1) 
	end
	return(overlaps)
end


"""
	isin(a::GenomicFeatures.Interval{S}, GenomicFeatures.Interval{T})

Check if interval a is fully contained in interval b.
"""
function isin(a::GenomicFeatures.Interval{S},b::GenomicFeatures.Interval{T}) where {S,T}
	return (a.seqname == b.seqname && a.first >= b.first && a.last <= b.last)
end


"""
	isin(a,b)

Linearly check if interval a is fully contained in collection b.
"""
function isin(a::GenomicFeatures.Interval{S},
				 b::GenomicFeatures.IntervalCollection{T}) where {S,T}
	
	for interval in b
		isin(a, interval) && return(true::Bool)
	end
	return(false::Bool)
end


"""
	isin(a,b)
	
Linearly check if all intervals in collection a are contained in collection b. 
"""
function isin(a::GenomicFeatures.IntervalCollection{S},
				 b::GenomicFeatures.IntervalCollection{T}) where {S,T}
	
	for interval in a
		isin(interval, b) || return(false::Bool)
	end
	return(true::Bool)
end


"""
	dist(a, b)

Return the minimum unsigend distance between Intervals a and b.
returns missing if the two intervals do not have the same seqname. 
"""
function dist(a::GenomicFeatures.Interval{S},b::GenomicFeatures.Interval{T}) where {S, T}
	
	a.seqname != b.seqname && return missing
	return minimum(abs.([a.first-b.first, a.first-b.last, a.last-b.first, a.last-b.last]))
end


"""
	dist(a,b, f = minimum)

Linearly calculate all the distances between an interval a and a collection b,
return a value according to function f.

# arguments
	- a::GenomicFeatures.interval{T}
	- b::GenomicFeatures.IntervalCollection{S}
	- f::Function summarising function, Defaults to minimum

Summarise all the pairwise distance for a and each interval in b according to the supplied function f
Currently tested with `minimum`, `maximum`, `mean`, and `median`.
"""
function dist(a::GenomicFeatures.Interval{T}, b::GenomicFeatures.IntervalCollection{S},
			  f::Function=minimum) where {S,T}
	
	d = Vector{Union{Int, Float64}}(undef,0)
	for interval in b
		tmp_d = dist(a,interval)
		if ismissing(tmp_d)
			continue
		end
		push!(d, tmp_d)
	end
	# Workaround to deal with empty collecion
	isempty(d) && return(d)
	return f(d)
end


"""
	pass(x)

Return  x, avoid the need for the anon x -> return x
"""
pass(x::T) where {T} = x::T


"""
	dist(a, b, f, g)

Lineraly caluclates the distance between collection a and b.
Each distance between a and b is summarised according to f and then summarised according to g

# agruments
-	a::GenomicFeatures.IntervalCollection{T}
-	b::GenomicFeatures.IntervalCollection{S}
-	f = pass
-	g = minimum

Calculates the distance between each intervals in a and collection b according to function f, defaults to `minimum`.
The overal distance between all elements of a and b (i.e the array of all distance produced by f
is summarised according to the supplied function g. defaults to `pass` to return the full array. 
Currently g is also tested with, `minimum`, `maximum`, `mean`, `median`.
Currently f is only tested explicltly with `minimum`. `maximum`, `mean`, and `median` should work
"""
function dist(a::GenomicFeatures.IntervalCollection{T}, 
			  b::GenomicFeatures.IntervalCollection{S}, 
			  g::Function=GenomePermutations.pass, 
			  f::Function=minimum) where {S,T}
	
	d = Vector{Union{Integer, Float64}}(undef,0)
	for interval in a
		tmp_d  =  dist(interval,b, f)
		# deal with ampty collections
		isempty(tmp_d) && continue 
		push!(d,tmp_d)
	end
	return g(d)
end


"""
	generatedistribution(collection)

Generate a Distributions.MixtureModel from the collection 

The model is composed of a Mixture model of individual discrete uniform distributions
categorical 

Debating if this should be private?

# Note

At the moment it only supports the 
Generating a distribution from overlapping is undocumented
(I assume it would add up the probability but this is untested)
"""
function generatedistribution(collection)
	if length(collecion.trees > 1)
		throw(error("more than 1 sequence detected, Please use iter_getcollection to subset the	collection"))
	end	
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

	PRIVATE - DO NOT USE!

Generate a new random intercval 
Optionally prevents intervals overlapping collection col 

Under the hood?

# Note 

This is a private function. Nothing is modified in place but there is no mechanism 
to ensure the distribution, the regions and, interval belong to the same thing 

generates a new random interval from a given distribution of possible start locations 
ensuring it is contained in the given regions
Optionally intervals in the collection can be disallowed with allow_overlap = false, default: true. 
At this point this function is intened to be PRIVATE as it does not deal with collection of different sequences.
"""
function _randominterval(interval::GenomicFeatures.Interval{T}, 
						 distribution::Distributions.Sampleable, regions::GenomicFeatures.IntervalCollection{S};
						 collection::GenomicFeatures.IntervalCollection{T} = GenomicFeatures.IntervalCollection{T}(), 
						 allow_overlap::Bool = true, max_tries::Integer = 1000) where {S,T}

	il = interval.last - interval.first
	new = GenomicFeatures.Interval{T}
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
	error("""
		  Error: can't randomise regions with $max_tries tries, 
		  Please consider increasing max_iter or use randomisation strategy 2,
		  Please also check the intervals come from the same sequnce used to generate
		  the distribution from the regions
		  """)	
end


"""
	randomiseregions(collecion, distribution, regions, )

	This again should be private and undocumented 
Generates randomised regions from an interval collection according to a start postion  distribution. It guarantees that the new intervals are contained in the given regions.
Note this functions assumes that the interval collection contains intervals belonging to a single sequence that is the same sequence used for the distribution. 
The distribution should be derived from the regions, else the randomiser will (safely) fail
New intervals overlapping with other intervals in the nascient collection are disallowed by default. To allow overlaps use `allow_overlap=true`
"""
function randomiseregions(collection::GenomicFeatures.IntervalCollection{T}, 
						  distribution::Distributions.Sampleable, regions; allow_overlap = false, 
						  max_tries::Integer = 1000) where {T}

	new_collection = GenomicFeatures.IntervalCollection{T}()
	for i in collection
		push!(new_collection, GenomePermutations._randominterval(i, distribution, regions;
		collection = new_collection, allow_overlap, max_tries))
	end 
	return new_collection
end 


"""
	randomisegenome(col, regions; allow_overlap, max_tries)

# arguments

-
-
-



Randomises a collection, chromosome by chromosome along the specified regions, returning the new collection
Note: it will skip any sequences in the collection not found in the target regions. if no sequences are found returns an empty collection.
"""
function randomisegenome(col::GenomicFeatures.IntervalCollection{T}, 
						 regions::GenomicFeatures.IntervalCollection{S}; 
						 allow_overlap=false, max_tries::Integer = 1000) where {S,T}
	
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
			for interval in randomiseregions(a[seq],D[seq],r[seq]; allow_overlap, max_tries)
				push!(ran, interval)
			end
		end
	end
	return ran::GenomicFeatures.IntervalCollection{T}
end


"""
A structure to save RegioneR - style simple P tests
"""
struct SimplePTest
	
	iterations::Integer
	p_val::Float64
	alternative::String
	min_p_val::Float64
	function SimplePTest(iterations, p_val, alternative)
		new(iterations, p_val, alternative, (1/(iterations+1)))
	end
end
 

"""
Calculates the p value for the alternative hypothesis that the observed value is either more or less than what we expect from the Vecor of random values.
If no alternative is provided, an alternative will be automatically selected 
on the basis of the data. 
Please see regioneR for full details
"""
function simpleP(obs, ran, iterations; alternative = "Auto")
	
	if alternative == "Less"
		p = (length(filter(x-> x <= obs, ran)) + 
		1)/(length(ran) + 1)
	elseif alternative == "More" 
		p = ((length(filter(x-> x >= obs, ran)) + 
		1)/(length(ran) + 1))
	else
		# the median should be the same as mean for large numbers but is more robust for small numbers
		obs <= Statistics.median(ran) && return(simpleP(obs, ran, iterations; alternative = "Less"))
		return(simpleP(obs, ran, iterations, alternative ="More"))
	end
	return(SimplePTest(iterations, p, alternative))
end 


"""
A "simple" structre and constructor to store results of any permutation test
"""
struct PermTestResult{S <: Union{Real, Vector{<:Real}}, N <: Union{Real, Vector{<:Real}}, T <: Any}
	
	iterations::Integer
	obs::S 
	ran::Vector{N}
	test::T
	tested_regions::Union{String, Nothing}
	randomised_regions::Union{String, Nothing}
	eval_function::Function
	random_strategy::String
end 


"""
TBD (to be documented)
"""
function overlappermtest(col1, col2, regions, iterations; bed_names = (nothing,nothing), 
					     allow_overlap=false, max_tries::Integer = 1000)
		 
	# count base overlaps 
	obs = countoverlapping(col1,col2)
	if !(isin(col2,regions))  
		@warn "Intervals to randomise are not a subset of target regions for the randomisation"
	end
	if !(keys(col1.trees) ⊆ keys(regions.trees) && keys(col2.trees) ⊆ keys(regions.trees))
		throw(ErrorException("Sequences in either collecion does not match selected regions"))
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

	for i in 1:iterations
		tmprand = Array{Int}(undef,0)
		for seq in keys(col2.trees)
			# check if the sequence is in collection 1 too
			seq in keys(col1.trees) || continue 

			# run chromosome by chromosome.
			rand_b = randomiseregions(b[seq],D[seq],r[seq]; allow_overlap = allow_overlap,  
									  max_tries = max_tries)
			push!(tmprand, countoverlapping(a[seq], rand_b))
		end 
		ran[i] = sum(tmprand)
	end 
	# get p value for overlap test

	test = simpleP(obs, ran, iterations)
		
	return PermTestResult(iterations, obs, ran, test, bed_names[1], bed_names[2], 
						  countoverlapping, "guess and check")
end


"""
dist(a, b, f, g)

# agruments
-	a
-	b
-	f = pass
-	g = minimum

A function to generate n random genome intervals in 


	This is actually safe enough not to need being private.
	What's the point of a private function we never call anyways? 
"""
function _randomgenome(regions, n::Integer, C::Distributions.Categorical, 
	L::Distributions.DiscreteDistribution, names::Vector{String}; allow_overlap = false)

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
		push!(col, _randominterval(int, S[int.seqname], regions; collection = col, allow_overlap = allow_overlap))
	end 
	return(col)
end


"""
Runs a flexible permuation test, using a custsom evaluation function f and,
in the future, an optional custom summary function g
"""
function permtest(col1, col2, regions, iterations, f::Function=GenomePermutations.dist;
	bed_names = (nothing,nothing), allow_overlap=false, max_tries::Integer= 1000)

	# count base overlaps 
	obs = f(col1, col2)
	
	# should spin this off into a function, rather than copy paste from the other function
	if !(isin(col2, regions))  
		@warn "Intervals to randomise are not a subset of target regions for the randomisation"
	end
	if !(keys(col1.trees) ⊆ keys(regions.trees) && keys(col2.trees) ⊆ keys(regions.trees))
		throw(ErrorException("Sequences in either collecion does not match selected regions"))
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
	for i in 1:iterations
		
		tmprand = Array{typeof(obs)}(undef,0)
		for seq in keys(col2.trees)
			# check if the sequence is in collection 1 too
			seq in keys(col1.trees) || continue 

			# run chromosome by chromosome.
			rand_b = randomiseregions(b[seq],D[seq],r[seq]; allow_overlap = allow_overlap,
									  max_tries = max_tries)
			push!(tmprand, f(rand_b, a[seq]))
		end 
		ran[i] = reduce(vcat, tmprand)
	end 


	# get p value for overlap test, if we only have 1 observaton per iteration
	if length(obs) == 1
		test = simpleP(obs, ran, iterations)
	else
		ran = map(Statistics.mean, ran)
		test = HypothesisTests.UnequalVarianceZTest(obs, ran)
	end	
		
	return PermTestResult(iterations, obs, ran, test, bed_names[1], bed_names[2], f,
						  "guess and check") 

end                  

end
