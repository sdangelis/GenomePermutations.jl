"""
	anyoverlapping(a::GenomicFeatures.Interval{T},
	b::GenomicFeatures.IntervalCollection{S})

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

See Also [`GenomePermutations.alloverlapping`](@ref), [`GenomePermutations.countoverlapping`](@ref).
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
a = GenomicFeatures.IntervalCollection([
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

See Also [`GenomePermutations.anyoverlapping`](@ref), [`GenomePermutations.countoverlapping`](@ref).
"""
function alloverlapping(a::GenomicFeatures.IntervalCollection{T}, b::GenomicFeatures.IntervalCollection{S}) where {T,S}
	for interval in a
		anyoverlapping(interval, b) || return(false::Bool)
	end
	return(true::Bool)
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

2-element Vector{Int64}:
 1
 2
```

#Note 
As seen above, countoverlapping(a,b) != countoverlapping(b,a) as one interval in a 
can overlap with >1 intervals in b.

See Also [`GenomePermutations.anyoverlapping`](@ref), [`GenomePermutations.alloverlapping`](@ref).
"""
function countoverlapping(
	a::GenomicFeatures.IntervalCollection{T},
	b::GenomicFeatures.IntervalCollection{S}
) where {T,S}

	overlaps::Int = 0
	for interval in a
		anyoverlapping(interval, b) && (overlaps += 1) 
	end
	return(overlaps)
end


"""
	isin(a::GenomicFeatures.Interval{S}, GenomicFeatures.Interval{T})

Check if interval a is fully contained in interval b.

```jldoctest
using GenomicFeatures 
a = GenomicFeatures.Interval("chr1", 5, 10)
b = GenomicFeatures.Interval("chr1", 1, 15)

isin(a, b) 
	
# output

true
```

```jldoctest
using GenomicFeatures 
a = GenomicFeatures.Interval("chr1", 20, 30)
b = GenomicFeatures.Interval("chr1", 1, 15)

isin(a, b) 
	
# output

false
```

```jldoctest
using GenomicFeatures 
a = GenomicFeatures.Interval("chr1", 10, 20)
b = GenomicFeatures.Interval("chr1", 1, 15)

isin(a, b) 
	
# output

false
```
"""
function isin(a::GenomicFeatures.Interval{S},b::GenomicFeatures.Interval{T}) where {S,T}
	return (a.seqname == b.seqname && a.first >= b.first && a.last <= b.last)
end


"""
	isin(a::GenomicFeatures.Interval{T}, b::GenomicFeatures.IntervalCollection{S})

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
	isin(a:::GenomicFeatures.IntervalCollection{T}, b:::GenomicFeatures.IntervalCollection{S})
	
Linearly check if all intervals in collection a are contained in collection b. 

# Note

At this point there are 2 layers of linear search so the time complexity is n^2.
"""
function isin(a::GenomicFeatures.IntervalCollection{T},
				 b::GenomicFeatures.IntervalCollection{S}) where {T,S}
	
	for interval in a
		isin(interval, b) || return(false::Bool)
	end
	return(true::Bool)
end
