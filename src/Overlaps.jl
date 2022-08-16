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

Linearly count how many intervals in collection a overlap with any interval in collection b.

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

See Also [`GenomePermutations.anyoverlapping`](@ref), [`GenomePermutations.alloverlapping`](@ref), 
[`GenomePermutations.countoverlaps`](@ref).
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
	countoverlaps(a::GenomicFeatures.IntervalCollection{T},b::GenomicFeatures.IntervalCollection{S})

linearly counts how many intervals in collection b overlaps with collection a  

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
 2
 2
```

#Note 
As seen above, countoverlapping(a,b) == countoverlapping(b,a) as we count the total number 
of overlaps for each interval

See Also [`GenomePermutations.anyoverlapping`](@ref), [`GenomePermutations.alloverlapping`](@ref), 
[`GenomePermutations.countoverlaps`](@ref)
"""
function countoverlaps(
	a::GenomicFeatures.IntervalCollection{T},
	b::GenomicFeatures.IntervalCollection{S}
) where {T,S}

	overlaps::Int = 0
	for interval_a in a
		for interval_b in b
			GenomicFeatures.isoverlapping(interval_a, interval_b) && (overlaps += 1)
		end 
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

See Also [`GenomePermutations.isin`](@ref), [`GenomePermutations.anyin`](@ref).
"""
function isin(a::GenomicFeatures.Interval{S},b::GenomicFeatures.Interval{T}) where {S,T}
	return (a.seqname == b.seqname && a.first >= b.first && a.last <= b.last)
end


"""
	isin(a::GenomicFeatures.Interval{T}, b::GenomicFeatures.IntervalCollection{S})

Linearly check if interval a is fully contained in any interval of collection b.
	
```jldoctest
using GenomicFeatures 
a = GenomicFeatures.Interval("chr1", 5, 10)
b = GenomicFeatures.IntervalCollection([
	GenomicFeatures.Interval("chr1", 5, 10),
	GenomicFeatures.Interval("chr1", 25, 35 ),
	GenomicFeatures.Interval("chr1", 40, 65)
	])

isin(a, b) 
		
# output

true
```

```jldoctest
using GenomicFeatures 
a = GenomicFeatures.Interval("chr1", 30, 50)
b = GenomicFeatures.IntervalCollection([
	GenomicFeatures.Interval("chr1", 5, 10),
	GenomicFeatures.Interval("chr1", 25, 35),
	GenomicFeatures.Interval("chr1", 35, 65)
	])

isin(a, b) 
		
# output

false
```

#Note 

`isin(a, b)` strictly checks if a is contained in any single interval of b, 
rather than whether a is contained within any combination of intervals in collection b

See Also [`GenomePermutations.isin`](@ref), [`GenomePermutations.anyin`](@ref).
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
	
Linearly check if *all* intervals in collection a are contained in collection b. 

```jldoctest

	using GenomicFeatures 
	a = GenomicFeatures.IntervalCollection([
		GenomicFeatures.Interval("chr1", 6, 9),
		GenomicFeatures.Interval("chr1", 45, 50)
		])
	b = GenomicFeatures.IntervalCollection([
		GenomicFeatures.Interval("chr1", 5, 10),
		GenomicFeatures.Interval("chr1", 25, 35 ),
		GenomicFeatures.Interval("chr1", 40, 65)
		])
	
	[isin(a, b), isin(b, a)]
	
	# output
	
	2-element Vector{Bool}:
	 true
	 false
	```
	
# Note

At this point there are 2 layers of linear search so the time complexity is n^2.

See Also [`GenomePermutations.isin`](@ref), [`GenomePermutations.anyin`](@ref).
"""
function isin(a::GenomicFeatures.IntervalCollection{T},
				 b::GenomicFeatures.IntervalCollection{S}) where {T,S}
	
	for interval in a
		isin(interval, b) || return(false::Bool)
	end
	return(true::Bool)
end

"""
	anyin(a:::GenomicFeatures.IntervalCollection{T}, b:::GenomicFeatures.IntervalCollection{S})
	
Linearly check if *any* interval in collection a are contained in collection b. 

```jldoctest

using GenomicFeatures 
a = GenomicFeatures.IntervalCollection([
	GenomicFeatures.Interval("chr1", 25, 27),
	GenomicFeatures.Interval("chr1", 100, 150)
	])
b = GenomicFeatures.IntervalCollection([
	GenomicFeatures.Interval("chr1", 5, 10),
	GenomicFeatures.Interval("chr1", 25, 35 ),
	GenomicFeatures.Interval("chr1", 40, 65)
	])
	
anyin(a, b)
	
# output

true
```

# Note

At this point there are 2 layers of linear search so the time complexity is n^2

See Also [`GenomePermutations.isin`](@ref), [`GenomePermutations.anyin`](@ref).
"""
function anyin(a::GenomicFeatures.IntervalCollection{T},
				 b::GenomicFeatures.IntervalCollection{S}) where {T,S}
	
	for interval in a
		isin(interval, b) && return(true::Bool)
	end
	return(false::Bool)
end
