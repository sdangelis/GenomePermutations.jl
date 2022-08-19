"""
dist(a, b)

Return the unsigend distance between Intervals a and b.
returns missing if the two intervals do not have the same seqname. 

```jldoctest

using GenomicFeatures 
a = GenomicFeatures.Interval("chr1", 10, 20)
b = GenomicFeatures.Interval("chr1", 30, 65)
    
[dist(a, b), dist(b, a)]

# output

2-element Vector{Int}:
 10
 10
```
"""
function dist(a::GenomicFeatures.Interval{S},b::GenomicFeatures.Interval{T}) where {S, T}

    a.seqname != b.seqname && return missing 
    return minimum(abs.([a.first-b.first, a.first-b.last, a.last-b.first, a.last-b.last]))
end


"""
dist(a,b)

Linearly calculate the minimum distances between an interval a and a collection b,
returning a value summarised according to function f.

# Arguments

- `a::GenomicFeatures.interval{T}`: interval to calculate distances from.
- `b::GenomicFeatures.IntervalCollection{S}`: collection to calculate distance to.

```jldoctest
using GenomicFeatures
using Statistics

a = GenomicFeatures.Interval("chr1", 25, 35)

b = GenomicFeatures.IntervalCollection([
	GenomicFeatures.Interval("chr1", 5, 10),
	GenomicFeatures.Interval("chr1", 25, 35 ),
	GenomicFeatures.Interval("chr1", 40, 100)
	])

dist(a, b)
  0

"""
function dist(
    a::GenomicFeatures.Interval{T}, 
    b::GenomicFeatures.IntervalCollection{S}
) where {T, S}

    d = Vector{Int}(undef,0)
    for interval in b
        tmp_d = dist(a,interval)
        if ismissing(tmp_d)
            continue
        end
        push!(d, tmp_d)
    end
    # Workaround to deal with empty collection
    isempty(d) && return(d)
    return minimum(d)
end


"""
    dist(a, b; f = x->x)

Lineraly caluclates the distance between each item in collection a and collection b,
summarised according to f.

# Agruments

- `a::GenomicFeatures.IntervalCollection{T}`: collection a.
- `b::GenomicFeatures.IntervalCollection{S}`: collectionb b.
- `f::Function = x->x`: g summarises the distances returned by f for all
items in a

# Note

Both g and can be replaced by user-defined functions.
`f` must take a Vector{Int} but can return any type.


```jldoctest
using GenomicFeatures
using Statistics

a = GenomicFeatures.IntervalCollection([
	GenomicFeatures.Interval("chr1", 15, 20),
	GenomicFeatures.Interval("chr1", 150, 200)
	])
b = GenomicFeatures.IntervalCollection([
	GenomicFeatures.Interval("chr1", 5, 10),
	GenomicFeatures.Interval("chr1", 25, 35 ),
	GenomicFeatures.Interval("chr1", 40, 100)
	])

[dist(a,b, x->x), dist(a,b, maximum), dist(a,b, maximum), dist(a, b, median), dist(a,b, x->x, mean)]

# output

4-element Vector{Float64}:
  5.0
 50.0
 27.5
 27.5
```
"""
function dist(
    a::GenomicFeatures.IntervalCollection{T}, 
    b::GenomicFeatures.IntervalCollection{S}, 
    g::Function = x -> x, 
    f::Function = minimum
) where {S,T}

    d = Vector{Union{Int, Float64}}(undef,0)
    for interval in a
        tmp_d  =  dist(interval,b, f)
        # deal with ampty collections
        isempty(tmp_d) && continue 
        push!(d, tmp_d)
    end
    return g(d)
end