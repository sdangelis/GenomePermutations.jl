"""
dist(a, b)

Return the minimum unsigend distance between Intervals a and b.
returns missing if the two intervals do not have the same seqname. 

```jldoctest

using GenomicFeatures 
a = GenomicFeatures.Interval("chr1", 10, 20)
b = GenomicFeatures.Interval("chr1", 30, 65)
    
[dist(a, b), dist(b, a)]

# output

2-element Vector{Int64}:
 10
 10
```
"""
function dist(a::GenomicFeatures.Interval{S},b::GenomicFeatures.Interval{T}) where {S, T}

    a.seqname != b.seqname && return missing 
    return minimum(abs.([a.first-b.first, a.first-b.last, a.last-b.first, a.last-b.last]))
end


"""
dist(a,b, f = minimum)

Linearly calculate all the distances between an interval a and a collection b,
returning a value summarised according to function f.

# Arguments

- `a::GenomicFeatures.interval{T}`: interval to calculate distances from.
- `b::GenomicFeatures.IntervalCollection{S}`: collection to calculate all distances to.
- `f::Function`: f Summarise all the pairwise distance for a and each interval in b 
according to the supplied function f 
Currently tested with `minimum`, `maximum`, `mean`, and `median`.


```jldoctest
using GenomicFeatures
using Statistics

a = GenomicFeatures.Interval("chr1", 25, 35)

b = GenomicFeatures.IntervalCollection([
	GenomicFeatures.Interval("chr1", 5, 10),
	GenomicFeatures.Interval("chr1", 25, 35 ),
	GenomicFeatures.Interval("chr1", 40, 100)
	])

[dist(a,b, minimum), dist(a,b, maximum), dist(a,b, mean), dist(a,b, median)]
# output

4-element Vector{Float64}:
  0.0
 15.0
  6.666666666666667
  5.0
```
"""
function dist(
    a::GenomicFeatures.Interval{T}, 
    b::GenomicFeatures.IntervalCollection{S},
    f::Function = minimum
) where {T, S}

    d = Vector{Union{Int, Float64}}(undef,0)
    for interval in b
        tmp_d = dist(a,interval)
        if ismissing(tmp_d)
            continue
        end
        push!(d, tmp_d)
    end
    # Workaround to deal with empty collection
    isempty(d) && return(d)
    return f(d)
end


"""
    dist(a, b, f, g)

Lineraly caluclates the distance between collection a and b.
Each distance between a and b is summarised according to f and then summarised 
according to g.

# Agruments

- `a::GenomicFeatures.IntervalCollection{T}`: collection a.
- `b::GenomicFeatures.IntervalCollection{S}`: collectionb b.
- `f::Function = x -> x`: f summarisees the result of the distances 
between each interval of a and all intervals in b.
tested with `x -> x`, `minimum`, 'maximum`, `mean`, `median`
- `g::Function = minimum`: g summarises the distances returned by f for all
items in a. Currently tested with `minimum`. `maximum`, `mean`, and `median`.

# Note

Both g and f can be replaced by user-defined functions. f must take a  
Vector{Union{Int, Float64}} and return Union{Int, Float64}.
g must take a Vector{Union{Int, Float64}} but can return any type.
Due to the fact that the result of f is stored in d before being summarised, it 
is not possible to change f's return type  -
even if we use a function in g that could acept a vector of that type. Indeed
this is the reason for the need to Union{Int, Float64} as the return type of f.
but types & their optimisation in Julia can be funny. 

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

[dist(a,b, x->x, minimum), dist(a,b, x->x, maximum), dist(a, x->x, mean), dist(a,b, x->x,median)]

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