"""
    dist(a, b)

Return the unsigend distance between the two breakpoints of Intervals a and b or 
missing if the two intervals do not have the same seqname. 

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

See Also [`GenomePermutations.featuredist`](@ref)
"""
function dist(a::GenomicFeatures.Interval{S},b::GenomicFeatures.Interval{T}) where {S, T}

    a.seqname != b.seqname && return missing 
    return minimum(abs.([a.first-b.first, a.first-b.last, a.last-b.first, a.last-b.last]))
end

"""
    featuredist(interval feature)

returns the feature distance between the breakpoints of an interval 

Unilke normal distance, the feature distance takes into account overlaps between intervals and figures, 
returning 0 if an interval breakpoint 

See Also [`GenomePermutations.dist`](@ref)
"""
function featuredist(interval::GenomicFeatures.Interval{S}, feature::GenomicFeatures.Interval{T}) where {S, T}
    if GenomicFeatures.isoverlapping(interval, feature)
        return 0
    else 
        return dist(interval, feature)
    end 
end

loopdoc = """
    _loopdist(a, b, f)

_loopdist loops a distance function f that takes 2 intervals over a and b, handling corner cases.

# Aeguments 
- a::GenomicFeatures.Interval or GenomicFeatures.IntervalCollection:  
- b::GenomicFeatures.IntervalCollection
- f::Function

# Note: 
this is a private function that enables `featuredist` and `dist` to share the same loop and
corner case handling core while being flexible enough to plug any function that takes 
two `GenomicFeatures.Interval` objects as arguments and returns an `Int` or `missing`

See Also [`GenomePermutations.dist`](@ref), [`GenomePermutations.featuredist`](@ref)
"""


"$loopdoc"
function _loopdist(a::GenomicFeatures.Interval, b::GenomicFeatures.IntervalCollection, f::Function)
    d = typemax(Int)
    for interval in b
        tmp_d = f(a, interval)
        if ismissing(tmp_d)
            continue
        elseif tmp_d < d
            d = tmp_d
        end
    end
    # Workaround to deal with empty collection
    d == typemax(Int) && return missing 
    return d
end

"$loopdoc"
function GenomePermutations._loopdist(a::GenomicFeatures.IntervalCollection, b::GenomicFeatures.IntervalCollection, f)
    d = Vector{Int}(undef,0)
    for interval in a 
        tmp_d = _loopdist(interval,b, f)
        ismissing(tmp_d) && continue
        push!(d, tmp_d)
    end 
    isempty(d) && return missing
    return(d)
end

function dist(
    a::GenomicFeatures.Interval{T}, 
    b::GenomicFeatures.IntervalCollection{S},
) where {T, S}
    GenomePermutations._loopdist(a, b, dist)
end 

function featuredist(
    interval::GenomicFeatures.Interval{T}, 
    features::GenomicFeatures.IntervalCollection{S}) where {T, S}
    GenomePermutations._loopdist(interval, features, featuredist)
    end

function dist(
    a::GenomicFeatures.IntervalCollection{T}, 
    b::GenomicFeatures.IntervalCollection{S},
) where {T, S}
    GenomePermutations._loopdist(a, b, dist)
end 
    
function featuredist(
    intervals::GenomicFeatures.IntervalCollection{T}, 
    features::GenomicFeatures.IntervalCollection{S}
) where {T, S}
    GenomePermutations._loopdist(intervals, features, featuredist)
end