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