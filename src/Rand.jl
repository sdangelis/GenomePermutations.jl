"""
"""
abstract type GenomeRand end

"""
"""
struct ByChr{T<: GenomeRand} 
    # we need a concrete type here
end     


"""
"""
struct StartMixture <: GenomeRand 
    # To Do:
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

const CircularRand = ByChr{CircularRegions}
