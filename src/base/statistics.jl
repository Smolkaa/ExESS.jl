############################################################################################
#::. STRUCTS
############################################################################################
abstract type AbstractDistribution end


############################################################################################
#::. FUNCTIONS
# TODO: finish docstrings
############################################################################################
"""
# ExESS.jl -- `Base.rand` Extension
"""
Base.rand(S::Type, d::AbstractDistribution) = S(rand(d))


"""
    [1] cdf([S::Type{<:AbstractFloat}], d::AbstractDistribution, [l], u)
"""
cdf(S::Type{<:AbstractFloat}, d::AbstractDistribution, args...) = S(cdf(d, args...))


"""
    [1] pdf([S::Type{<:AbstractFloat}], d::AbstractDistribution, x)
"""
pdf(S::Type{<:AbstractFloat}, d::AbstractDistribution, args...) = S(pdf(d, args...))


############################################################################################
#::. EXPORTS
############################################################################################
export cdf, pdf