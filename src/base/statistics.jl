############################################################################################
#::. STRUCTS
############################################################################################
abstract type AbstractDistribution end


############################################################################################
#::. FUNCTIONS
# TODO: finish docstrings
############################################################################################
"""
    [1] cdf([S::Type{<:AbstractFloat}], d::AbstractDistribution, [l], u)
"""
cdf(S::Type{<:AbstractFloat}, d::AbstractDistribution, args...) = S(cdf(d, args...))


"""
# ExESS.jl -- `Statistics.mean` Extension

    [1] mean([S::Type{<:AbstractFloat}], d::AbstractDistribution)
"""
Statistics.mean(S::Type{<:AbstractFloat}, d::AbstractDistribution) = S(mean(d))


"""
    [1] mode([S::Type{<:AbstractFloat}], d::AbstractDistribution)
"""
mode(S::Type{<:AbstractFloat}, d::AbstractDistribution) = S(mode(d))


"""
    [1] pdf([S::Type{<:AbstractFloat}], d::AbstractDistribution, x)
"""
pdf(S::Type{<:AbstractFloat}, d::AbstractDistribution, args...) = S(pdf(d, args...))


"""
# ExESS.jl -- `Base.rand` Extension

    [1] rand([S::Type{<:AbstractFloat}], d::AbstractDistribution)
"""
Base.rand(S::Type, d::AbstractDistribution) = S(rand(d))


"""
    [1] rms([S::Type{<:AbstractFloat}], d::AbstractDistribution)
"""
rms(S::Type{<:AbstractFloat}, d::AbstractDistribution) = S(rms(d))


############################################################################################
#::. EXPORTS
############################################################################################
export cdf, pdf, rms