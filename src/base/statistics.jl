############################################################################################
#::. STRUCTS
############################################################################################
abstract type AbstractDistribution end


############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    [1] cdf([S::Type{<:AbstractFloat}], d::AbstractDistribution, [l], u)

Computes the cumulative distribution function (CDF) of the distribution `d` between the
lower limit `l` and upper limit `u`. If only the upper limit `u` is provided, the lower
limit is assumed to be the minimum value of the distribution. Additionally accepts a type
`S` to cast the result to a specific floating point type.

Note that the format of the limits `l` and `u` depend on the distribution type.
"""
cdf(S::Type{<:AbstractFloat}, d::AbstractDistribution, args...) = S(cdf(d, args...))


"""
    [1] mean([S::Type{<:AbstractFloat}], d::AbstractDistribution)

`Statistics.mean` ExtensionComputes the mean of the distribution `d`. Additionally accepts a 
type `S` to cast the result to a specific floating point type.
"""
Statistics.mean(S::Type{<:AbstractFloat}, d::AbstractDistribution) = S(mean(d))


"""
    [1] mode([S::Type{<:AbstractFloat}], d::AbstractDistribution)

Computes the mode of the distribution `d`. Additionally accepts a type `S` to cast the
result to a specific floating point type.
"""
mode(S::Type{<:AbstractFloat}, d::AbstractDistribution) = S(mode(d))


"""
    [1] pdf([S::Type{<:AbstractFloat}], d::AbstractDistribution, x)

Computes the probability density function (PDF) of the distribution `d` at the value `x`. 
Additionally accepts a type `S` to cast the result to a specific floating point type.

Note that the format of the value `x` depends on the distribution type.
"""
pdf(S::Type{<:AbstractFloat}, d::AbstractDistribution, args...) = S(pdf(d, args...))


"""
    [1] rand([S::Type{<:AbstractFloat}], d::AbstractDistribution, [N::Integer])

`Base.rand` Extension. Generates a random value from the distribution `d`. Additionally 
accepts a type `S` to cast the result to a specific floating point type, and an integer `N` 
to generate `N` random samples automatically.

Note that the out format of the random sample depends on the distribution type.
"""
Base.rand(S::Type{<:AbstractFloat}, d::AbstractDistribution) = S(rand(d))
Base.rand(d::AbstractDistribution, N::Integer) = [rand(d) for _ in 1:N]
Base.rand(S::Type{<:AbstractFloat}, d::AbstractDistribution, N::Integer) = S.(rand(d, N))


"""
    [1] rms([S::Type{<:AbstractFloat}], d::AbstractDistribution)

Computes the root mean square (RMS) of the distribution `d`. Additionally accepts a type
`S` to cast the result to a specific floating point type.
"""
rms(S::Type{<:AbstractFloat}, d::AbstractDistribution) = S(rms(d))


############################################################################################
#::. EXPORTS
############################################################################################
export cdf, mode, pdf, rms