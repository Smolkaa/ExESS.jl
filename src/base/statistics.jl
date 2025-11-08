############################################################################################
#::. STRUCTS
############################################################################################
abstract type AbstractDistribution end


############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    cdf([S], d, [l], u) -> Real

Computes the cumulative distribution function (CDF) of the distribution `d` between the
lower limit `l` and upper limit `u`. If only the upper limit `u` is provided, the lower
limit is assumed to be the minimum value of the distribution. Additionally accepts a type
`S` to cast the result to a specific floating point type. The format of the limits `l` and
`u` depend on the distribution type.

# Arguments
- `d::AbstractDistribution`: distribution to compute the CDF for
- `l`: lower limit of the CDF computation (optional)
- `u`: upper limit of the CDF computation
- `S::Type{<:AbstractFloat}`: floating point type to cast the result to (optional)
"""
cdf(S::Type{<:AbstractFloat}, d::AbstractDistribution, args...) = S(cdf(d, args...))


"""
    mean([S], d) -> Real

`Statistics.mean` extension. Computes the mean of the distribution `d`. Additionally accepts
a type `S` to cast the result to a specific floating point type.

# Arguments
- `d::AbstractDistribution`: distribution to compute the mean for
- `S::Type{<:AbstractFloat}`: floating point type to cast the result to (optional)
"""
Statistics.mean(S::Type{<:AbstractFloat}, d::AbstractDistribution) = S(mean(d))


"""
    mode([S], d) -> Real

Computes the mode of the distribution `d`. Additionally accepts a type `S` to cast the
result to a specific floating point type.

# Arguments
- `d::AbstractDistribution`: distribution to compute the mode for
- `S::Type{<:AbstractFloat}`: floating point type to cast the result to (optional)
"""
mode(S::Type{<:AbstractFloat}, d::AbstractDistribution) = S(mode(d))


"""
    pdf([S], d, x) -> Real

Computes the probability density function (PDF) of the distribution `d` at the value `x`. 
Additionally accepts a type `S` to cast the result to a specific floating point type. The 
format of the value `x` depends on the distribution type.

# Arguments
- `d::AbstractDistribution`: distribution to compute the PDF for
- `x`: value at which to compute the PDF
- `S::Type{<:AbstractFloat}`: floating point type to cast the result to (optional)
"""
pdf(S::Type{<:AbstractFloat}, d::AbstractDistribution, args...) = S(pdf(d, args...))


"""
    rand([S], d, [N])

`Base.rand` extension. Generates a random value from the distribution `d`. Additionally 
accepts a type `S` to cast the result to a specific floating point type, and an integer `N` 
to generate `N` random samples automatically.

# Arguments
- `d::AbstractDistribution`: distribution to generate the random sample from
- `N::Integer`: number of random samples to generate (optional)
- `S::Type{<:AbstractFloat}`: floating point type to cast the result to (optional)
"""
Base.rand(S::Type{<:AbstractFloat}, d::AbstractDistribution) = S(rand(d))
Base.rand(d::AbstractDistribution, N::Integer) = [rand(d) for _ in 1:N]
Base.rand(S::Type{<:AbstractFloat}, d::AbstractDistribution, N::Integer) = S.(rand(d, N))


"""
    rms([S], d) -> Real

Computes the root mean square (RMS) of the distribution `d`. Additionally accepts a type
`S` to cast the result to a specific floating point type.

# Arguments
- `d::AbstractDistribution`: distribution to compute the RMS for
- `S::Type{<:AbstractFloat}`: floating point type to cast the result to (optional)
"""
rms(S::Type{<:AbstractFloat}, d::AbstractDistribution) = S(rms(d))


############################################################################################
#::. EXPORTS
############################################################################################
export cdf, mode, pdf, rms