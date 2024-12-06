############################################################################################
#::. STRUCTS
############################################################################################
"""
    PSDSpeedDistribution(T, m, [v0=575], [k=1.7f0])

Custom struct defining a Weibull distribution to describe the speed distribution of
particles released by photon-stimulated-desorption (PSD), based on the works of Gamborino
et al. (2018). The distribution is defined by the surface temperature `T`, the particle's
mass `m`, the offset speed `v0`, and the shape parameter `k`.

The distribution is zero for speeds lower than `v0`.

Defined following functions are defined for the distribution: `pdf`, `cdf`, `rand`,
`Statistics.mean`, `mode`.

# Arguments
- `T::Real`: Surface temperature (K)
- `m::Real`: Mass of the particle (kg)
- `v0::Real=575`: Offset speed (m s-1); default value based on Gamborino (2018)
- `k::Real=1.7f0`: Shape parameter (-); default value based on Gamborino (2018)

# References
- Gamborino, D., & Wurz, P. (2018). Velocity distribution function of Na released by photons
  from planetary surfaces. Planetary and Space Science, 159, 97-104.
  DOI: 10.1016/j.pss.2018.04.021
"""
struct PSDSpeedDistribution{S<:AbstractFloat} <: AbstractDistribution
    T::S; m::S; v0::S; k::S
end
function PSDSpeedDistribution(T::Real, m::Real, v0::Real, k::Real)
    return PSDSpeedDistribution(promote(T, m, v0, k)...)
end
function PSDSpeedDistribution(T::Integer, m::Integer, v0::Integer, k::Integer)
    return PSDSpeedDistribution(float(T), float(m), float(v0), float(k))
end
function PSDSpeedDistribution(T::Real, m::Real)
    return PSDSpeedDistribution(promote(T, m, 575, 1.7f0)...)
end


############################################################################################
#::. FUNCTIONS
############################################################################################
function cdf(d::PSDSpeedDistribution{S}, l::Real, u::Real) where {S<:AbstractFloat}
    if l <= d.v0; l = v0; end
    if u <= d.v0; return zero(S); end

    a = sqrt(d.m / (3 * BOLTZMANN_CONSTANT * d.T))
    g = gamma(1 + 1/d.k)

    return S(1 - exp(-((x-d.v0) * a * g)^d.k))
end


function Statistics.mean(d::PSDSpeedDistribution{S}) where {S<:AbstractFloat}
    return S(sqrt(d.m / (3 * BOLTZMANN_CONSTANT * d.T)) + d.v0)
end


function mode(d::PSDSpeedDistribution{S}) where {S<:AbstractFloat}
    if d.k <= 1; return d.v0; end

    a = sqrt(d.m / (3 * BOLTZMANN_CONSTANT * d.T))
    g = gamma(1 + 1/d.k)

    return S(d.v0 + ((d.k - 1)/d.k)^(1/d.k) / (a*g))
end


function pdf(d::PSDSpeedDistribution{S}, x::Real) where {S<:AbstractFloat}
    if x <= d.v0; return zero(S); end

    a = sqrt(d.m / (3 * BOLTZMANN_CONSTANT * d.T))
    g = gamma(1 + 1/d.k)

    return S(d.k * a * g * ((x-d.v0) * a * g)^(d.k-1) * exp(-((x-d.v0) * a * g)^d.k))
end


function Base.rand(d::PSDSpeedDistribution{S}) where {S<:AbstractFloat}
    a = sqrt(d.m / (3 * BOLTZMANN_CONSTANT * d.T))
    g = gamma(1 + 1/d.k)

    return S(d.v0 + (-log(1 - rand()))^(1/d.k) / (a*g))
end


# rms



############################################################################################
#::. EXPORTS
############################################################################################
export
    PSDSpeedDistribution
