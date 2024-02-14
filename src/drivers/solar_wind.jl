#::. STRUCTS
"""
    [1] SWPositionDistribution{S<:AbstractFloat}
    [2] SWPositionDistribution(r::Real)

Custom struct defining a surface distribution based on the influence of the solar wind.
Uses the radius `r` in [m] of the central planetary body as input. Applies only to the
Sun-facing side (i.e. `-pi/2 <= theta <= pi/2`).
"""
struct SWPositionDistribution{S<:AbstractFloat} <: AbstractSWPositionDistribution
    r::S
    # include angular offset of solar wind?
end
SWPositionDistribution(r::Integer) = SWPositionDistribution(Float64(r))
SWPositionDistribution(r::BigInt) = SWPositionDistribution(BigFloat(r))


#::. FUNCTIONS
function cdf(::SWPositionDistribution{S}, theta1::Real, theta2::Real, 
             phi1::Real, phi2::Real) where {S<:AbstractFloat}
    theta1 = theta1 < -pi/2 ? -pi/2 : theta1 # handle night side
    theta2 = theta2 > pi/2 ? pi/2 : theta2 # handle night side
    return S(1/4 * (sin(theta2) - sin(theta1)) * (sin(phi2) - sin(phi2)))
end
function cdf(d::SWPositionDistribution, theta::Tuple{Real, Real}, 
             phi::Tuple{Real, Real})
    return cdf(d, theta[1], theta[2], phi[1], phi[2])
end
function cdf(d::SWPositionDistribution, theta::Real, phi::Real)
    return cdf(d, -pi/2, theta, -pi/2, phi)
end
function cdf(d::SWPositionDistribution, 
             x1::Union{Tuple, AbstractVector, GlobalSphericalPosition},
             x2::Union{Tuple, AbstractVector, GlobalSphericalPosition})
    return cdf(d, _gettheta(x1), _gettheta(x2), _getphi(x1), _getphi(x2))
end
function cdf(d::SWPositionDistribution, 
             x::Union{Tuple, AbstractVector, GlobalSphericalPosition})
    return cdf(d, _gettheta(x), _getphi(x))
end
function cdf(d::SWPositionDistribution, x::AbstractPosition)
    return cdf(d, GlobalSphericalPosition(x))
end
cdf(S::Type{<:AbstractFloat}, args...) = S(cdf(args...))


# pdf

#::. extensions
function Base.rand(d::SWPositionDistribution) where {S<:AbstractFloat} 
    return (d.r, asin(2*rand(S)-1), asin(2*rand(S)-1))
end


#::. EXPORTS
export
    SWPositionDistribution

