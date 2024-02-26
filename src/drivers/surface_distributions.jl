############################################################################################
#::. STRUCTS
############################################################################################
"""
    [1] EqualSurfaceDistribution{S<:AbstractFloat}
    [2] EqualSurfaceDistribution(r::Real)

Custom struct defining an equal surface distribution based on the surface area of a sphere. 
Uses the radius `r` in [m] of the central planetary body as input.
"""
struct EqualSurfaceDistribution{S<:AbstractFloat} <: AbstractDistribution; r::S; end
EqualSurfaceDistribution(r::Integer) = EqualSurfaceDistribution(Float64(r))
EqualSurfaceDistribution(r::BigInt) = EqualSurfaceDistribution(BigFloat(r))


"""
    [1] SolarSurfaceDistribution{S<:AbstractFloat}
    [2] SolarSurfaceDistribution(r::Real)

Custom struct defining a surface distribution based on the influence of the solar wind.
Uses the radius `r` in [m] of the central planetary body as input. Applies only to the
Sun-facing side (i.e. subsolar longitude in `-pi/2 <= lon <= pi/2`).
"""
struct SolarSurfaceDistribution{S<:AbstractFloat} <: AbstractDistribution; r::S; end
SolarSurfaceDistribution(r::Integer) = SolarSurfaceDistribution(Float64(r))
SolarSurfaceDistribution(r::BigInt) = SolarSurfaceDistribution(BigFloat(r))



############################################################################################
#::. FUNCTIONS

# internal union for simplified type handling
_TVGSP = Union{Tuple{Real, Real, Real}, AbstractVector{<:Real}, GlobalSphericalPosition}
############################################################################################
function cdf(::EqualSurfaceDistribution{S}, l::Tuple{Real, Real, Real}, 
             u::Tuple{Real, Real, Real}) where {S<:AbstractFloat}
    lon1, lon2, lat1, lat2 = l[2], u[2], l[3], u[3]
    return S(1/4/pi * (lon2 - lon1) * (sin(lat2) - sin(lat1)))
end
cdf(d::EqualSurfaceDistribution, l::_TVGSP, u::_TVGSP) = cdf(d, _get(l), _get(u))
cdf(d::EqualSurfaceDistribution, u::_TVGSP) = cdf(d, (u[1], -pi, -pi/2), u)

function cdf(::SolarSurfaceDistribution{S}, l::Tuple{Real, Real, Real}, 
             u::Tuple{Real, Real, Real}) where {S<:AbstractFloat}

    lon1 = l[2] < -pi/2 ? -pi/2 : l[2] # handle night side
    lon2 = u[2] >  pi/2 ?  pi/2 : u[2] # handle night side
    lat1, lat2 = l[3], u[3] # extract latitudes
    return S(1/4 * (sin(lon2) - sin(lon1)) * (sin(lat2) - sin(lat1)))
end
cdf(d::SolarSurfaceDistribution, l::_TVGSP, u::_TVGSP) = cdf(d, _get(l), _get(u))
cdf(d::SolarSurfaceDistribution, u::_TVGSP) = cdf(d, (u[1], -pi/2, -pi/2), u)



function pdf(::EqualSurfaceDistribution{S}, 
             x::Tuple{Real, Real, Real}) where {S<:AbstractFloat}
    return S(1/4/pi * cos(x[3]))
end
pdf(d::EqualSurfaceDistribution, x::_TVGSP) = pdf(d, _get(x))

function pdf(::SolarSurfaceDistribution{S}, 
             x::Tuple{Real, Real, Real}) where {S<:AbstractFloat}
    lon, lat = x[2], x[3]
    if lon < -pi/2 || lon > pi/2; return zero(S); end # handle night side
    return S(1/4 * cos(lon) * cos(lat))
end
pdf(d::SolarSurfaceDistribution, x::_TVGSP) = pdf(d, _get(x))


function Base.rand(d::EqualSurfaceDistribution{S}) where {S<:AbstractFloat} 
    return (d.r, S(rand()*2pi - pi), asin(2*rand(S)-1))
end
function Base.rand(d::SolarSurfaceDistribution{S}) where {S<:AbstractFloat} 
    return (d.r, asin(2*rand(S)-1), asin(2*rand(S)-1))
end


############################################################################################
#::. EXPORTS
############################################################################################
export EqualSurfaceDistribution, SolarSurfaceDistribution

