############################################################################################
#::. STRUCTS
############################################################################################
"""
    [1] SWSurfaceDistribution{S<:AbstractFloat}
    [2] SWSurfaceDistribution(r::Real)

Custom struct defining a surface distribution based on the influence of the solar wind.
Uses the radius `r` in [m] of the central planetary body as input. Applies only to the
Sun-facing side (i.e. `-pi/2 <= theta <= pi/2`).
"""
struct SWSurfaceDistribution{S<:AbstractFloat} <: AbstractDistribution; r::S; end
SWSurfaceDistribution(r::Integer) = SWSurfaceDistribution(Float64(r))
SWSurfaceDistribution(r::BigInt) = SWSurfaceDistribution(BigFloat(r))



############################################################################################
#::. FUNCTIONS

# internal union for simplified type handling
_TVGSP = Union{Tuple{Real, Real, Real}, AbstractVector{<:Real}, GlobalSphericalPosition}
############################################################################################
function cdf(::SWSurfaceDistribution{S}, l::Tuple{Real, Real, Real}, 
             u::Tuple{Real, Real, Real}) where {S<:AbstractFloat}

    lon1 = l[2] < -pi/2 ? -pi/2 : l[2] # handle night side
    lon2 = u[2] >  pi/2 ?  pi/2 : u[2] # handle night side
    lat1, lat2 = l[3], u[3] # extract latitudes
    return S(1/4 * (sin(lon2) - sin(lon1)) * (sin(lat2) - sin(lat1)))
end
cdf(d::SWSurfaceDistribution, l::_TVGSP, u::_TVGSP) = cdf(d, _get(l), _get(u))
cdf(d::SWSurfaceDistribution, u::_TVGSP) = cdf(d, (u[1], -pi/2, -pi/2), u)


function pdf(::SWSurfaceDistribution{S}, x::Tuple{Real, Real, Real}) where {S<:AbstractFloat}
    lon, lat = x[2], x[3]
    if lon < -pi/2 || lon > pi/2; return zero(S); end # handle night side
    return S(1/4 * cos(lon) * cos(lat))
end
pdf(d::SWSurfaceDistribution, x::_TVGSP) = pdf(d, _get(x))


function Base.rand(d::SWSurfaceDistribution{S}) where {S<:AbstractFloat} 
    return (d.r, asin(2*rand(S)-1), asin(2*rand(S)-1))
end


############################################################################################
#::. EXPORTS
############################################################################################
export
    SWSurfaceDistribution

