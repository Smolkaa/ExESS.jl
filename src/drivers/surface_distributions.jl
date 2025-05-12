############################################################################################
#::. STRUCTS
############################################################################################
"""
    EqualSurfaceDistribution{S<:AbstractFloat}
    EqualSurfaceDistribution(r)

Custom struct defining an equal surface distribution based on the surface area of a sphere.
If sampled, returns a uniformly distributed random point on the surface of a sphere in 
spherical coordinates (radius, longitude, latitude).

# Arguments
- `r::Real`: Radius of the central planetary body (m)

# Usage
The distribution is to be used with the statistic functions: 
- `cdf(..., l, u)` to evaluate the cumulative distribution function between two points `l` 
  and `u`, where `l` and `u` are `Tuple`s or `Vector`s in spherical coordinates (radius, 
  longitude, latitude) or `GlobalSphericalPosition` or `GlobalCartesianPosition`.
- `pdf(..., x)` to evaluate the probability density function at a given point `x`, where `x` 
  is a `Tuple` or `Vector` in spherical coordinates (radius, longitude, latitude) or 
  `GlobalSphericalPosition` or `GlobalCartesianPosition`.
- `rand([S], ..., [N])` to draw a uniformly random point on a spherical surface in 
  spherical coordinates (radius, longitude, latitude). This function can be called with an 
  optional type `S` to cast the output to that type. If `N` is provided, it will return
  `N` samples in a `Vector` of `S`.
"""
struct EqualSurfaceDistribution{S<:AbstractFloat} <: AbstractDistribution; r::S; end
EqualSurfaceDistribution(r::Integer) = EqualSurfaceDistribution(Float64(r))
EqualSurfaceDistribution(r::BigInt) = EqualSurfaceDistribution(BigFloat(r))


"""
    SolarSurfaceDistribution{S<:AbstractFloat}
    SolarSurfaceDistribution(r)

Custom struct defining a surface distribution based on the influence of the solar wind.
If sampled, returns a solar-wind distributed random point on the Sun-facing surface 
(i.e. subsolar longitude in `-pi/2 <= lon <= pi/2`) of a sphere in  spherical coordinates 
(radius, longitude, latitude).

# Arguments
- `r::Real`: Radius of the central planetary body (m)

# Usage
The distribution is to be used with the statistic functions:
- `cdf(..., l, u)` to evaluate the cumulative distribution function between two points `l` 
  and `u`, where `l` and `u` are `Tuple`s or `Vector`s in spherical coordinates (radius, 
  longitude, latitude) or `GlobalSphericalPosition` or `GlobalCartesianPosition`.
- `pdf(..., x)` to evaluate the probability density function at a given point `x`, where `x`
  is a `Tuple` or `Vector` in spherical coordinates (radius, longitude, latitude) or 
  `GlobalSphericalPosition` or `GlobalCartesianPosition`.
- `rand([S], ..., [N])` to draw a uniformly random point on a solar-wind distributed
  spherical surface in spherical coordinates (radius, longitude, latitude). This function 
  can be called with an optional type `S` to cast the output to that type. If `N` is provided,
  it will return `N` samples in a `Vector` of `S`.
"""
struct SolarSurfaceDistribution{S<:AbstractFloat} <: AbstractDistribution; r::S; end
SolarSurfaceDistribution(r::Integer) = SolarSurfaceDistribution(Float64(r))
SolarSurfaceDistribution(r::BigInt) = SolarSurfaceDistribution(BigFloat(r))



############################################################################################
#::. FUNCTIONS

# internal union for simplified type handling
_TVGSP = Union{Tuple{Real, Real, Real}, AbstractVector{<:Real}, GlobalSphericalPosition}
_SDs = Union{EqualSurfaceDistribution, SolarSurfaceDistribution}
############################################################################################
function cdf(::EqualSurfaceDistribution{S}, l::Tuple{Real, Real, Real},
             u::Tuple{Real, Real, Real}) where {S<:AbstractFloat}
    lon1, lon2, lat1, lat2 = l[2], u[2], l[3], u[3]
    return S(1/4/pi * (lon2 - lon1) * (sin(lat2) - sin(lat1)))
end
cdf(d::EqualSurfaceDistribution, l::_TVGSP, u::_TVGSP) = cdf(d, Tuple(l), Tuple(u))
cdf(d::EqualSurfaceDistribution, u::_TVGSP) = cdf(d, (u[1], -pi, -pi/2), u)

function cdf(::SolarSurfaceDistribution{S}, l::Tuple{Real, Real, Real},
             u::Tuple{Real, Real, Real}) where {S<:AbstractFloat}

    lon1 = l[2] < -pi/2 ? -pi/2 : l[2] # handle night side
    lon2 = u[2] >  pi/2 ?  pi/2 : u[2] # handle night side
    lat1, lat2 = l[3], u[3] # extract latitudes
    return S(1/4 * (sin(lon2) - sin(lon1)) * (sin(lat2) - sin(lat1)))
end
cdf(d::SolarSurfaceDistribution, l::_TVGSP, u::_TVGSP) = cdf(d, Tuple(l), Tuple(u))
cdf(d::SolarSurfaceDistribution, u::_TVGSP) = cdf(d, (u[1], -pi/2, -pi/2), u)



function pdf(::EqualSurfaceDistribution{S},
             x::Tuple{Real, Real, Real}) where {S<:AbstractFloat}
    return S(1/4/pi * cos(x[3]))
end
pdf(d::EqualSurfaceDistribution, x::_TVGSP) = pdf(d, Tuple(x))

function pdf(::SolarSurfaceDistribution{S},
             x::Tuple{Real, Real, Real}) where {S<:AbstractFloat}
    lon, lat = x[2], x[3]
    if lon < -pi/2 || lon > pi/2; return zero(S); end # handle night side
    return S(1/4 * cos(lon) * cos(lat))
end
pdf(d::SolarSurfaceDistribution, x::_TVGSP) = pdf(d, Tuple(x))


function Base.rand(d::EqualSurfaceDistribution{S}) where {S<:AbstractFloat}
    return (d.r, S(rand()*2pi - pi), asin(2*rand(S)-1))
end
function Base.rand(d::SolarSurfaceDistribution{S}) where {S<:AbstractFloat}
    return (d.r, asin(2*rand(S)-1), asin(2*rand(S)-1))
end

# overwriting additional `Base.rand` methods for multivariate surface distributions `_SDs`
Base.rand(S::Type{<:AbstractFloat}, d::_SDs) = S.(rand(d))
Base.rand(S::Type{<:GlobalSphericalPosition}, d::_SDs) = S(rand(d))
Base.rand(S::Type{<:GlobalCartesianPosition}, d::_SDs) = S(rand(GlobalSphericalPosition, d))
function Base.rand(S::Type{<:AbstractFloat}, d::_SDs, N::Integer) # performance improvement
    SAMPLES = Vector{NTuple{3, S}}(undef, N)
    for i in 1:N; SAMPLES[i] = rand(S, d); end
    return SAMPLES
end
Base.rand(S::Type{<:GlobalSphericalPosition}, d::_SDs, N::Integer) = S.(rand(d, N))
Base.rand(S::Type{<:GlobalCartesianPosition}, d::_SDs, N::Integer) = S.(rand(GlobalSphericalPosition, d, N))


############################################################################################
#::. EXPORTS
############################################################################################
export EqualSurfaceDistribution, SolarSurfaceDistribution
