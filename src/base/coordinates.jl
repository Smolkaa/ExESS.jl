#::. STRUCTS
abstract type AbstractCoordinate end
abstract type AbstractPosition{T<:AbstractFloat} <: AbstractCoordinate end
abstract type AbstractVelocity{T<:AbstractFloat} <: AbstractCoordinate end

struct GlobalCartesianPosition{T<:AbstractFloat} <: AbstractPosition{T}
    x::T
    y::T
    z::T
end
function GlobalCartesianPosition(x::Real, y::Real, z::Real)
    GlobalCartesianPosition(promote(x, y, z)...)
end
function GlobalCartesianPosition(x::Integer, y::Integer, z::Integer)
    GlobalCartesianPosition(promote(x, y, z, 1.0)[1:3]...)
end

struct LocalCartesianPosition{T<:AbstractFloat} <: AbstractPosition{T}
    x::T
    y::T
    z::T
end
function LocalCartesianPosition(x::Real, y::Real, z::Real)
    LocalCartesianPosition(promote(x, y, z)...)
end
function LocalCartesianPosition(x::Integer, y::Integer, z::Integer)
    LocalCartesianPosition(promote(x, y, z, 1.0)[1:3]...)
end

struct GlobalSphericalPosition{T<:AbstractFloat} <: AbstractPosition{T}
    r::T
    theta::T
    phi::T
end
function GlobalSphericalPosition(r::Real, theta::Real, phi::Real)
    GlobalSphericalPosition(promote(r, theta, phi)...)
end
function GlobalSphericalPosition(r::Integer, theta::Integer, phi::Integer)
    GlobalSphericalPosition(promote(r, theta, phi, 1.0)[1:3]...)
end

struct GlobalCartesianVelocity{T<:AbstractFloat} <: AbstractVelocity{T}
    x::T
    y::T
    z::T
end
function GlobalCartesianVelocity(x::Real, y::Real, z::Real)
    GlobalCartesianVelocity(promote(x, y, z)...)
end
function GlobalCartesianVelocity(x::Integer, y::Integer, z::Integer)
    GlobalCartesianVelocity(promote(x, y, z, 1.0)[1:3]...)
end

struct LocalCartesianVelocity{T<:AbstractFloat} <: AbstractVelocity{T}
    x::T # (E)
    y::T # (N)
    z::T
end
function LocalCartesianVelocity(x::Real, y::Real, z::Real)
    LocalCartesianVelocity(promote(x, y, z)...)
end
function LocalCartesianVelocity(x::Integer, y::Integer, z::Integer)
    LocalCartesianVelocity(promote(x, y, z, 1.0)[1:3]...)
end

struct GlobalSphericalVelocity{T<:AbstractFloat} <: AbstractVelocity{T}
    r::T
    theta::T
    phi::T
end
function GlobalSphericalVelocity(r::Real, theta::Real, phi::Real)
    GlobalSphericalVelocity(promote(r, theta, phi)...)
end
function GlobalSphericalVelocity(r::Integer, theta::Integer, phi::Integer)
    GlobalSphericalVelocity(promote(r, theta, phi, 1.0)[1:3]...)
end


#::. (useful) unions
const AbstractGlobalCoordinate = Union{
    GlobalCartesianPosition, 
    GlobalSphericalPosition, 
    GlobalCartesianVelocity, 
    GlobalSphericalVelocity,}
const AbstractLocalCoordinate = Union{LocalCartesianPosition, LocalCartesianVelocity}

const AbstractGlobalPosition = Union{GlobalCartesianPosition, GlobalSphericalPosition}
const AbstractGlobalVelocity = Union{GlobalCartesianVelocity, GlobalSphericalVelocity}

const AbstractLocalPosition = Union{LocalCartesianPosition}
const AbstractLocalVelocity = Union{LocalCartesianVelocity}

const AbstractCartesianCoordinates = Union{
    GlobalCartesianPosition, 
    LocalCartesianPosition, 
    GlobalCartesianVelocity, 
    LocalCartesianVelocity,}
const AbstractSphericalCoordinates = Union{GlobalSphericalPosition, GlobalSphericalVelocity}


#::. FUNCTIONS
"""
    [1] GlobalCartesianPosition(x::Real, y::Real, z::Real)
    [2] GlobalCartesianPosition(X::Tuple)
    [3] GlobalCartesianPosition(X::AbstractVector)
    [4] GlobalCartesianPosition(x::AbstractGlobalPosition)

Three dimensional, global, cartesian position vector [1].
Converts tuples, vector or other position type into `GlobalCartesianPosition` [2,3,4].

In a planetary context, global cartesian coordinates are defined as follows:
- the x-axis points towards the Sun (from the central body's center to the subsolar point)
- the y-axis points towards the positive longitude (east) 
- the z-axis points towards the positive latitude (north), completing the coordinate system
"""
GlobalCartesianPosition(X::Tuple) = GlobalCartesianPosition(X...)
GlobalCartesianPosition(X::AbstractVector) = GlobalCartesianPosition(X...)
GlobalCartesianPosition(x::GlobalCartesianPosition) = x
function GlobalCartesianPosition(x::GlobalSphericalPosition{T}) where {T<:AbstractFloat}
    return GlobalCartesianPosition{T}(
        x.r * cos(x.theta) * cos(x.phi),
        x.r * sin(x.theta) * cos(x.phi),
        x.r * sin(x.phi),
    )
end


"""
    [1] LocalCartesianPosition(x::Real, y::Real, z::Real)
    [2] LocalCartesianPosition(X::Tuple)
    [3] LocalCartesianPosition(X::AbstractVector)
    [4] LocalCartesianPosition(x::AbstractLocalPosition)

Three dimensional, local, cartesian position vector [1].
Converts tuples, vector or other position type into `LocalCartesianPosition` [2,3,4].

Local cartesian coordinates are defined as follows:
- the x- and y-axis are spanning the tangent plane at the given position
- the z-axis is pointing upwards
"""
LocalCartesianPosition(X::Tuple) = LocalCartesianPosition(X...)
LocalCartesianPosition(X::AbstractVector) = LocalCartesianPosition(X...)
LocalCartesianPosition(x::LocalCartesianPosition) = x


"""
    [1] GlobalSphericalPosition(r::Real, theta::Real, phi::Real)
    [2] GlobalSphericalPosition(X::Tuple)
    [3] GlobalSphericalPosition(X::AbstractVector)
    [4] GlobalSphericalPosition(x::AbstractGlobalPosition)

Three dimensional, global, spherical position vector [1].
Converts tuples, vector or other position type into `GlobalSphericalPosition` [2,3,4].

In a planetary context, global spherical coordinates are defined as follows:
- the r-axis points towards the Sun (from the central body's center to the subsolar point)
- the theta-angle describes positive longitudes (east)
- the phi-angle describes positive latitudes (north)

The subsolar point is defined as `R, 0, 0` where `R` is the radius of the central body.
"""
GlobalSphericalPosition(X::Tuple) = GlobalSphericalPosition(X...)
GlobalSphericalPosition(X::AbstractVector) = GlobalSphericalPosition(X...)
GlobalSphericalPosition(x::GlobalSphericalPosition) = x
function GlobalSphericalPosition(x::GlobalCartesianPosition{T}) where {T<:AbstractFloat}
    r     = norm(x)
    theta = atan(x.y / x.x)
    phi   = asin(x.z / r)
    if isnan(theta); theta = 0;    end
    if x.x < 0;      theta += pi;  end
    if theta > pi;   theta -= 2pi; end
    return GlobalSphericalPosition{T}(r, theta, phi)
end


"""
    [1] GlobalCartesianVelocity(x::Real, y::Real, z::Real)
    [3] GlobalCartesianVelocity(V::Tuple)
    [4] GlobalCartesianVelocity(V::AbstractVector)
    [5] GlobalCartesianVelocity(v::GlobalCartesianVelocity)
    [6] GlobalCartesianVelocity(x::AbstractGlobalPosition, v::AbstractVelocity)

Three dimensional velocity vector in a global cartesian coordinate system [1].
Converts tuples, vector or other velocity type into `GlobalCartesianVelocity` [2,3,4,5].

In a planetary context, global cartesian coordinates are defined as follows:
- the x-axis points towards the Sun (from the central body's center to the subsolar point)
- the y-axis points towards the positive longitude (east)
- the z-axis points towards the positive latitude (north), completing the coordinate system

**Notes**
- type-stability cannot be guaranteed for integer subtypes of `AbstractGlobalPosition` in [5]
"""
GlobalCartesianVelocity(V::Tuple) = GlobalCartesianVelocity(V...)
GlobalCartesianVelocity(V::AbstractVector) = GlobalCartesianVelocity(V...)
GlobalCartesianVelocity(v::GlobalCartesianVelocity) = v
GlobalCartesianVelocity(::AbstractPosition, v::GlobalCartesianVelocity) = v
function GlobalCartesianVelocity(x::AbstractPosition, v::AbstractVelocity)
    return GlobalCartesianVelocity(GlobalSphericalPosition(x), v)
end
function GlobalCartesianVelocity(x::GlobalSphericalPosition, v::LocalCartesianVelocity)
    return GlobalCartesianVelocity(
      - v.x*sin(x.theta) - v.y*cos(x.theta)*sin(x.phi) + v.z*cos(x.theta)*cos(x.phi),
        v.x*cos(x.theta) - v.y*sin(x.theta)*sin(x.phi) + v.z*sin(x.theta)*cos(x.phi),
                           v.y*             cos(x.phi) + v.z*             sin(x.phi),
    )
end
function GlobalCartesianVelocity(x::GlobalSphericalPosition, v::GlobalSphericalVelocity)
    return GlobalCartesianVelocity(
        v.r*cos(x.theta)*cos(x.phi) - v.theta*x.r*sin(x.theta)*cos(x.phi) - v.phi*x.r*cos(x.theta)*sin(x.phi),
        v.r*sin(x.theta)*cos(x.phi) + v.theta*x.r*cos(x.theta)*cos(x.phi) - v.phi*x.r*sin(x.theta)*sin(x.phi),
        v.r             *sin(x.phi)                                       + v.phi*x.r             *cos(x.phi),
    )
end


"""
    [1] LocalCartesianVelocity(x::Real, y::Real, z::Real)
    [2] LocalCartesianVelocity(V::Tuple)
    [3] LocalCartesianVelocity(V::AbstractVector)
    [4] LocalCartesianVelocity(v::LocalCartesianVelocity)
    [5] LocalCartesianVelocity(x::AbstractGlobalPosition, v::AbstractVelocity)

Three dimensional velocity vector in a local cartesian coordinate system [1].
Converts tuples, vector or other velocity type into `LocalCartesianVelocity` [2,3,4,5].

In planetary contexts, local cartesian coordinates are defined as follows:
- the x-axis is tangent to the surface, pointing towards positive longitude (east)
- the y-axis is tangent to the surface, pointing towards positive latitude (north)
- the z-axis is pointing radially outwards

**Notes**
- type-stability cannot be guaranteed for integer subtypes of `AbstractGlobalPosition` in [5]
"""
LocalCartesianVelocity(V::Tuple) = LocalCartesianVelocity(V...)
LocalCartesianVelocity(V::AbstractVector) = LocalCartesianVelocity(V...)
LocalCartesianVelocity(v::LocalCartesianVelocity) = v
LocalCartesianVelocity(::AbstractPosition, v::LocalCartesianVelocity) = v
LocalCartesianVelocity(x::AbstractPosition, v::AbstractVelocity) = LocalCartesianVelocity(GlobalSphericalPosition(x), v)
function LocalCartesianVelocity(x::GlobalCartesianPosition, v::GlobalCartesianVelocity)
    xy, xyz = sqrt(x.x^2 + x.y^2), norm(x)
    sintheta, costheta = x.y / xy, x.x / xy
    sinphi, cosphi = x.z / xyz, xy / xyz
    return LocalCartesianVelocity(
       -v.x*sintheta        + v.y*costheta,
       -v.x*costheta*sinphi - v.y*sintheta*sinphi + v.z*cosphi,
        v.x*costheta*cosphi + v.y*sintheta*cosphi + v.z*sinphi
    )
end
function LocalCartesianVelocity(x::GlobalSphericalPosition, v::GlobalCartesianVelocity)
    return LocalCartesianVelocity(
       -v.x*sin(x.theta)            + v.y*cos(x.theta),
       -v.x*cos(x.theta)*sin(x.phi) - v.y*sin(x.theta)*sin(x.phi) + v.z*cos(x.phi),
        v.x*cos(x.theta)*cos(x.phi) + v.y*sin(x.theta)*cos(x.phi) + v.z*sin(x.phi)
    )
end
function LocalCartesianVelocity(x::GlobalSphericalPosition, v::GlobalSphericalVelocity)
    return LocalCartesianVelocity(x.r*cos(x.phi)*v.theta, x.r*v.phi, v.r)
end


"""
    [1] GlobalSphericalVelocity(r::Real, theta::Real, phi::Real)
    [2] GlobalSphericalVelocity(V::Tuple)
    [3] GlobalSphericalVelocity(V::AbstractVector)
    [4] GlobalSphericalVelocity(v::GlobalSphericalVelocity)
    [5] GlobalSphericalVelocity(x::AbstractGlobalPosition, v::AbstractVelocity)

Three dimensional velocity vector in a global spherical coordinate system [1].
Converts tuples, vector or other velocity type into `GlobalSphericalVelocity` [2,3,4,5].

In a planetary context, global spherical coordinates are defined as follows:
- the r-axis points towards the Sun (from the central body's center to the subsolar point)
- the theta-angle describes positive longitudes (east)
- the phi-angle describes positive latitudes (north)

**Notes**
- type-stability cannot be guaranteed for integer subtypes of `AbstractGlobalPosition` in [5]
"""
GlobalSphericalVelocity(V::Tuple)= GlobalSphericalVelocity(V...)
GlobalSphericalVelocity(V::AbstractVector)= GlobalSphericalVelocity(V...)
GlobalSphericalVelocity(v::GlobalSphericalVelocity) = v
GlobalSphericalVelocity(::AbstractPosition, v::GlobalSphericalVelocity) = v
GlobalSphericalVelocity(x::AbstractPosition, v::AbstractVelocity) = GlobalSphericalVelocity(GlobalSphericalPosition(x), v)
function GlobalSphericalVelocity(x::GlobalSphericalPosition, v::GlobalCartesianVelocity)
    return GlobalSphericalVelocity(
        v.x*cos(x.theta)*cos(x.phi)     + v.y*sin(x.theta)*cos(x.phi)      + v.z*sin(x.phi),
       -v.x*sin(x.theta)/x.r/cos(x.phi) + v.y*cos(x.theta)/x.r/cos(x.phi),
       -v.x*cos(x.theta)*sin(x.phi)/x.r - v.y*sin(x.theta)*sin(x.phi)/x.r  + v.z*cos(x.phi)/x.r,
    )
end
function GlobalSphericalVelocity(x::GlobalSphericalPosition, v::LocalCartesianVelocity)
    return GlobalSphericalVelocity(v.z, v.x / x.r / cos(x.phi), v.y / x.r)
end



#::. utility
_get(a::T) where {T<:AbstractCartesianCoordinates} = (a.x, a.y, a.z)
_get(a::T) where {T<:AbstractSphericalCoordinates} = (a.r, a.theta, a.phi)
_get(a::Tuple) = a
_get(a::AbstractVector) = Tuple(a)

_getx(a::T) where {T<:AbstractCartesianCoordinates} = a.x
_getx(a::Union{Tuple, AbstractVector}) = a[1]

_gety(a::T) where {T<:AbstractCartesianCoordinates} = a.y
_gety(a::Union{Tuple, AbstractVector}) = a[2]

_getz(a::T) where {T<:AbstractCartesianCoordinates} = a.z
_getz(a::Union{Tuple, AbstractVector}) = a[3]

_getr(a::T) where {T<:AbstractSphericalCoordinates} = a.r
_getr(a::Union{Tuple, AbstractVector}) = a[1]

_gettheta(a::T) where {T<:AbstractSphericalCoordinates} = a.theta
_gettheta(a::Union{Tuple, AbstractVector}) = a[2]

_getphi(a::T) where {T<:AbstractSphericalCoordinates} = a.phi
_getphi(a::Union{Tuple, AbstractVector}) = a[3]


"""
    [1] azimuth([S::Type{<:AbstractFloat}], v::Tuple)
    [2] azimuth([S::Type{<:AbstractFloat}], v::AbstractVector)
    [3] azimuth([S::Type{<:AbstractFloat}], v::AbstractCartesianCoordinates)

Calculate the azimuth angle (in [rad]) of the cartesian vector `c`.
"""
azimuth(c::Union{Tuple, AbstractVector}) = sgn(c[2]) * acos( c[1] / sqrt(c[1]^2 + c[2]^2) )
azimuth(c::AbstractCartesianCoordinates) = sgn(c.y) * acos( c.x / sqrt(c.x^2 + c.y^2) )
azimuth(S::Type{<:AbstractFloat}, args...) = S(azimuth(args...))


"""
    [1] elevation([S::Type{<:AbstractFloat}], v::Tuple)
    [2] elevation([S::Type{<:AbstractFloat}], v::AbstractVector)
    [3] elevation([S::Type{<:AbstractFloat}], v::AbstractCartesianCoordinates)

Calculate the elevation angle (in [rad]) of the cartesian vector `c`.
"""
elevation(c::Union{Tuple, AbstractVector}) = atan(c[3] / sqrt(c[1]^2 + c[2]^2))
elevation(c::AbstractCartesianCoordinates) = atan(c.z / sqrt(c.x^2 + c.y^2))
elevation(S::Type{<:AbstractFloat}, args...) = S(elevation(args...))


"""
    [1] speed([S::Type], v::Tuple)
    [2] speed([S::Type], v::AbstractVector)
    [3] speed([S::Type], v::LocalCartesianVelocity)
    [4] speed([S::Type], v::GobalCartesianVelocity)

Calculate the speed (in [m s-1]) of the cartesian velocity vector `v` (in [m s-1]).
"""
speed(v::Union{Tuple, AbstractVector, LocalCartesianVelocity, GlobalCartesianVelocity}) = norm(v)
speed(S::Type{<:AbstractFloat}, args...) = S(speed(args...))


"""
    [1] zenith([S::Type{<:Real}], v::Tuple)
    [2] zenith([S::Type{<:Real}], v::AbstractVector)
    [3] zenith([S::Type{<:Real}], v::AbstractCartesianCoordinates)

Calculate the zenith angle (in [rad]) of the the cartesian vector `c`.

**Notes**
- the zenith angle is the same-sign `pi/2`-inversion of the elevation angle
"""
zenith(c::Union{Tuple, AbstractVector}) = atan(sqrt(c[1]^2 + c[2]^2) / c[3])
zenith(c::AbstractCartesianCoordinates) = atan(sqrt(c.x^2 + c.y^2) / c.z)
zenith(S::Type{<:AbstractFloat}, args...) = S(zenith(args...))



#::. EXTENSIONS
Base.:+(a::AbstractCoordinate) = a
Base.:+(a::T, b::Tuple) where {T<:AbstractCoordinate} = T(_get(a) .+ b ...)
Base.:+(a::T, b::AbstractVector) where {T<:AbstractCoordinate} = T(vec(a) .+ b ...)
Base.:+(a::Tuple, b::AbstractCoordinate) = b + a
Base.:+(a::AbstractVector, b::AbstractCoordinate) = b + a
Base.:+(a::GlobalCartesianPosition, b::GlobalCartesianPosition) = a + _get(b)
Base.:+(a::GlobalCartesianPosition, b::GlobalSphericalPosition) = a + GlobalCartesianPosition(b)
Base.:+(a::LocalCartesianPosition, b::LocalCartesianPosition) = a + _get(b)
Base.:+(a::GlobalSphericalPosition, b::GlobalSphericalPosition) = GlobalSphericalPosition(GlobalCartesianPosition(a) + GlobalCartesianPosition(b))
Base.:+(a::GlobalSphericalPosition, b::GlobalCartesianPosition) = GlobalSphericalPosition(GlobalCartesianPosition(a) + b)
Base.:+(a::GlobalCartesianVelocity, b::GlobalCartesianVelocity) = a + _get(b)
Base.:+(a::LocalCartesianVelocity, b::LocalCartesianVelocity) = a + _get(b)


Base.:-(a::T) where {T<:AbstractCoordinate} = T( .- _get(a) ...)
Base.:-(a::AbstractCoordinate, b::Any) = a + (-b)
Base.:-(a::AbstractCoordinate, b::Tuple) = a + (.-b)
Base.:-(a::AbstractVector, b::AbstractCoordinate) = b + (-a)
Base.:-(a::Tuple, b::AbstractCoordinate) = b + (.-a)


Base.:*(a::T, b::Real) where {T<:AbstractCoordinate} = T(_get(a) .* b ... )
Base.:*(a::T, b::Tuple) where {T<:AbstractCoordinate} = T(_get(a) .* b ... )
Base.:*(a::T, b::AbstractVector) where {T<:AbstractCoordinate} = T(vec(a) .* b ... )
Base.:*(a::Real, b::AbstractCoordinate) = b * a
Base.:*(a::Tuple, b::AbstractCoordinate) = b * a
Base.:*(a::AbstractVector, b::AbstractCoordinate) = b * a

Base.:*(a::GlobalCartesianPosition, b::GlobalCartesianPosition) = norm(_get(a) .* _get(b))
Base.:*(a::GlobalCartesianPosition, b::GlobalSphericalPosition) = a * GlobalCartesianPosition(b)
Base.:*(a::LocalCartesianPosition, b::LocalCartesianPosition) = norm(_get(a) .* _get(b))
Base.:*(a::GlobalSphericalPosition, b::Real) = GlobalSphericalPosition(a.r * b, a.theta, a.phi)
Base.:*(a::GlobalSphericalPosition, b::GlobalCartesianPosition) = GlobalCartesianPosition(a) * b
Base.:*(a::GlobalSphericalPosition, b::GlobalSphericalPosition) = GlobalCartesianPosition(a) * GlobalCartesianPosition(b)
Base.:*(a::GlobalCartesianVelocity, b::GlobalCartesianVelocity) = norm(_get(a) .* _get(b))
Base.:*(a::LocalCartesianVelocity, b::LocalCartesianVelocity) = norm(_get(a) .* _get(b))
Base.:*(a::GlobalSphericalVelocity, b::Real) = GlobalSphericalVelocity(a.r * b, a.theta, a.phi)


Base.:/(a::AbstractCoordinate, b::Real) = a * (1/b)
Base.:/(a::AbstractCoordinate, b::Tuple) = a * (1 ./ b)
Base.:/(a::AbstractCoordinate, b::AbstractVector) = a * (1 ./ b)


Base.eltype(a::T) where {T<:AbstractCoordinate} = typeof(_get(a)[1])


Base.isapprox(a::T, b::T; kwargs...) where {T<:AbstractCoordinate} = isapprox(_get(a), _get(b); kwargs...)


Base.isequal(a::T, b::T; kwargs...) where {T<:AbstractCoordinate} = isapprox(_get(a), _get(b); kwargs...)


Base.length(::T) where {T<:AbstractCoordinate} = 3


Base.size(::T) where {T<:AbstractCoordinate} = (3,)


Base.vec(a::T) where {T<:AbstractCartesianCoordinates} = [a.x, a.y, a.z]
Base.vec(a::T) where {T<:AbstractSphericalCoordinates} = [a.r, a.theta, a.phi]

LinearAlgebra.norm(a::T) where {T<:AbstractCartesianCoordinates} = sqrt(a.x^2 + a.y^2 + a.z^2)
LinearAlgebra.norm(x::GlobalSphericalPosition) = x.r


LinearAlgebra.dot(a::T, b::S) where {T<:AbstractCoordinate, S<:AbstractCoordinate} = a * b


#::. EXPORTS
export 
    GlobalCartesianPosition, 
    LocalCartesianPosition, 
    GlobalSphericalPosition, 
    GlobalCartesianVelocity, 
    LocalCartesianVelocity, 
    GlobalSphericalVelocity,

    azimuth, elevation, speed, zenith