############################################################################################
#::. STRUCTS
############################################################################################
abstract type AbstractEVector end
abstract type AbstractPosition{T<:AbstractFloat} <: AbstractEVector end
abstract type AbstractVelocity{T<:AbstractFloat} <: AbstractEVector end

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

struct LocalSphericalVelocity{T<:AbstractFloat} <: AbstractVelocity{T}
    r::T
    theta::T
    phi::T
end
function LocalSphericalVelocity(r::Real, theta::Real, phi::Real)
    LocalSphericalVelocity(promote(r, theta, phi)...)
end
function LocalSphericalVelocity(r::Integer, theta::Integer, phi::Integer)
    LocalSphericalVelocity(promote(r, theta, phi, 1.0)[1:3]...)
end


############################################################################################
#::. INTERNAL UNIONS
############################################################################################
const AbstractGlobalVector = Union{
    GlobalCartesianPosition,
    GlobalSphericalPosition,
    GlobalCartesianVelocity,}
const AbstractLocalVector = Union{
    LocalCartesianPosition,
    LocalCartesianVelocity,
    LocalSphericalVelocity}

const AbstractGlobalPosition = Union{GlobalCartesianPosition, GlobalSphericalPosition}
const AbstractGlobalVelocity = Union{GlobalCartesianVelocity}

const AbstractLocalPosition = Union{LocalCartesianPosition}
const AbstractLocalVelocity = Union{LocalCartesianVelocity, LocalSphericalVelocity}

const AbstractCartesianVector = Union{
    GlobalCartesianPosition,
    LocalCartesianPosition,
    GlobalCartesianVelocity,
    LocalCartesianVelocity,}
const AbstractSphericalVector = Union{
    GlobalSphericalPosition,
    LocalSphericalVelocity}


############################################################################################
#::. FUNCTIONS
############################################################################################
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
function GlobalCartesianVelocity(x::AbstractPosition, v::LocalSphericalVelocity)
    return GlobalCartesianVelocity(x, LocalCartesianVelocity(v))
end
function GlobalCartesianVelocity(x::GlobalSphericalPosition, v::LocalCartesianVelocity)
    return GlobalCartesianVelocity(
      - v.x*sin(x.theta) - v.y*cos(x.theta)*sin(x.phi) + v.z*cos(x.theta)*cos(x.phi),
        v.x*cos(x.theta) - v.y*sin(x.theta)*sin(x.phi) + v.z*sin(x.theta)*cos(x.phi),
                           v.y*             cos(x.phi) + v.z*             sin(x.phi),
    )
end



"""
    [1] LocalCartesianVelocity(x::Real, y::Real, z::Real)
    [2] LocalCartesianVelocity(V::Tuple)
    [3] LocalCartesianVelocity(V::AbstractVector)
    [4] LocalCartesianVelocity(v::LocalCartesianVelocity)
    [5] LocalCartesianVelocity(v::LocalSphericalVelocity)
    [6] LocalCartesianVelocity(x::AbstractGlobalPosition, v::AbstractVelocity)

Three dimensional velocity vector in a local cartesian coordinate system [1].
Converts tuples, vector or other velocity type into `LocalCartesianVelocity` [2,3,4,5,6].

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
function LocalCartesianVelocity(v::LocalSphericalVelocity{T}) where {T<:AbstractFloat}
    return LocalCartesianVelocity{T}(
        v.r * cos(v.theta) * cos(v.phi),
        v.r * sin(v.theta) * cos(v.phi),
        v.r * sin(v.phi),
    )
end
function LocalCartesianVelocity(::AbstractPosition, v::LocalSphericalVelocity)
    return LocalCartesianVelocity(v)
end




"""
    [1] LocalSphericalVelocity(r::Real, theta::Real, phi::Real)
    [2] LocalSphericalVelocity(V::Tuple)
    [3] LocalSphericalVelocity(V::AbstractVector)
    [4] LocalSphericalVelocity(v::LocalCartesianlVelocity)
    [5] LocalSphericalVelocity(v::LocalSphericalVelocity)
    [6] LocalSphericalVelocity(x::AbstractPosition, v::AbstractVelocity)

Three dimensional velocity vector in a local spherical coordinate system [1].
Converts tuples, vector or other velocity type into `LocalSphericalVelocity` [2,3,4,5,6].

In planetary contexts, local spherical coordinates are defined as follows:
- the r-axis is along the velocity vector, i.e. is the speed of the velocity vector
- the theta-angle describes the azimuth angle (positive longitudes, east), between the
  local x axis and the projection of the velocity vector onto the x-y plane
- the phi-angle describes the elevation angle (positive latitudes, north), between the
  velocity vector and the x-y plane
"""
LocalSphericalVelocity(V::Tuple) = LocalSphericalVelocity(V...)
LocalSphericalVelocity(V::AbstractVector) = LocalSphericalVelocity(V...)
LocalSphericalVelocity(v::LocalSphericalVelocity) = v
LocalSphericalVelocity(::AbstractPosition, v::LocalSphericalVelocity) = v
function LocalSphericalVelocity(x::AbstractPosition, v::AbstractVelocity)
    return LocalSphericalVelocity(LocalCartesianVelocity(x, v))
end
function LocalSphericalVelocity(v::LocalCartesianVelocity{T}) where {T<:AbstractFloat}
    r     = norm(v)
    theta = atan(v.y / v.x)
    phi   = asin(v.z / r)
    if isnan(theta); theta = 0;    end
    if v.x < 0;      theta += pi;  end
    if theta > pi;   theta -= 2pi; end
    return LocalSphericalVelocity{T}(r, theta, phi)
end


############################################################################################
#::. UTILITY FUNCTIONS
############################################################################################
_get(a::T) where {T<:AbstractCartesianVector} = (a.x, a.y, a.z)
_get(a::T) where {T<:AbstractSphericalVector} = (a.r, a.theta, a.phi)
_get(a::Tuple) = a
_get(a::AbstractVector) = Tuple(a)

_getx(a::T) where {T<:AbstractCartesianVector} = a.x
_getx(a::Union{Tuple, AbstractVector}) = a[1]

_gety(a::T) where {T<:AbstractCartesianVector} = a.y
_gety(a::Union{Tuple, AbstractVector}) = a[2]

_getz(a::T) where {T<:AbstractCartesianVector} = a.z
_getz(a::Union{Tuple, AbstractVector}) = a[3]

_getr(a::T) where {T<:AbstractSphericalVector} = a.r
_getr(a::Union{Tuple, AbstractVector}) = a[1]

_gettheta(a::T) where {T<:AbstractSphericalVector} = a.theta
_gettheta(a::Union{Tuple, AbstractVector}) = a[2]

_getphi(a::T) where {T<:AbstractSphericalVector} = a.phi
_getphi(a::Union{Tuple, AbstractVector}) = a[3]


"""
    [1] azimuth([S::Type{<:AbstractFloat}], x::Tuple)
    [2] azimuth([S::Type{<:AbstractFloat}], x::AbstractVector)
    [3] azimuth([S::Type{<:AbstractFloat}], x::AbstractEVector)

Calculate the azimuth angle (in [rad]) of the vector `x`.
"""
azimuth(x::Union{Tuple, AbstractVector}) = sgn(x[2]) * acos( x[1] / sqrt(x[1]^2 + x[2]^2) )
azimuth(x::AbstractCartesianVector) = sgn(x.y) * acos( x.x / sqrt(x.x^2 + x.y^2) )
azimuth(x::AbstractSphericalVector) = x.theta
azimuth(S::Type{<:AbstractFloat}, args...) = S(azimuth(args...))


"""
    [1] elevation([S::Type{<:AbstractFloat}], x::Tuple)
    [2] elevation([S::Type{<:AbstractFloat}], x::AbstractVector)
    [3] elevation([S::Type{<:AbstractFloat}], x::AbstractEVector)

Calculate the elevation angle (in [rad]) of the cartesian vector `x`.
"""
elevation(x::Union{Tuple, AbstractVector}) = atan(x[3] / sqrt(x[1]^2 + x[2]^2))
elevation(x::AbstractCartesianVector) = atan(x.z / sqrt(x.x^2 + x.y^2))
elevation(x::AbstractSphericalVector) = x.phi
elevation(S::Type{<:AbstractFloat}, args...) = S(elevation(args...))


"""
    [1] speed([S::Type], v::Tuple)
    [2] speed([S::Type], v::AbstractVector)
    [3] speed([S::Type], v::LocalCartesianVelocity)
    [4] speed([S::Type], v::GobalCartesianVelocity)
    [5] speed([S::Type], v::LocalSphericalVelocity)

Calculate the speed (in [m s-1]) of the velocity vector `v` (in [m s-1]).
"""
speed(v::Union{Tuple, AbstractVector, LocalCartesianVelocity, GlobalCartesianVelocity}) = norm(v)
speed(v::LocalSphericalVelocity) = v.r
speed(S::Type{<:AbstractFloat}, args...) = S(speed(args...))


"""
    [1] zenith([S::Type{<:Real}], c::Tuple)
    [2] zenith([S::Type{<:Real}], c::AbstractVector)
    [3] zenith([S::Type{<:Real}], c::AbstractEVector)

Calculate the zenith angle (in [rad]) of the the vector `x`.

**Notes**
- the zenith angle is the same-sign `pi/2`-inversion of the elevation angle
"""
zenith(x::Union{Tuple, AbstractVector}) = atan(sqrt(x[1]^2 + x[2]^2) / x[3])
zenith(x::AbstractCartesianVector) = atan(sqrt(x.x^2 + x.y^2) / x.z)
zenith(x::AbstractSphericalVector) = pi/2 - x.phi
zenith(S::Type{<:AbstractFloat}, args...) = S(zenith(args...))



############################################################################################
#::. EXTENSIONS
############################################################################################
Base.:+(a::AbstractEVector) = a
Base.:+(a::T, b::Tuple) where {T<:AbstractEVector} = T(_get(a) .+ b ...)
Base.:+(a::T, b::AbstractVector) where {T<:AbstractEVector} = T(vec(a) .+ b ...)
Base.:+(a::Tuple, b::AbstractEVector) = b + a
Base.:+(a::AbstractVector, b::AbstractEVector) = b + a
Base.:+(a::T, b::T) where {T<:AbstractCartesianVector} = T(a.x + b.x, a.y + b.y, a.z + b.z)
Base.:+(a::GlobalCartesianPosition, b::GlobalSphericalPosition) = a + GlobalCartesianPosition(b)
Base.:+(a::GlobalSphericalPosition, b::GlobalSphericalPosition) = GlobalSphericalPosition(GlobalCartesianPosition(a) + GlobalCartesianPosition(b))
Base.:+(a::GlobalSphericalPosition, b::GlobalCartesianPosition) = GlobalSphericalPosition(GlobalCartesianPosition(a) + b)
Base.:+(a::LocalCartesianVelocity, b::LocalSphericalVelocity) = a + LocalCartesianVelocity(b)
Base.:+(a::LocalSphericalVelocity, b::LocalSphericalVelocity) = LocalSphericalVelocity(LocalCartesianVelocity(a) + LocalCartesianVelocity(b))
Base.:+(a::LocalSphericalVelocity, b::LocalCartesianVelocity) = LocalSphericalVelocity(LocalCartesianVelocity(a) + b)


Base.:-(a::T) where {T<:AbstractEVector} = T( .- _get(a) ...)
Base.:-(a::AbstractEVector, b::Any) = a + (-b)
Base.:-(a::AbstractEVector, b::Tuple) = a + (.-b)
Base.:-(a::AbstractVector, b::AbstractEVector) = b + (-a)
Base.:-(a::Tuple, b::AbstractEVector) = b + (.-a)


Base.:*(a::T, b::Real) where {T<:AbstractEVector} = T(_get(a) .* b ... )
Base.:*(a::T, b::Tuple) where {T<:AbstractEVector} = T(_get(a) .* b ... )
Base.:*(a::T, b::AbstractVector) where {T<:AbstractEVector} = T(vec(a) .* b ... )
Base.:*(a::Real, b::AbstractEVector) = b * a
Base.:*(a::Tuple, b::AbstractEVector) = b * a
Base.:*(a::AbstractVector, b::AbstractEVector) = b * a
Base.:*(a::GlobalCartesianPosition, b::GlobalCartesianPosition) = norm(_get(a) .* _get(b))
Base.:*(a::GlobalCartesianPosition, b::GlobalSphericalPosition) = a * GlobalCartesianPosition(b)
Base.:*(a::LocalCartesianPosition, b::LocalCartesianPosition) = norm(_get(a) .* _get(b))
Base.:*(a::GlobalSphericalPosition, b::Real) = GlobalSphericalPosition(a.r * b, a.theta, a.phi)
Base.:*(a::GlobalSphericalPosition, b::GlobalCartesianPosition) = GlobalCartesianPosition(a) * b
Base.:*(a::GlobalSphericalPosition, b::GlobalSphericalPosition) = GlobalCartesianPosition(a) * GlobalCartesianPosition(b)
Base.:*(a::GlobalCartesianVelocity, b::GlobalCartesianVelocity) = norm(_get(a) .* _get(b))
Base.:*(a::LocalCartesianVelocity, b::LocalCartesianVelocity) = norm(_get(a) .* _get(b))
Base.:*(a::LocalCartesianVelocity, b::LocalSphericalVelocity) = a * LocalCartesianVelocity(b)
Base.:*(a::LocalSphericalVelocity, b::Real) = LocalSphericalVelocity(a.r * b, a.theta, a.phi)
Base.:*(a::LocalSphericalVelocity, b::LocalCartesianVelocity) = LocalCartesianVelocity(a) * b
Base.:*(a::LocalSphericalVelocity, b::LocalSphericalVelocity) = LocalCartesianVelocity(a) * LocalCartesianVelocity(b)



Base.:/(a::AbstractEVector, b::Real) = a * (1/b)
Base.:/(a::AbstractEVector, b::Tuple) = a * (1 ./ b)
Base.:/(a::AbstractEVector, b::AbstractVector) = a * (1 ./ b)


Base.eltype(a::T) where {T<:AbstractEVector} = typeof(_get(a)[1])


Base.isapprox(a::T, b::T; kwargs...) where {T<:AbstractEVector} = isapprox(_get(a), _get(b); kwargs...)


Base.isequal(a::T, b::T; kwargs...) where {T<:AbstractEVector} = isapprox(_get(a), _get(b); kwargs...)


Base.length(::T) where {T<:AbstractEVector} = 3


Base.size(::T) where {T<:AbstractEVector} = (3,)


Base.vec(a::T) where {T<:AbstractCartesianVector} = [a.x, a.y, a.z]
Base.vec(a::T) where {T<:AbstractSphericalVector} = [a.r, a.theta, a.phi]


LinearAlgebra.dot(a::T, b::S) where {T<:AbstractEVector, S<:AbstractEVector} = a * b


LinearAlgebra.norm(a::T) where {T<:AbstractCartesianVector} = sqrt(a.x^2 + a.y^2 + a.z^2)
LinearAlgebra.norm(x::AbstractSphericalVector) = x.r


############################################################################################
#::. EXPORTS
############################################################################################
export
    GlobalCartesianPosition,
    LocalCartesianPosition,
    GlobalSphericalPosition,
    GlobalCartesianVelocity,
    LocalCartesianVelocity,
    LocalSphericalVelocity,

    azimuth, elevation, speed, zenith
