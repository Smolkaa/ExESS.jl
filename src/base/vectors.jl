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


############################################################################################
#::. INTERNAL UNIONS
############################################################################################
const AbstractGlobalVector = Union{
    GlobalCartesianPosition,
    GlobalSphericalPosition,
    GlobalCartesianVelocity,}
const AbstractLocalVector = Union{
    LocalCartesianPosition,
    LocalCartesianVelocity,}

const AbstractGlobalPosition = Union{GlobalCartesianPosition, GlobalSphericalPosition}
const AbstractGlobalVelocity = Union{GlobalCartesianVelocity}

const AbstractLocalPosition = Union{LocalCartesianPosition}
const AbstractLocalVelocity = Union{LocalCartesianVelocity}

const AbstractCartesianVector = Union{
    GlobalCartesianPosition,
    LocalCartesianPosition,
    GlobalCartesianVelocity,
    LocalCartesianVelocity,}
const AbstractSphericalVector = Union{
    GlobalSphericalPosition,}


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



############################################################################################
#::. INTERNAL UTILITY FUNCTIONS
############################################################################################
_getx(a::T) where {T<:AbstractCartesianVector} = a.x
_getx(a::GlobalSphericalPosition) = _getx(GlobalCartesianPosition(a))
_getx(a::Union{Tuple, AbstractVector}) = a[1]

_gety(a::T) where {T<:AbstractCartesianVector} = a.y
_gety(a::GlobalSphericalPosition) = _gety(GlobalCartesianPosition(a))
_gety(a::Union{Tuple, AbstractVector}) = a[2]

_getz(a::T) where {T<:AbstractCartesianVector} = a.z
_getz(a::GlobalSphericalPosition) = _getz(GlobalCartesianPosition(a))
_getz(a::Union{Tuple, AbstractVector}) = a[3]

_getr(a::T) where {T<:AbstractSphericalVector} = a.r
_getr(a::GlobalCartesianPosition) = _getr(GlobalSphericalPosition(a))
_getr(a::Union{Tuple, AbstractVector}) = a[1]

_gettheta(a::T) where {T<:AbstractSphericalVector} = a.theta
_gettheta(a::GlobalCartesianPosition) = _gettheta(GlobalSphericalPosition(a))
_gettheta(a::Union{Tuple, AbstractVector}) = a[2]

_getphi(a::T) where {T<:AbstractSphericalVector} = a.phi
_getphi(a::GlobalCartesianPosition) = _getphi(GlobalSphericalPosition(a))
_getphi(a::Union{Tuple, AbstractVector}) = a[3]



############################################################################################
#::. EXPORTED UTILITY FUNCTIONS
############################################################################################
"""
    azimuth([S], x)

Calculate the azimuth angle in (rad) of the tuple/vector `x`. The optional argument `S` can
be provided to cast the result to a specific type.
"""
azimuth(x::Union{Tuple, AbstractVector}) = sgn(x[2]) * acos( x[1] / sqrt(x[1]^2 + x[2]^2) )
azimuth(x::AbstractCartesianVector) = sgn(x.y) * acos( x.x / sqrt(x.x^2 + x.y^2) )
azimuth(x::AbstractSphericalVector) = x.theta
azimuth(S::Type{<:AbstractFloat}, args...) = S(azimuth(args...))


"""
    elevation([S], x)

Calculate the elevation angle in (rad) of the tuple/vector `x`. The optional argument `S`
can be provided to cast the result to a specific type.
"""
elevation(x::Union{Tuple, AbstractVector}) = atan(x[3] / sqrt(x[1]^2 + x[2]^2))
elevation(x::AbstractCartesianVector) = atan(x.z / sqrt(x.x^2 + x.y^2))
elevation(x::AbstractSphericalVector) = x.phi
elevation(S::Type{<:AbstractFloat}, args...) = S(elevation(args...))



"""
    rotate_x(x, a)

Rotates a cartesian tuple/vector `x` around the x-axis by the angle `a` in (rad).
"""
function rotate_x(x::Tuple, a::Real)
    return (x[1], x[2]*cos(a) - x[3]*sin(a), x[2]*sin(a) + x[3]*cos(a))
end
rotate_x(x::AbstractVector, a::Real) = [rotate_x(Tuple(x), a)...]
rotate_x(x::T, a::Real) where {T<:AbstractCartesianVector} = T(rotate_x(Tuple(x), a)...)



"""
    rotate_y(x, a)

Rotates a cartesian tuple/vector `x` around the y-axis by the angle `a` in (rad).
"""
function rotate_y(x::Tuple, a::Real)
    return (x[1]*cos(a) + x[3]*sin(a), x[2], -x[1]*sin(a) + x[3]*cos(a))
end
rotate_y(x::AbstractVector, a::Real) = [rotate_y(Tuple(x), a)...]
rotate_y(x::T, a::Real) where {T<:AbstractCartesianVector} = T(rotate_y(Tuple(x), a)...)



"""
    rotate_z(x, a)

Rotates a cartesian tuple/vector `x` around the z-axis by the angle `a` in (rad).
"""
function rotate_z(x::Tuple, a::Real)
    return (x[1]*cos(a) - x[2]*sin(a), x[1]*sin(a) + x[2]*cos(a), x[3])
end
rotate_z(x::AbstractVector, a::Real) = [rotate_z(Tuple(x), a)...]
rotate_z(x::T, a::Real) where {T<:AbstractCartesianVector} = T(rotate_z(Tuple(x), a)...)



"""
    speed([S], v)

Calculate the speed in (m s-1) of the velocity vector `v` in (m s-1). The optional argument
`S` can be provided to cast the result to a specific type.
"""
speed(v::Union{Tuple, AbstractVector, LocalCartesianVelocity, GlobalCartesianVelocity}) = norm(v)
speed(S::Type{<:AbstractFloat}, args...) = S(speed(args...))


"""
    zenith([S], x)

Calculate the zenith angle in (rad) of the the vector `x`. The optional argument `S` can
be provided to cast the result to a specific type.

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
Base.:+(a::T, b::Tuple) where {T<:AbstractEVector} = T(Tuple(a) .+ b ...)
Base.:+(a::T, b::AbstractVector) where {T<:AbstractEVector} = T(vec(a) .+ b ...)
Base.:+(a::Tuple, b::AbstractEVector) = b + a
Base.:+(a::AbstractVector, b::AbstractEVector) = b + a
Base.:+(a::T, b::T) where {T<:AbstractCartesianVector} = T(a.x + b.x, a.y + b.y, a.z + b.z)
Base.:+(a::GlobalCartesianPosition, b::GlobalSphericalPosition) = a + GlobalCartesianPosition(b)
Base.:+(a::GlobalSphericalPosition, b::GlobalSphericalPosition) = GlobalSphericalPosition(GlobalCartesianPosition(a) + GlobalCartesianPosition(b))
Base.:+(a::GlobalSphericalPosition, b::GlobalCartesianPosition) = GlobalSphericalPosition(GlobalCartesianPosition(a) + b)


Base.:-(a::T) where {T<:AbstractEVector} = T( .- Tuple(a) ...)
Base.:-(a::AbstractEVector, b::Any) = a + (-b)
Base.:-(a::AbstractEVector, b::Tuple) = a + (.-b)
Base.:-(a::AbstractVector, b::AbstractEVector) = b + (-a)
Base.:-(a::Tuple, b::AbstractEVector) = b + (.-a)


Base.:*(a::T, b::Real) where {T<:AbstractEVector} = T(Tuple(a) .* b ... )
Base.:*(a::T, b::Tuple) where {T<:AbstractEVector} = T(Tuple(a) .* b ... )
Base.:*(a::T, b::AbstractVector) where {T<:AbstractEVector} = T(vec(a) .* b ... )
Base.:*(a::Real, b::AbstractEVector) = b * a
Base.:*(a::Tuple, b::AbstractEVector) = b * a
Base.:*(a::AbstractVector, b::AbstractEVector) = b * a
Base.:*(a::GlobalCartesianPosition, b::GlobalCartesianPosition) = norm(Tuple(a) .* Tuple(b))
Base.:*(a::GlobalCartesianPosition, b::GlobalSphericalPosition) = a * GlobalCartesianPosition(b)
Base.:*(a::LocalCartesianPosition, b::LocalCartesianPosition) = norm(Tuple(a) .* Tuple(b))
Base.:*(a::GlobalSphericalPosition, b::Real) = GlobalSphericalPosition(a.r * b, a.theta, a.phi)
Base.:*(a::GlobalSphericalPosition, b::GlobalCartesianPosition) = GlobalCartesianPosition(a) * b
Base.:*(a::GlobalSphericalPosition, b::GlobalSphericalPosition) = GlobalCartesianPosition(a) * GlobalCartesianPosition(b)
Base.:*(a::GlobalCartesianVelocity, b::GlobalCartesianVelocity) = norm(Tuple(a) .* Tuple(b))
Base.:*(a::LocalCartesianVelocity, b::LocalCartesianVelocity) = norm(Tuple(a) .* Tuple(b))



Base.:/(a::AbstractEVector, b::Real) = a * (1/b)
Base.:/(a::AbstractEVector, b::Tuple) = a * (1 ./ b)
Base.:/(a::AbstractEVector, b::AbstractVector) = a * (1 ./ b)


Base.eltype(a::T) where {T<:AbstractEVector} = eltype(Tuple(a))


Base.isapprox(a::T, b::T; kwargs...) where {T<:AbstractEVector} = isapprox(Tuple(a), Tuple(b); kwargs...)
Base.isapprox(a::AbstractEVector, b::Tuple; kwargs...) = isapprox(Tuple(a), b; kwargs...)
Base.isapprox(a::Tuple, b::AbstractEVector; kwargs...) = isapprox(a, Tuple(b); kwargs...)
Base.isapprox(a::AbstractEVector, b::AbstractVector; kwargs...) = isapprox(vec(a), b; kwargs...)
Base.isapprox(a::AbstractVector, b::AbstractEVector; kwargs...) = isapprox(a, vec(b); kwargs...)


Base.isequal(a::T, b::T; kwargs...) where {T<:AbstractEVector} = isapprox(Tuple(a), Tuple(b); kwargs...)
Base.isequal(a::AbstractEVector, b::Tuple; kwargs...) = isequal(Tuple(a), b; kwargs...)
Base.isequal(a::Tuple, b::AbstractEVector; kwargs...) = isequal(a, Tuple(b); kwargs...)
Base.isequal(a::AbstractEVector, b::AbstractVector; kwargs...) = isequal(vec(a), b; kwargs...)
Base.isequal(a::AbstractVector, b::AbstractEVector; kwargs...) = isequal(a, vec(b); kwargs...)


Base.length(::T) where {T<:AbstractEVector} = 3


Base.size(::T) where {T<:AbstractEVector} = (3,)


Base.Tuple(a::T) where {T<:AbstractCartesianVector} = (a.x, a.y, a.z)
Base.Tuple(a::T) where {T<:AbstractSphericalVector} = (a.r, a.theta, a.phi)


Base.vec(a::T) where {T<:AbstractCartesianVector} = [a.x, a.y, a.z]
Base.vec(a::T) where {T<:AbstractSphericalVector} = [a.r, a.theta, a.phi]


LinearAlgebra.dot(a::T, b::S) where {T<:AbstractEVector, S<:AbstractEVector} = a * b


LinearAlgebra.norm(a::T) where {T<:AbstractCartesianVector} = sqrt(a.x^2 + a.y^2 + a.z^2)
LinearAlgebra.norm(x::AbstractSphericalVector) = x.r
# LinearAlgebra.normalize is automatically defined now.


############################################################################################
#::. EXPORTS
############################################################################################
export
    GlobalCartesianPosition,
    LocalCartesianPosition,
    GlobalSphericalPosition,
    GlobalCartesianVelocity,
    LocalCartesianVelocity,

    azimuth,
    elevation,
    rotate_x, rotate_y, rotate_z,
    speed,
    zenith
