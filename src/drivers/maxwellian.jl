############################################################################################
#::. STRUCTS
############################################################################################
"""
    [1] MBAzimuthDistribution{S<:AbstractFloat}
    [2] MBAzimuthDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann azimuth angle distribution. Uses the 
temperature `T` in K and the mass `m` in kg as inputs.

The lower and upper bounds of the distribution are `-pi` and `pi`, respectively.

**Defined Methods**: `rand`, `cdf`, `pdf`
"""
struct MBAzimuthDistribution{S<:AbstractFloat} <: AbstractDistribution; T::S; m::S; end
MBAzimuthDistribution(T::Real, m::Real) = MBAzimuthDistribution(promote(T, m)...)
MBAzimuthDistribution(T::Integer, m::Integer) = MBAzimuthDistribution(promote(T, m, 1.0)[1:2]...)


"""
    [1] MBFluxAzimuthDistribution{S<:AbstractFloat}
    [2] MBFluxAzimuthDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann flux azimuth angle distribution. Uses the 
temperature `T` in K and the mass `m` in kg as inputs.

The lower and upper bounds of the distribution are `-pi` and `pi`, respectively.

**Defined Methods**: `rand`, `cdf`, `pdf`
"""
struct MBFluxAzimuthDistribution{S<:AbstractFloat} <: AbstractDistribution; T::S; m::S; end
MBFluxAzimuthDistribution(T::Real, m::Real) = MBFluxAzimuthDistribution(promote(T, m)...)
MBFluxAzimuthDistribution(T::Integer, m::Integer) = MBFluxAzimuthDistribution(promote(T, m, 1.0)[1:2]...)


"""
    [1] MBElevationDistribution{S<:AbstractFloat}
    [2] MBElevationDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann elevation angle distribution. Uses the 
temperature `T` in K and the mass `m` in kg as inputs.

The lower and upper bounds of the distribution are `-pi/2` and `pi/2`, respectively.

**Defined Methods**: `rand`, `cdf`, `pdf`
"""
struct MBElevationDistribution{S<:AbstractFloat} <: AbstractDistribution; T::S; m::S; end
MBElevationDistribution(T::Real, m::Real) = MBElevationDistribution(promote(T, m)...)
MBElevationDistribution(T::Integer, m::Integer) = MBElevationDistribution(promote(T, m, 1.0)[1:2]...)


"""
    [1] MBFluxElevationDistribution{S<:AbstractFloat}
    [2] MBFluxElevationDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann flux elevation angle distribution. Uses the 
temperature `T` in K and the mass `m` in kg as inputs.

The lower and upper bounds of the distribution are `0` and `pi/2`, respectively.

**Defined Methods**: `rand`, `cdf`, `pdf`
"""
struct MBFluxElevationDistribution{S<:AbstractFloat} <: AbstractDistribution; T::S; m::S; end
MBFluxElevationDistribution(T::Real, m::Real) = MBFluxElevationDistribution(promote(T, m)...)
MBFluxElevationDistribution(T::Integer, m::Integer) = MBFluxElevationDistribution(promote(T, m, 1.0)[1:2]...)


"""
    [1] MBSpeedDistribution{S<:AbstractFloat}
    [2] MBSpeedDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann speed distribution. Uses the temperature
`T` in K and the mass `m` in kg as inputs.

The lower and upper bounds of the distribution are `0` and `Inf`, respectively.

**Defined Methods**: `rand`, `cdf`, `pdf`
"""
struct MBSpeedDistribution{S<:AbstractFloat} <: AbstractDistribution; T::S; m::S; end
MBSpeedDistribution(T::Real, m::Real) = MBSpeedDistribution(promote(T, m)...)
MBSpeedDistribution(T::Integer, m::Integer) = MBSpeedDistribution(promote(T, m, 1.0)[1:2]...)


"""
    [1] MBFluxSpeedDistribution{S<:AbstractFloat}
    [2] MBFluxSpeedDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann flux speed distribution. Uses the 
temperature `T` in K and the mass `m` in kg as inputs.

The lower and upper bounds of the distribution are `0` and `Inf`, respectively.

**Defined Methods**: `rand`, `cdf`, `pdf`
"""
struct MBFluxSpeedDistribution{S<:AbstractFloat} <: AbstractDistribution; T::S; m::S; end
MBFluxSpeedDistribution(T::Real, m::Real) = MBFluxSpeedDistribution(promote(T, m)...)
MBFluxSpeedDistribution(T::Integer, m::Integer) = MBFluxSpeedDistribution(promote(T, m, 1.0)[1:2]...)


"""
    [1] MBVelocityDistribution{S<:AbstractFloat}
    [2] MBVelocityDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann velocity distribution. Uses the temperature
`T` in K and the mass `m` in kg as inputs.

The lower and upper bounds of the distribution are `(-Inf, -Inf, -Inf)` and 
`(Inf, Inf, Inf)`, respectively.

**Defined Methods**: `rand`, `cdf`, `pdf`
"""
struct MBVelocityDistribution{S<:AbstractFloat} <: AbstractDistribution; T::S; m::S; end
MBVelocityDistribution(T::Real, m::Real) = MBVelocityDistribution(promote(T, m)...)
MBVelocityDistribution(T::Integer, m::Integer) = MBVelocityDistribution(promote(T, m, 1.0)[1:2]...)


"""
    [1] MBFluxVelocityDistribution{S<:AbstractFloat}
    [2] MBFluxVelocityDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann flux velocity distribution. Uses the 
temperature `T` in K and the mass `m` in kg as inputs.

The lower and upper bounds of the distribution are `(-Inf, -Inf, -Inf)` and 
`(Inf, Inf, Inf)`, respectively.

**Defined Methods**: `rand`, `cdf`, `pdf`
"""
struct MBFluxVelocityDistribution{S<:AbstractFloat} <: AbstractDistribution; T::S; m::S; end
MBFluxVelocityDistribution(T::Real, m::Real) = MBFluxVelocityDistribution(promote(T, m)...)
MBFluxVelocityDistribution(T::Integer, m::Integer) = MBFluxVelocityDistribution(promote(T, m, 1.0)[1:2]...)



############################################################################################
#::. FUNCTIONS

# internal union for simplified type handling
_TVLCV = Union{Tuple{Real, Real, Real}, AbstractVector{<:Real}, LocalCartesianVelocity}
############################################################################################
function cdf(::MBAzimuthDistribution{S}, l::Real, u::Real) where {S<:AbstractFloat} 
    return S((u - l)/(2*pi))
end
cdf(d::MBAzimuthDistribution, u::Real) = cdf(d, -pi, u)

function cdf(::MBFluxAzimuthDistribution{S}, l::Real, u::Real) where {S<:AbstractFloat} 
    return S((u - l)/(2*pi))
end
cdf(d::MBFluxAzimuthDistribution, u::Real) = cdf(d, -pi, u)

function cdf(::MBElevationDistribution{S}, l::Real, u::Real) where {S<:AbstractFloat} 
    return S(0.5 * (sin(u)-sin(l)))
end
cdf(d::MBElevationDistribution, u::Real) = cdf(d, -pi/2, u)

function cdf(::MBFluxElevationDistribution{S}, l::Real, u::Real) where {S<:AbstractFloat} 
    return S(0.5 * (cos(2*l) - cos(2*u)))
end
cdf(d::MBFluxElevationDistribution, u::Real) = cdf(d, 0, u)

function cdf(d::MBSpeedDistribution{S}, l::Real, u::Real) where {S<:AbstractFloat}
    if d.T*d.m == 0; return zero(S); end
    a = d.m / (2 * BOLTZMANN_CONSTANT * d.T)
    c1 = (sqrt(pi) * erf(sqrt(a) * l)) / (4 * a^(3/2)) - (l * exp(-a*l^2)) / (2*a)
    c2 = (sqrt(pi) * erf(sqrt(a) * u)) / (4 * a^(3/2)) - (u * exp(-a*u^2)) / (2*a)
    return S(4 * pi * (a/pi)^(3/2) * (c2 - c1))
end
cdf(d::MBSpeedDistribution, u::Real) = cdf(d, 0, u)
function cdf(d::MBFluxSpeedDistribution{S}, l::Real, u::Real) where {S<:AbstractFloat}
    if d.T*d.m == 0; return zero(S); end
    a = d.m / (2 * BOLTZMANN_CONSTANT * d.T)
    c1 = (exp(-a*l^2) * (a*l^2 + 1)) / (2 * a^2)
    c2 = (exp(-a*u^2) * (a*u^2 + 1)) / (2 * a^2)
    return S(- 2 * a^2 * (c2 - c1))
end
cdf(d::MBFluxSpeedDistribution, u::Real) = cdf(d, 0, u)



function cdf(d::MBVelocityDistribution{S}, l::_TVLCV, u::_TVLCV) where {S<:AbstractFloat}
    return S(cdf(MBSpeedDistribution(d.T, d.m), speed(l), speed(u)) *
             cdf(MBElevationDistribution(d.T, d.m), elevation(l), elevation(u)) *
             cdf(MBAzimuthDistribution(d.T, d.m), azimuth(l), azimuth(u)))
end
function cdf(d::MBFluxVelocityDistribution{S}, l::_TVLCV, u::_TVLCV) where {S<:AbstractFloat}
    return S(cdf(MBFluxSpeedDistribution(d.T, d.m), speed(l), speed(u)) *
             cdf(MBFluxElevationDistribution(d.T, d.m), elevation(l), elevation(u)) *
             cdf(MBFluxAzimuthDistribution(d.T, d.m), azimuth(l), azimuth(u)))
end



function Statistics.mean(d::MBSpeedDistribution{S}) where {S<:AbstractFloat} 
    return S(sqrt(8 * BOLTZMANN_CONSTANT * d.T / (pi * d.m)))
end
function Statistics.mean(d::MBFluxSpeedDistribution{S}) where {S<:AbstractFloat} 
    return S(sqrt(9 * pi * BOLTZMANN_CONSTANT * d.T / (8 * d.m)))
end
# TODO: add mean value calculation for the other distributions



function mode(d::MBSpeedDistribution{S}) where {S<:AbstractFloat}
    return S(sqrt(2 * BOLTZMANN_CONSTANT * d.T / d.m))
end
function mode(d::MBFluxSpeedDistribution{S}) where {S<:AbstractFloat}
    return S(sqrt(3 * BOLTZMANN_CONSTANT * d.T / d.m))
end
# TODO: add mode value calculation for the other distributions



pdf(::MBAzimuthDistribution{S}, x::Real) where {S<:AbstractFloat} = S(1/(2*pi))
pdf(::MBFluxAzimuthDistribution{S}, x::Real) where {S<:AbstractFloat} = S(1/(2*pi))

pdf(::MBElevationDistribution{S}, x) where {S<:AbstractFloat} = S(0.5 * cos(x))
pdf(::MBFluxElevationDistribution{S}, x::Real) where {S<:AbstractFloat} = S(2*sin(x)*cos(x))

function pdf(d::MBSpeedDistribution{S}, x::Real) where {S<:AbstractFloat}
    if d.T == 0; return zero(S); end
    a = d.m / (2 * BOLTZMANN_CONSTANT * d.T)
    return S((a/pi)^(3//2) * 4*x^2*pi * exp(- a*x^2))
end
function pdf(d::MBFluxSpeedDistribution{S}, x::Real) where {S<:AbstractFloat}
    if d.T == 0; return zero(S); end
    a = d.m / (2 * BOLTZMANN_CONSTANT * d.T)
    return S(a^2 * x^3 * exp(- a*x^2))
end

function pdf(d::MBVelocityDistribution{S}, x::_TVLCV) where {S<:AbstractFloat}
    return S(pdf(MBSpeedDistribution(d.T, d.m), speed(x)) *
             pdf(MBElevationDistribution(d.T, d.m), elevation(x)) *
             pdf(MBAzimuthDistribution(d.T, d.m), azimuth(x)))
end
function pdf(d::MBFluxVelocityDistribution{S}, x::_TVLCV) where {S<:AbstractFloat}
    return S(pdf(MBFluxSpeedDistribution(d.T, d.m), speed(x)) *
             pdf(MBFluxElevationDistribution(d.T, d.m), elevation(x)) *
             pdf(MBFluxAzimuthDistribution(d.T, d.m), azimuth(x)))
end



Base.rand(::MBAzimuthDistribution{S}) where {S<:AbstractFloat} = S(rand()*2pi - pi)
Base.rand(::MBFluxAzimuthDistribution{S}) where {S<:AbstractFloat} = S(rand()*2pi - pi)

Base.rand(::MBElevationDistribution{S}) where {S<:AbstractFloat} = S(asin(2*rand() - 1))
Base.rand(::MBFluxElevationDistribution{S}) where {S<:AbstractFloat} = S(acos(sqrt(1-rand())))

function Base.rand(d::MBSpeedDistribution{S}) where {S<:AbstractFloat}
    return S(norm(rand(MBVelocityDistribution(d.T, d.m))))
end
function Base.rand(d::MBFluxSpeedDistribution{S}) where {S<:AbstractFloat}
    return S(norm(rand(MBFluxVelocityDistribution(d.T, d.m))))
end

function Base.rand(d::MBVelocityDistribution{S}) where {S<:AbstractFloat}
    a = sqrt(BOLTZMANN_CONSTANT * d.T / d.m)
    return S.((a*randn(), a*randn(), a*randn()))
end
Base.rand(S::Type{<:AbstractFloat}, d::MBVelocityDistribution) = S.(rand(d))
function Base.rand(d::MBFluxVelocityDistribution{S}) where {S<:AbstractFloat}
    a = sqrt(BOLTZMANN_CONSTANT * d.T / d.m)
    return S.(( a*randn(), a*randn(), a * sqrt(-2*log(1-rand())) ))
end
Base.rand(S::Type{<:AbstractFloat}, d::MBFluxVelocityDistribution) = S.(rand(d))



function rms(d::MBSpeedDistribution{S}) where {S<:AbstractFloat}
    return S(sqrt(3 * BOLTZMANN_CONSTANT * d.T / d.m))
end
function rms(d::MBFluxSpeedDistribution{S}) where {S<:AbstractFloat}
    return S(sqrt(4 * BOLTZMANN_CONSTANT * d.T / d.m))
end
# TODO: add rms value calculation for the other distributions



############################################################################################
#::. EXPORTS
############################################################################################
export 
    MBAzimuthDistribution, 
    MBElevationDistribution, 
    MBSpeedDistribution, 
    MBVelocityDistribution,

    MBFluxAzimuthDistribution, 
    MBFluxElevationDistribution, 
    MBFluxSpeedDistribution, 
    MBFluxVelocityDistribution