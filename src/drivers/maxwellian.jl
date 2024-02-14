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
    [1] MBSpeedDistribution{S<:AbstractFloat}
    [2] MBSpeedDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann speed distribution. Uses the temperature
`T` in K and the mass `m` in kg as inputs.

The lower and upper bounds of the distribution are `0` and `Inf`, respectively.

**Defined Methods**: `rand`, `cdf`, `pdf`
"""
struct MBSpeedDistribution{S<:AbstractFloat}; T::S; m::S; end
MBSpeedDistribution(T::Real, m::Real) = MBSpeedDistribution(promote(T, m)...)
MBSpeedDistribution(T::Integer, m::Integer) = MBSpeedDistribution(promote(T, m, 1.0)[1:2]...)


"""
    [1] MBVelocityDistribution{S<:AbstractFloat}
    [2] MBVelocityDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann velocity distribution. Uses the temperature
`T` in K and the mass `m` in kg as inputs.

The lower and upper bounds of the distribution are `(-Inf, -Inf, -Inf)` and 
`(Inf, Inf, Inf)`, respectively.

**Defined Methods**: `rand`, `cdf`, `pdf`
"""
struct MBVelocityDistribution{S<:AbstractFloat}; T::S; m::S; end
MBVelocityDistribution(T::Real, m::Real) = MBVelocityDistribution(promote(T, m)...)
MBVelocityDistribution(T::Integer, m::Integer) = MBVelocityDistribution(promote(T, m, 1.0)[1:2]...)



############################################################################################
#::. FUNCTIONS
############################################################################################
Base.rand(::MBAzimuthDistribution{S}) where {S<:AbstractFloat} = S(rand()*2pi - pi)

Base.rand(::MBElevationDistribution{S}) where {S<:AbstractFloat} = S(asin(2*rand() - 1))

function Base.rand(d::MBSpeedDistribution{S}) where {S<:AbstractFloat}
    return S(norm(rand(MBVelocityDistribution(d.T, d.m))))
end

function Base.rand(d::MBVelocityDistribution{S}) where {S<:AbstractFloat}
    a = sqrt(BOLTZMANN_CONSTANT * d.T / d.m)
    return S.((a*randn(), a*randn(), a*randn()))
end


cdf(::MBAzimuthDistribution{S}, l, u) where {S<:AbstractFloat} = S((u - l)/(2*pi))
cdf(d::MBAzimuthDistribution, u) = cdf(d, -pi, u)

cdf(::MBElevationDistribution{S}, l, u) where {S<:AbstractFloat} = S(0.5 * (sin(u)-sin(l)))
cdf(d::MBElevationDistribution, u) = cdf(d, -pi/2, u)

function cdf(d::MBSpeedDistribution{S}, l, u) where {S<:AbstractFloat}
    if d.T*d.m == 0; return zero(S); end
    a = d.m / (2 * BOLTZMANN_CONSTANT * d.T)
    c1 = (sqrt(pi) * erf(sqrt(a) * l)) / (4 * a^(3/2)) - (l * exp(-a*l^2)) / (2*a)
    c2 = (sqrt(pi) * erf(sqrt(a) * u)) / (4 * a^(3/2)) - (u * exp(-a*u^2)) / (2*a)
    return S(4 * pi * (a/pi)^(3/2) * (c2 - c1))
end
cdf(d::MBSpeedDistribution, u) = cdf(d, 0, u)

function cdf(d::MBVelocityDistribution{S}, l, u) where {S<:AbstractFloat}
    return S(cdf(MBSpeedDistribution(d.T, d.m), speed(l), speed(u)) *
             cdf(MBElevationDistribution(d.T, d.m), elevation(l), elevation(u)) *
             cdf(MBAzimuthDistribution(d.T, d.m), azimuth(l), azimuth(u)))
end
# cdf(d::MBVelocityDistribution, u) = cdf(d, (0,0,0), v)


pdf(::MBAzimuthDistribution{S}, x) where {S<:AbstractFloat} = S(1/(2*pi))

pdf(::MBElevationDistribution{S}, x) where {S<:AbstractFloat} = S(0.5 * cos(x))

function pdf(d::MBSpeedDistribution{S}, x) where {S<:AbstractFloat}
    if d.T == 0; return zero(S); end
    a = d.m / (2 * BOLTZMANN_CONSTANT * d.T)
    return S((a/pi)^(3//2) * 4*v^2*pi * exp(- a*x^2))
end

function pdf(d::MBVelocityDistribution{S}, x) where {S<:AbstractFloat}
    return S(pdf(MBSpeedDistribution(d.T, d.m), speed(x)) *
             pdf(MBElevationDistribution(d.T, d.m), elevation(x)) *
             pdf(MBAzimuthDistribution(d.T, d.m), azimuht(x)))
end


############################################################################################
#::. EXPORTS
############################################################################################
export 
    MBAzimuthDistribution, 
    MBElevationDistribution, 
    MBSpeedDistribution, 
    MBVelocityDistribution