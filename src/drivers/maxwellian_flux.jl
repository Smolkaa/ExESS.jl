############################################################################################
#::. STRUCTS
############################################################################################
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
    [1] MBFluxSpeedDistribution{S<:AbstractFloat}
    [2] MBFluxSpeedDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann flux speed distribution. Uses the 
temperature `T` in K and the mass `m` in kg as inputs.

The lower and upper bounds of the distribution are `0` and `Inf`, respectively.

**Defined Methods**: `rand`, `cdf`, `pdf`
"""
struct MBFluxSpeedDistribution{S<:AbstractFloat}; T::S; m::S; end
MBFluxSpeedDistribution(T::Real, m::Real) = MBFluxSpeedDistribution(promote(T, m)...)
MBFluxSpeedDistribution(T::Integer, m::Integer) = MBFluxSpeedDistribution(promote(T, m, 1.0)[1:2]...)


"""
    [1] MBFluxVelocityDistribution{S<:AbstractFloat}
    [2] MBFluxVelocityDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann flux velocity distribution. Uses the 
temperature `T` in K and the mass `m` in kg as inputs.

The lower and upper bounds of the distribution are `(-Inf, -Inf, -Inf)` and 
`(Inf, Inf, Inf)`, respectively.

**Defined Methods**: `rand`, `cdf`, `pdf`
"""
struct MBFluxVelocityDistribution{S<:AbstractFloat}; T::S; m::S; end
MBFluxVelocityDistribution(T::Real, m::Real) = MBFluxVelocityDistribution(promote(T, m)...)
MBFluxVelocityDistribution(T::Integer, m::Integer) = MBFluxVelocityDistribution(promote(T, m, 1.0)[1:2]...)



############################################################################################
#::. FUNCTIONS
############################################################################################
Base.rand(::MBFluxAzimuthDistribution{S}) where {S<:AbstractFloat} = S(rand()*2pi - pi)

Base.rand(::MBFluxElevationDistribution{S}) where {S<:AbstractFloat} = S(acos(sqrt(1 - rand())))

function Base.rand(d::MBFluxSpeedDistribution{S}) where {S<:AbstractFloat}
    return S(norm(rand(MBFluxVelocityDistribution(d.T, d.m))))
end

function Base.rand(d::MBFluxVelocityDistribution{S}) where {S<:AbstractFloat}
    a = sqrt(BOLTZMANN_CONSTANT * d.T / d.m)
    return S.(( a*randn(), a*randn(), a * sqrt(-2*log(1-rand())) ))
end


cdf(::MBFluxAzimuthDistribution{S}, l, u) where {S<:AbstractFloat} = S((u - l)/(2*pi))
cdf(d::MBFluxAzimuthDistribution, u) = cdf(d, -pi, u)

cdf(::MBFluxElevationDistribution{S}, l, u) where {S<:AbstractFloat} = S(0.5 * (cos(2*l) - cos(2*u)))
cdf(d::MBFluxElevationDistribution, u) = cdf(d, 0, u)

function cdf(d::MBFluxSpeedDistribution{S}, l, u) where {S<:AbstractFloat}
    if d.T*d.m == 0; return zero(S); end
    a = sd.m / (2 * BOLTZMANN_CONSTANT * sd.T)
    c1 = (exp(-a*l^2) * (a*l^2 + 1)) / (2 * a^2)
    c2 = (exp(-a*u^2) * (a*u^2 + 1)) / (2 * a^2)
    return S(- 2 * a^2 * (c2 - c1))
end
cdf(d::MBFluxSpeedDistribution, u) = cdf(d, 0, u)

function cdf(d::MBFluxVelocityDistribution{S}, l, u) where {S<:AbstractFloat}
    return S(cdf(MBFluxSpeedDistribution(d.T, d.m), speed(l), speed(u)) *
             cdf(MBFluxElevationDistribution(d.T, d.m), elevation(l), elevation(u)) *
             cdf(MBFluxAzimuthDistribution(d.T, d.m), azimuth(l), azimuth(u)))
end
# cdf(d::MBFluxVelocityDistribution, u) = cdf(d, (0,0,0), v)


pdf(::MBFluxAzimuthDistribution{S}, x) where {S<:AbstractFloat} = S(1/(2*pi))

pdf(::MBFluxElevationDistribution{S}, x) where {S<:AbstractFloat} = S(2 * sin(x) * cos(x))

function pdf(d::MBFluxSpeedDistribution{S}, x) where {S<:AbstractFloat}
    if d.T == 0; return zero(S); end
    a = sd.m / (2 * BOLTZMANN_CONSTANT * sd.T)
    return S(a^2 * x^3 * exp(- a*x^2))
end

function pdf(d::MBFluxVelocityDistribution{S}, x) where {S<:AbstractFloat}
    return S(pdf(MBFluxSpeedDistribution(d.T, d.m), speed(x)) *
             pdf(MBFluxElevationDistribution(d.T, d.m), elevation(x)) *
             pdf(MBFluxAzimuthDistribution(d.T, d.m), azimuht(x)))
end


############################################################################################
#::. EXPORTS
############################################################################################
export 
    MBFluxAzimuthDistribution, 
    MBFluxElevationDistribution, 
    MBFluxSpeedDistribution, 
    MBFluxVelocityDistribution