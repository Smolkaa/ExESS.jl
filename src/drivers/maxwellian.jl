############################################################################################
#::. STRUCTS
############################################################################################
"""
    MBAzimuthDistribution(T, m)

Custom struct defining a (3D) Maxwell-Boltzmann azimuth angle distribution, based on the 
temperature `T` in (K) and the mass `m` in (kg).

# Arguments
- `T::Real`: Temperature (K)
- `m::Real`: Mass (kg)

# Notes
- The lower and upper bounds of the distribution are `-pi` and `pi`, respectively.
- Defined Methods: `rand`, `cdf`, `pdf`
"""
struct MBAzimuthDistribution{S<:AbstractFloat} <: AbstractDistribution; T::S; m::S; end
MBAzimuthDistribution(T::Real, m::Real) = MBAzimuthDistribution(promote(T, m)...)
MBAzimuthDistribution(T::Integer, m::Integer) = MBAzimuthDistribution(promote(T, m, 1.0)[1:2]...)


"""
    MBFluxAzimuthDistribution(T, m)

Custom struct defining a (3D) Maxwell-Boltzmann flux azimuth angle distribution, based on
the temperature `T` in (K) and the mass `m` in (kg).

# Arguments
- `T::Real`: Temperature (K)
- `m::Real`: Mass (kg)

# Notes
- The lower and upper bounds of the distribution are `-pi` and `pi`, respectively.
- Defined Methods: `rand`, `cdf`, `pdf`
"""
struct MBFluxAzimuthDistribution{S<:AbstractFloat} <: AbstractDistribution; T::S; m::S; end
MBFluxAzimuthDistribution(T::Real, m::Real) = MBFluxAzimuthDistribution(promote(T, m)...)
MBFluxAzimuthDistribution(T::Integer, m::Integer) = MBFluxAzimuthDistribution(promote(T, m, 1.0)[1:2]...)


"""
    MBElevationDistribution(T, m)

Custom struct defining a (3D) Maxwell-Boltzmann elevation angle distribution, based on the
temperature `T` in (K) and the mass `m` in (kg).

# Arguments
- `T::Real`: Temperature (K)
- `m::Real`: Mass (kg)

# Notes
- The lower and upper bounds of the distribution are `-pi/2` and `pi/2`, respectively.
- Defined Methods: `rand`, `cdf`, `pdf`
"""
struct MBElevationDistribution{S<:AbstractFloat} <: AbstractDistribution; T::S; m::S; end
MBElevationDistribution(T::Real, m::Real) = MBElevationDistribution(promote(T, m)...)
MBElevationDistribution(T::Integer, m::Integer) = MBElevationDistribution(promote(T, m, 1.0)[1:2]...)


"""
    MBFluxElevationDistribution(T, m)

Custom struct defining a (3D) Maxwell-Boltzmann flux elevation angle distribution, based on
the temperature `T` in (K) and the mass `m` in (kg).

# Arguments
- `T::Real`: Temperature (K)
- `m::Real`: Mass (kg)

# Notes
- The lower and upper bounds of the distribution are `0` and `pi/2`, respectively.
- Defined Methods: `rand`, `cdf`, `pdf`
"""
struct MBFluxElevationDistribution{S<:AbstractFloat} <: AbstractDistribution; T::S; m::S; end
MBFluxElevationDistribution(T::Real, m::Real) = MBFluxElevationDistribution(promote(T, m)...)
MBFluxElevationDistribution(T::Integer, m::Integer) = MBFluxElevationDistribution(promote(T, m, 1.0)[1:2]...)


"""
    MBSpeedDistribution(T, m)

Custom struct defining a (3D) Maxwell-Boltzmann speed distribution, based on the
temperature `T` in (K) and the mass `m` in (kg).

# Arguments
- `T::Real`: Temperature (K)
- `m::Real`: Mass (kg)

# Notes
- The lower and upper bounds of the distribution are `0` and `Inf`, respectively.
- Defined Methods: `rand`, `cdf`, `pdf`
"""
struct MBSpeedDistribution{S<:AbstractFloat} <: AbstractDistribution; T::S; m::S; end
MBSpeedDistribution(T::Real, m::Real) = MBSpeedDistribution(promote(T, m)...)
MBSpeedDistribution(T::Integer, m::Integer) = MBSpeedDistribution(promote(T, m, 1.0)[1:2]...)


"""
    MBFluxSpeedDistribution(T, m)

Custom struct defining a (3D) Maxwell-Boltzmann flux speed distribution, based on the
temperature `T` in (K) and the mass `m` in (kg).

# Arguments
- `T::Real`: Temperature (K)
- `m::Real`: Mass (kg)

# Notes
- The lower and upper bounds of the distribution are `0` and `Inf`, respectively.
- Defined Methods: `rand`, `cdf`, `pdf`
"""
struct MBFluxSpeedDistribution{S<:AbstractFloat} <: AbstractDistribution; T::S; m::S; end
MBFluxSpeedDistribution(T::Real, m::Real) = MBFluxSpeedDistribution(promote(T, m)...)
MBFluxSpeedDistribution(T::Integer, m::Integer) = MBFluxSpeedDistribution(promote(T, m, 1.0)[1:2]...)


"""
    MBVelocityDistribution(T, m)

Custom struct defining a (3D) Maxwell-Boltzmann velocity distribution, based on the
temperature `T` in (K) and the mass `m` in (kg).

# Arguments
- `T::Real`: Temperature (K)
- `m::Real`: Mass (kg)

# Notes
- The lower and upper bounds of the distribution are `(-Inf, -Inf, 0)` and
  `(Inf, Inf, Inf)`, respectively.
- Defined Methods: `rand`, `cdf`, `pdf`
"""
struct MBVelocityDistribution{S<:AbstractFloat} <: AbstractDistribution; T::S; m::S; end
MBVelocityDistribution(T::Real, m::Real) = MBVelocityDistribution(promote(T, m)...)
MBVelocityDistribution(T::Integer, m::Integer) = MBVelocityDistribution(promote(T, m, 1.0)[1:2]...)


"""
    MBFluxVelocityDistribution(T, m)

Custom struct defining a (3D) Maxwell-Boltzmann flux velocity distribution, based on the
temperature `T` in (K) and the mass `m` in (kg).

# Arguments
- `T::Real`: Temperature (K)
- `m::Real`: Mass (kg)

# Notes
- The lower and upper bounds of the distribution are `(-Inf, -Inf, 0)` and
  `(Inf, Inf, Inf)`, respectively.
- Defined Methods: `rand`, `cdf`, `pdf`
"""
struct MBFluxVelocityDistribution{S<:AbstractFloat} <: AbstractDistribution; T::S; m::S; end
MBFluxVelocityDistribution(T::Real, m::Real) = MBFluxVelocityDistribution(promote(T, m)...)
MBFluxVelocityDistribution(T::Integer, m::Integer) = MBFluxVelocityDistribution(promote(T, m, 1.0)[1:2]...)



############################################################################################
#::. FUNCTIONS

# internal union for simplified type handling
_TVLCV = Union{Tuple{Real, Real, Real}, AbstractVector{<:Real}, LocalCartesianVelocity}
_VDs = Union{MBVelocityDistribution, MBFluxVelocityDistribution}
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




Statistics.mean(d::MBElevationDistribution{S}) where {S<:AbstractFloat} = S((pi-2)/2)
Statistics.mean(d::MBFluxElevationDistribution{S}) where {S<:AbstractFloat} = S(pi/4)
function Statistics.mean(d::MBSpeedDistribution{S}) where {S<:AbstractFloat}
    return S(sqrt(8 * BOLTZMANN_CONSTANT * d.T / (pi * d.m)))
end
function Statistics.mean(d::MBFluxSpeedDistribution{S}) where {S<:AbstractFloat}
    return S(sqrt(9 * pi * BOLTZMANN_CONSTANT * d.T / (8 * d.m)))
end
# TODO: add mean value calculation for the other distributions



mode(::MBAzimuthDistribution{S}) where {S<:AbstractFloat} = zero(S)
mode(::MBFluxAzimuthDistribution{S}) where {S<:AbstractFloat} = zero(S)
# mode(::MBElevationDistribution{S}) where {S<:AbstractFloat} = zero(S)
# mode(::MBFluxElevationDistribution{S}) where {S<:AbstractFloat} = S(pi/4)
function mode(d::MBSpeedDistribution{S}) where {S<:AbstractFloat}
    return S(sqrt(2 * BOLTZMANN_CONSTANT * d.T / d.m))
end
function mode(d::MBFluxSpeedDistribution{S}) where {S<:AbstractFloat}
    return S(sqrt(3 * BOLTZMANN_CONSTANT * d.T / d.m))
end
# function mode(d::MBVelocityDistribution{S}) where {S<:AbstractFloat}
#     return (mode(MBSpeedDistribution(d.T, d.m)), zero(S), zero(S))
# end
# function mode(d::MBFluxVelocityDistribution{S}) where {S<:AbstractFloat}
#     v = mode(MBFluxSpeedDistribution(d.T, d.m)) * S(inv(sqrt(2)))
#     return (v, zero(S), v)
# end



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
    return S.((a*randn(), a*randn(), abs(a*randn())))
end
function Base.rand(d::MBFluxVelocityDistribution{S}) where {S<:AbstractFloat}
    a = sqrt(BOLTZMANN_CONSTANT * d.T / d.m)
    return S.(( a*randn(), a*randn(), a * sqrt(-2*log(1-rand())) ))
end

# overwriting additional `Base.rand` methods for multivariate maxwellians
Base.rand(S::Type{<:AbstractFloat}, d::_VDs) = S.(rand(d))
Base.rand(S::Type{<:LocalCartesianVelocity}, d::_VDs) = S(rand(d))
function Base.rand(S::Type{<:AbstractFloat}, d::_VDs, N::Integer) # perormance improvment
    SAMPLES = Vector{NTuple{3, S}}(undef, N)
    for i in 1:N; SAMPLES[i] = rand(S, d); end
    return SAMPLES
end
Base.rand(S::Type{<:LocalCartesianVelocity}, d::_VDs, N::Integer) = S.(rand(d, N))



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
