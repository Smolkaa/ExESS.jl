#::. STRUCTS
abstract type AbstractMBDistribution end
abstract type AbstractMBAzimuthDistribution{S<:AbstractFloat} <: AbstractMBDistribution end
abstract type AbstractMBElevationDistribution{S<:AbstractFloat} <: AbstractMBDistribution end
abstract type AbstractMBSpeedDistribution{S<:AbstractFloat} <: AbstractMBDistribution end
abstract type AbstractMBVelocityDistribution{S<:AbstractFloat} <: AbstractMBDistribution end


"""
    [1] MBAzimuthDistribution{S<:AbstractFloat}
    [2] MBAzimuthDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann azimuth angle distribution. Uses the temperature
`T` in K and the mass `m` in kg as inputs.
"""
struct MBAzimuthDistribution{S<:AbstractFloat} <: AbstractMBAzimuthDistribution{S}; T::S; m::S; end
MBAzimuthDistribution(T::Real, m::Real) = MBAzimuthDistribution(promote(T, m)...)
MBAzimuthDistribution(T::Integer, m::Integer) = MBAzimuthDistribution(promote(T, m, 1.0)[1:2]...)


"""
    [1] MBElevationDistribution{S<:AbstractFloat}
    [2] MBElevationDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann elevation angle distribution. Uses the temperature
`T` in K and the mass `m` in kg as inputs.
"""
struct MBElevationDistribution{S<:AbstractFloat} <: AbstractMBElevationDistribution{S}; T::S; m::S; end
MBElevationDistribution(T::Real, m::Real) = MBElevationDistribution(promote(T, m)...)
MBElevationDistribution(T::Integer, m::Integer) = MBElevationDistribution(promote(T, m, 1.0)[1:2]...)


"""
    [1] MBSpeedDistribution{S<:AbstractFloat}
    [2] MBSpeedDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann speed distribution. Uses the temperature
`T` in K and the mass `m` in kg as inputs.
"""
struct MBSpeedDistribution{S<:AbstractFloat} <: AbstractMBSpeedDistribution{S}; T::S; m::S; end
MBSpeedDistribution(T::Real, m::Real) = MBSpeedDistribution(promote(T, m)...)
MBSpeedDistribution(T::Integer, m::Integer) = MBSpeedDistribution(promote(T, m, 1.0)[1:2]...)


"""
    [1] MBVelocityDistribution{S<:AbstractFloat}
    [2] MBVelocityDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann velocity distribution. Uses the temperature
`T` in K and the mass `m` in kg as inputs.
"""
struct MBVelocityDistribution{S<:AbstractFloat} <: AbstractMBVelocityDistribution{S}; T::S; m::S; end
MBVelocityDistribution(T::Real, m::Real) = MBVelocityDistribution(promote(T, m)...)
MBVelocityDistribution(T::Integer, m::Integer) = MBVelocityDistribution(promote(T, m, 1.0)[1:2]...)


"""
    [1] MBFluxAzimuthDistribution{S<:AbstractFloat}
    [2] MBFluxAzimuthDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann flux azimuth angle distribution. Uses the 
temperature `T` in K and the mass `m` in kg as inputs.
"""
struct MBFluxAzimuthDistribution{S<:AbstractFloat} <: AbstractMBAzimuthDistribution{S}; T::S; m::S; end
MBFluxAzimuthDistribution(T::Real, m::Real) = MBFluxAzimuthDistribution(promote(T, m)...)
MBFluxAzimuthDistribution(T::Integer, m::Integer) = MBFluxAzimuthDistribution(promote(T, m, 1.0)[1:2]...)


"""
    [1] MBFluxElevationDistribution{S<:AbstractFloat}
    [2] MBFluxElevationDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann flux elevation angle distribution. Uses the 
temperature `T` in K and the mass `m` in kg as inputs.
"""
struct MBFluxElevationDistribution{S<:AbstractFloat} <: AbstractMBElevationDistribution{S}; T::S; m::S; end
MBFluxElevationDistribution(T::Real, m::Real) = MBFluxElevationDistribution(promote(T, m)...)
MBFluxElevationDistribution(T::Integer, m::Integer) = MBFluxElevationDistribution(promote(T, m, 1.0)[1:2]...)


"""
    [1] MBFluxSpeedDistribution{S<:AbstractFloat}
    [2] MBFluxSpeedDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann flux speed distribution. Uses the 
temperature `T` in K and the mass `m` in kg as inputs.
"""
struct MBFluxSpeedDistribution{S<:AbstractFloat} <: AbstractMBSpeedDistribution{S}; T::S; m::S; end
MBFluxSpeedDistribution(T::Real, m::Real) = MBFluxSpeedDistribution(promote(T, m)...)
MBFluxSpeedDistribution(T::Integer, m::Integer) = MBFluxSpeedDistribution(promote(T, m, 1.0)[1:2]...)


"""
    [1] MBFluxVelocityDistribution{S<:AbstractFloat}
    [2] MBFluxVelocityDistribution(T::Real, m::Real)

Custom struct defining a (3D) Maxwell-Boltzmann flux velocity distribution. Uses the 
temperature `T` in K and the mass `m` in kg as inputs.
"""
struct MBFluxVelocityDistribution{S<:AbstractFloat} <: AbstractMBVelocityDistribution{S}; T::S; m::S; end
MBFluxVelocityDistribution(T::Real, m::Real) = MBFluxVelocityDistribution(promote(T, m)...)
MBFluxVelocityDistribution(T::Integer, m::Integer) = MBFluxVelocityDistribution(promote(T, m, 1.0)[1:2]...)



#::. FUNCTIONS
"""
    [1] cdf([S::Type{<:AbstractFloat}], ad::AbstractMBAzimuthDistribution, theta1::Real, 
            [theta2::Real])
    [2] cdf([S::Type{<:AbstractFloat}], ad::AbstractMBAzimuthDistribution, 
            theta::Tuple{Real, Real})

Computes the cumulative distribution function of the azimuth angle of the azimuth
distribution `ad`, between the azimuth angles given in `theta` (in [rad]). Should 
only one angle be given, the lower boundary is assumed to be lower boundary of the
respective domain, `-pi`.
"""
function cdf(::AbstractMBAzimuthDistribution{S}, 
             theta1::Real, theta2::Real) where {S<:AbstractFloat}
    return S((theta2 - theta1)/(2*pi))
end
function cdf(ad::AbstractMBAzimuthDistribution, theta::Tuple{Real, Real})
    return cdf(ad, theta[1], theta[2])
end
function cdf(ad::AbstractMBAzimuthDistribution{S}, theta::Real) where {S<:AbstractFloat}
    return cdf(ad, S(-pi), theta)
end
cdf(S::Type{<:AbstractFloat}, ad::AbstractMBAzimuthDistribution, args...) = S(cdf(ad, args...))


"""
    [1] cdf([S::Type{<:AbstractFloat}], ed::AbstractMBElevationDistribution, phi1::Real, [phi2::Real])
    [2] cdf([S::Type{<:AbstractFloat}], ed::AbstractMBElevationDistribution, phi{Real, Real})

Computes the cumulative distribution function of the elevation angle of the elevation 
distribution `ed`, between the elevation angles given in `phi` (in [rad]). Should only one 
angle be given, the lower boundary is assumed to be lower boundary of the respective domain. 

**Notes**
- Domain of `MBVelocityDistribution` is `[-pi/2, pi/2]`, thus the lower boundary is `-pi/2`.
- Domain of `MBFluxVelocityDistribution` is `[0, pi/2]`, thus the lower boundary is `0`.
"""
function cdf(::MBElevationDistribution{S}, 
             phi1::Real, phi2::Real) where {S<:AbstractFloat}
    return S(0.5 * (sin(phi2) - sin(phi1)))
end
cdf(vd::MBElevationDistribution, phi::Real) = cdf(vd, -pi/2, phi)
function cdf(::MBFluxElevationDistribution{S}, 
             phi1::Real, phi2::Real) where {S<:AbstractFloat}
    return S(0.5 * (cos(2*phi1) - cos(2*phi2)))
end
cdf(vd::MBFluxElevationDistribution, phi::Real) = cdf(vd, 0, phi)
function cdf(ed::AbstractMBElevationDistribution, phi::Tuple{Real, Real})
    return cdf(ed, phi[1], phi[2])
end
cdf(S::Type{<:AbstractFloat}, ed::AbstractMBElevationDistribution, args...) = S(cdf(ed, args...))


"""
    [1] cdf([S::Type{<:AbstractFloat}], vd::AbstractMBSpeedDistribution, v1::Real, [v2::Real])
    [2] cdf([S::Type{<:AbstractFloat}], vd::AbstractMBSpeedDistribution, v::Tuple{Real, Real})

Computes the cumulative distribution function of the speed distribution `sd` between the 
upper and lower speed given in `v` (in [m s-1]) (if only one `v` is given, the lower boundary
is zero).
"""
function cdf(sd::MBSpeedDistribution{S}, v1::Real, v2::Real) where {S<:AbstractFloat}
    if sd.T*sd.m == 0; return zero(S); end
    a = sd.m / (2 * BOLTZMANN_CONSTANT * sd.T)
    c1 = (sqrt(pi) * erf(sqrt(a) * v1)) / (4 * a^(3/2)) - (v1 * exp(-a*v1^2)) / (2*a)
    c2 = (sqrt(pi) * erf(sqrt(a) * v2)) / (4 * a^(3/2)) - (v2 * exp(-a*v2^2)) / (2*a)
    return S(4 * pi * (a/pi)^(3/2) * (c2 - c1))
end
function cdf(sd::MBFluxSpeedDistribution{S}, v1::Real, v2::Real) where {S<:AbstractFloat}
    if sd.T*sd.m == 0; return zero(S); end
    a = sd.m / (2 * BOLTZMANN_CONSTANT * sd.T)
    c1 = (exp(-a*v1^2) * (a*v1^2 + 1)) / (2 * a^2)
    c2 = (exp(-a*v2^2) * (a*v2^2 + 1)) / (2 * a^2)
    return S(- 2 * a^2 * (c2 - c1))
end
cdf(sd::AbstractMBSpeedDistribution, v::Tuple{Real, Real}) = cdf(sd, v[1], v[2])
cdf(sd::AbstractMBSpeedDistribution, v::Real) = cdf(sd, 0, v)
cdf(S::Type{<:AbstractFloat}, sd::AbstractMBSpeedDistribution, args...) = S(cdf(sd, args...))


"""
    [1] cdf([S::Type{<:AbstractFloat}], vd::AbstractMBVelocityDistribution, 
            v1::Union{Tuple, AbstractVector, LocalCartesianVelocity}, 
            [v2::Union{Tuple, AbstractVector, LocalCartesianVelocity}])

Computes the cumulative distribution function of the velocity vector of the velocity 
distribution `vd`, between the velocity vectors `v1` and `v2` (in [m s-1]). Should only
one vector be given, the lower boundary is assumed to be the zero vector.

Note that if a tuple or vector is given, it is assumed to be given in cartesian coordinates
(see `LocalCartesianVelocity`).
"""
function cdf(vd::MBVelocityDistribution{S}, 
             v1::Union{Tuple, AbstractVector, LocalCartesianVelocity}, 
             v2::Union{Tuple, AbstractVector, LocalCartesianVelocity}) where {S<:AbstractFloat}
    return S(
        cdf(MBSpeedDistribution(vd.T, vd.m), speed(v1), speed(v2)) *
        cdf(MBElevationDistribution(vd.T, vd.m), elevation(v1), elevation(v2)) *
        cdf(MBAzimuthDistribution(vd.T, vd.m), azimuth(v1), azimuth(v2)))
end
function cdf(vd::MBFluxVelocityDistribution{S}, 
             v1::Union{Tuple, AbstractVector, LocalCartesianVelocity}, 
             v2::Union{Tuple, AbstractVector, LocalCartesianVelocity}) where {S<:AbstractFloat}
    return S(
        cdf(MBFluxSpeedDistribution(vd.T, vd.m), speed(v1), speed(v2)) *
        cdf(MBFluxElevationDistribution(vd.T, vd.m), elevation(v1), elevation(v2)) *
        cdf(MBFluxAzimuthDistribution(vd.T, vd.m), azimuth(v1), azimuth(v2)))
end
function cdf(vd::AbstractMBVelocityDistribution, 
             v::Union{Tuple, AbstractVector, LocalCartesianVelocity})
    return cdf(vd, (0,0,0), v)
end
cdf(S::Type{<:AbstractFloat}, vd::AbstractMBVelocityDistribution,  args...) = S(cdf(vd, args...))


"""
    [1] pdf([S::Type{<:AbstractFloat}], ad::AbstractMBAzimuthDistribution, theta::Real)

Calculates the probability density of the azimuth angle of the azimuth distribution `ad`,
at the azimuth angle `theta` (in [rad]). 
"""
function pdf(::AbstractMBAzimuthDistribution{S}, theta::Real) where {S<:AbstractFloat}
    return S(1/(2*pi))
end
pdf(S::Type{<:AbstractFloat}, ad::AbstractMBAzimuthDistribution, args...) = S(pdf(ad, args...))


"""
    [1] pdf([S::Type{<:AbstractFloat}], vd::AbstractMBElevationDistribution, phi::Real)

Calculates the probability density of the elevation angle of the elevation distribution `ed`,
at the elevation angle `phi` (in [rad]).
"""
function pdf(::MBElevationDistribution{S}, phi::Real) where {S<:AbstractFloat}
    return S(0.5 * cos(phi))
end
function pdf(::MBFluxElevationDistribution{S}, phi::Real) where {S<:AbstractFloat}
    return S(2 * sin(phi) * cos(phi))
end
pdf(S::Type{<:AbstractFloat}, ed::AbstractMBElevationDistribution, args...) = S(pdf(ed, args...))


"""
    [1] pdf([S::Type{<:AbstractFloat}], sd::AbstractVelocityDistribution, v::Real)

Calculates the probability density of the speed of the speed distribution `sd`, at the
speed `v` (in [m s-1]). 
"""
function pdf(sd::MBSpeedDistribution{S}, v::Real) where {S<:AbstractFloat}
    if sd.T == 0; return zero(S); end
    a = sd.m / (2 * BOLTZMANN_CONSTANT * sd.T)
    return S((a/pi)^(3//2) * 4*v^2*pi * exp(- a*v^2))
end
function pdf(sd::MBFluxSpeedDistribution{S}, v::Real) where {S<:AbstractFloat}
    if sd.T == 0; return zero(S); end
    a = sd.m / (2 * BOLTZMANN_CONSTANT * sd.T)
    return S(a^2 * v^3 * exp(- a*v^2))
end
pdf(S::Type{<:AbstractFloat}, sd::AbstractMBSpeedDistribution, args...) = S(pdf(sd, args...))


"""
    [1] pdf([S::Type{<:AbstractFloat}], vd::AbstractMBVelocityDistribution, 
            v::Union{Tuple, AbstractVector, LocalCartesianVelocity})

Calculates the probability density of the velocity vector `v` (in [m s-1]) of the velocity
distribution `vd`.
"""
function pdf(vd::MBVelocityDistribution{S}, 
             v::Union{Tuple, AbstractVector, LocalCartesianVelocity}) where {S<:AbstractFloat}
    return S(
        pdf(MBSpeedDistribution(vd.T, vd.m), speed(v)) *
        pdf(MBElevationDistribution(vd.T, vd.m), elevation(v)) *
        pdf(MBAzimuthDistribution(vd.T, vd.m), azimuht(v)))
end
function pdf(vd::MBFluxVelocityDistribution{S}, 
             v::Union{Tuple, AbstractVector, LocalCartesianVelocity}) where {S<:AbstractFloat}
    return S(
        pdf(MBFluxSpeedDistribution(vd.T, vd.m), speed(v)) *
        pdf(MBFluxElevationDistribution(vd.T, vd.m), elevation(v)) *
        pdf(MBFluxAzimuthDistribution(vd.T, vd.m), azimuht(v)))
end
pdf(S::Type{<:AbstractFloat}, vd::AbstractMBVelocityDistribution, args...) = S(pdf(vd, args...))



#::. extensions
# TODO: finish docstrings
"""
# ExESS.jl -- `Base.rand` Extension
"""
Base.rand(::MBAzimuthDistribution{S}) where {S<:AbstractFloat} = S(rand()*2pi - pi)


"""
# ExESS.jl -- `Base.rand` Extension
"""
Base.rand(::MBElevationDistribution{S}) where {S<:AbstractFloat} = S(asin(2*rand() - 1))


"""
# ExESS.jl -- `Base.rand` Extension
"""
function Base.rand(mbsd::MBSpeedDistribution{S}) where {S<:AbstractFloat}
    return S(norm(rand(MBVelocityDistribution(mbsd.T, mbsd.m))))
end


"""
# ExESS.jl -- `Base.rand` Extension
"""
function Base.rand(mbvd::MBVelocityDistribution{S}) where {S<:AbstractFloat}
    a = sqrt(BOLTZMANN_CONSTANT * mbvd.T / mbvd.m)
    return S.((a*randn(), a*randn(), a*randn()))
end


"""
# ExESS.jl -- `Base.rand` Extension
"""
Base.rand(::MBFluxAzimuthDistribution{S}) where {S<:AbstractFloat} = S(rand()*2pi - pi)


"""
# ExESS.jl -- `Base.rand` Extension
"""
Base.rand(::MBFluxElevationDistribution{S}) where {S<:AbstractFloat} = S(acos(sqrt(1 - rand())))


"""
# ExESS.jl -- `Base.rand` Extension
"""
function Base.rand(mbfsd::MBFluxSpeedDistribution{S}) where {S<:AbstractFloat}
    return S(norm(rand(MBFluxVelocityDistribution(mbfsd.T, mbfsd.m))))
end

"""
# ExESS.jl -- `Base.rand` Extension
"""
function Base.rand(mbfvd::MBFluxVelocityDistribution{S}) where {S<:AbstractFloat}
    a = sqrt(BOLTZMANN_CONSTANT * mbfvd.T / mbfvd.m)
    return S.(( a*randn(), a*randn(), a * sqrt(-2*log(1-rand())) ))
end


#::. EXPORTS
export 
    MBAzimuthDistribution, 
    MBElevationDistribution, 
    MBSpeedDistribution, 
    MBVelocityDistribution,

    MBFluxAzimuthDistribution, 
    MBFluxElevationDistribution, 
    MBFluxSpeedDistribution, 
    MBFluxVelocityDistribution,

    cdf, pdf