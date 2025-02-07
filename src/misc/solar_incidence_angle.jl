############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    solar_incidence_angle([S], lon, lat)
    solar_incidence_angle([S], x)
    solar_incidence_angle([S], grid)

Calculate the solar incidence angle based in selenocentric coordinates given as `lon` and
`lat` in radians. The function can also be called with a position `x` or a grid `grid`. The
optional type `S` can be used to convert the output to the desired type.

# Arguments
- `lon::Real` or `lon::AbstractVector`: Longitude(s) in the range (-π, π).
    -  OR Vectors of `lon`s, if `lat` is also a vector.
- `lat::Real` or `lat::AbstractVector`: Latitude(s) in the range (-π/2, π/2).
    -  OR Vectors of `lat`s, if `lon` is also a vector.
- `x::GlobalSphericalPosition` or `x::AbstractGlobalPosition`: Position of the point.
    - OR `AbstractVector`'s `x`s.
- `grid::AbstractGrid`: Spherical grid.
- (optional) `S::Type{<:AbstractFloat}`: Output type.
"""
solar_incidence_angle(lon::S, lat::S) where {S<:Real} = acos(cos(lon) * cos(lat))
solar_incidence_angle(lon::Real, lat::Real) = solar_incidence_angle(promote(lon, lat)...)
solar_incidence_angle(lons::AbstractVector, lats::AbstractVector) = solar_incidence_angle.(lons, lats)
solar_incidence_angle(xs::GlobalSphericalPosition) = solar_incidence_angle(xs.theta, xs.phi)
solar_incidence_angle(x::AbstractGlobalPosition) = solar_incidence_angle(GlobalSphericalPosition(x))
solar_incidence_angle(Xs::AbstractVector) = solar_incidence_angle.(Xs)
solar_incidence_angle(grid::AbstractGrid) = solar_incidence_angle(coords(grid))
solar_incidence_angle(S::Type{<:AbstractFloat}, args...) = S(solar_incidence_angle(args...))


"""
    local_solar_incidence_angle([S], lon, lat, slope, az, [decl])

Calculate the local solar incidence angle based in selenocentric coordinates given as `lon`
and `lat` in radians, the local slope `slope` in radians, and the azimuth `az` in radians.
The solar declination `decl` can also be provided to include seasonal effects. The optional
type `S` can be used to convert the output to the desired type.

# Arguments
- `lon::Real`: Longitude in the range (-π, π).
- `lat::Real`: Latitude in the range (-π/2, π/2).
- `slope::Real`: Local slope in radians.
- `az::Real`: Local slope azimuth in radians (0 towards south - negative towards east).
- `decl::Real`: Solar declination in radians.
- (optional) `S::Type{<:AbstractFloat}`: Output type.

Note that all inputs (but the type `S`) can be used as vectors.

# References
- Duffie, J. A., & Beckman, W. A. (2013). Solar Engineering of Thermal Processes
  (No. 4; Fourth Ed.). Hoboken: Wiley.
- Gallinger et al. (2024). Thermophysical Diversity of Young Lunar Crater Ejecta Revealed
  with LRO Diviner Observations. The Planetary Science Journal, 5(11), 261.
  DOI: 10.3847/psj/ad84e3
"""
function local_solar_incidence_angle(lon::Real, lat::Real, slope::Real, az::Real)
    S = eltype(promote(lon, lat, slope, az)) # for type-stability
    return S(acos(cos(lat)   * cos(slope)            * cos(lon) +
                  sin(lat)   * sin(slope) * cos(az)  * cos(lon) +
                               sin(slope) * sin(az)  * sin(lon) ))
end
function local_solar_incidence_angle(lon::Real, lat::Real, slope::Real, az::Real, decl::Real)
    S = eltype(promote(lon, lat, slope, az, decl)) # for type-stability
    return S(acos(sin(decl) * (sin(lat) * cos(slope) - cos(lat) * sin(slope) * cos(az)) +
                  cos(decl) * cos(local_solar_incidence_angle(lon, lat, slope, az))))
end
function local_solar_incidence_angle(x::GlobalSphericalPosition, args...)
    return local_solar_incidence_angle(x.theta, x.phi, args...)
end
function local_solar_incidence_angle(x::AbstractGlobalPosition, args...)
    return local_solar_incidence_angle(GlobalSphericalPosition(x), args...)
end
function local_solar_incidence_angle(S::Type{<:AbstractFloat}, args...)
    return S(local_solar_incidence_angle(args...))
end
function local_solar_incidence_angle(args...) # for vectorized inputs
    return local_solar_incidence_angle.(args...)
end


############################################################################################
#::. EXPORTS
############################################################################################
export solar_incidence_angle, local_solar_incidence_angle
