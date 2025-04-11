############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    lunar_surface_temperatures_HURLEY2015([S], x)
    lunar_surface_temperatures_HURLEY2015([S], lon, lat)
    lunar_surface_temperatures_HURLEY2015([S], grid)

Calculates the lunar surface temperatures based on the analytic formula given in Hurley et
al. 2015. The input parameters are in spherical, sub-solar coordinates, with the longitude
`lon` in the range (-π, π) and the latitude `lat` in the range (-π/2, π/2). Additionally,
the user can provide a type `S` to convert the output to the desired type.

# Arguments
- (optional) `S::Type{<:AbstractFloat}`: Output type.
- `x::GlobalSphericalPosition` or `x::Tuple{Real, Real, Real}` or `x::AbstractVector` or
  an `Abstractvector` with entries of the same: SSE coordinate(s).
- `lon::Real` or `lon::AbstractVector`: Longitude(s) in the range (-π, π).
- `lat::Real` or `lat::AbstractVector`: Latitude(s) in the range (-π/2, π/2).
- `grid::AbstractGrid`: Grid of points to evaluate the temperatures.

# References
- Hurley et al. (2015), An analytic function of lunar surface temperature for exospheric
  modeling, DOI: 10.1016/j.icarus.2014.08.043
"""
function lunar_surface_temperatures_HURLEY2015(lon::S, lat::S) where {S<:AbstractFloat}
    @assert -pi/2 <= lat <= pi/2 "Latitude must be in (-π/2, π/2)!"
    lon = pclamp(lon, -pi, pi)

    if abs(lon) >= pi/2
        a = (444.738, -448.937, 239.668, -63.8844, 8.34064, -0.423502)
        if lon < 0; lon += 2pi; end
        colat = pi/2 - lat
        return S(sum([a[i] * lon^(i-1) for i in 1:6]) + 35 * (sin(colat)-1))
    end
    return S(262*(cos(lon) * cos(lat))^(1/2) + 130)
end
function lunar_surface_temperatures_HURLEY2015(lon::Real, lat::Real)
    return lunar_surface_temperatures_HURLEY2015(promote(lon, lat)...)
end
function lunar_surface_temperatures_HURLEY2015(lon::Integer, lat::Integer)
    return lunar_surface_temperatures_HURLEY2015(float(lon), lat)
end
function lunar_surface_temperatures_HURLEY2015(lons::AbstractVector, lats::AbstractVector)
    return lunar_surface_temperatures_HURLEY2015.(lons, lats)
end
function lunar_surface_temperatures_HURLEY2015(x)
    return lunar_surface_temperatures_HURLEY2015(_gettheta(x), _getphi(x))
end
function lunar_surface_temperatures_HURLEY2015(X::AbstractVector)
    return lunar_surface_temperatures_HURLEY2015.(X)
end
function lunar_surface_temperatures_HURLEY2015(grid::AbstractGrid)
    return lunar_surface_temperatures_HURLEY2015(surfacecoords(grid))
end
function lunar_surface_temperatures_HURLEY2015(S::Type{<:AbstractFloat}, args...)
    return S.(lunar_surface_temperatures_HURLEY2015(args...))
end


############################################################################################
#::. EXPORTS
############################################################################################
export lunar_surface_temperatures_HURLEY2015
