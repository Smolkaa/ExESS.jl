############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    [1] lunar_surface_temperatures_BUTLER1997([S])
    [2] lunar_surface_temperatures_BUTLER1997([S], lon, lat)
    [3] lunar_surface_temperatures_BUTLER1997([S], x)
    [4] lunar_surface_temperatures_BUTLER1997([S], grid)

Calculates the lunar surface temperatures based on the analytic formula given in Butler
1997. The input parameters are in spherical, sub-solar coordinates, with the longitude
`lon` in the range (-π, π) and the latitude `lat` in the range (-π/2, π/2).  Additionally,
the user can provide a type `S` to convert the output to the desired type.

# Arguments
- (optional) `S::Type{<:AbstractFloat}`: Output type.
- `lon::Real` or `lon::AbstractVector`: Longitude(s) in the range (-π, π).
- `lat::Real` or `lat::AbstractVector`: Latitude(s) in the range (-π/2, π/2).
- `x::GlobalSphericalPosition` or `x::Tuple{Real, Real, Real}` or `x::AbstractVector` or
  an `Abstractvector` with entries of the same: SSE coordinate(s).
- `grid::AbstractGrid`: Grid of points to evaluate the temperatures.

# Reference(s)
- Butler (1997), The migration of volatiles on the surfaces of Mercury and the Moon.
  DOI: 10.1029/97JE01347
"""
function lunar_surface_temperatures_BUTLER1997(lon::S, lat::S) where {S<:AbstractFloat}
    @assert -pi/2 <= lat <= pi/2 "Latitude must be in (-π/2, π/2)!"
    lon = pclamp(lon, -pi, pi)

    return abs(lon) < pi/2 ? S(250*(cos(lon) * cos(lat))^(1/4) + 100) : S(100)
end
function lunar_surface_temperatures_BUTLER1997(lon::Real, lat::Real)
    return lunar_surface_temperatures_BUTLER1997(promote(lon, lat)...)
end
function lunar_surface_temperatures_BUTLER1997(lon::Integer, lat::Integer)
    return lunar_surface_temperatures_BUTLER1997(promote(lon, lat, 1.0)[1:2]...)
end
function lunar_surface_temperatures_BUTLER1997(lons::AbstractVector, lats::AbstractVector)
    return lunar_surface_temperatures_BUTLER1997.(lons, lats)
end
function lunar_surface_temperatures_BUTLER1997(x)
    return lunar_surface_temperatures_BUTLER1997(_gettheta(x), _getphi(x))
end
function lunar_surface_temperatures_BUTLER1997(X::AbstractVector)
    return lunar_surface_temperatures_BUTLER1997.(X)
end
function lunar_surface_temperatures_BUTLER1997(grid::AbstractGrid)
    return lunar_surface_temperatures_BUTLER1997(surfacecoords(grid))
end
function lunar_surface_temperatures_BUTLER1997(S::Type{<:AbstractFloat}, args...)
    return S.(lunar_surface_temperatures_BUTLER1997(args...))
end


############################################################################################
#::. EXPORTS
############################################################################################
export lunar_surface_temperatures_BUTLER1997
