############################################################################################
#::. FUNCTIONS

# internal union for simplified type handling
_TVGSP = Union{Tuple{Real, Real, Real}, AbstractVector{<:Real}, GlobalSphericalPosition}
############################################################################################
"""
    [1] lunar_surface_temperatures_BUTLER1997([S::Type{<:AbstractFloat}], theta::Real, 
                                              phi::Real)
    [2] lunar_surface_temperatures_BUTLER1997([S::Type{<:AbstractFloat}], 
                                              thetas::AbstractVector, phis::AbstractVector)
    [3] lunar_surface_temperatures_BUTLER1997([S::Type{<:AbstractFloat}], 
                                              x::Union{Tuple{Real, Real, Real}, 
                                                       AbstractVector{<:Real}, 
                                                       GlobalSphericalPosition})
    [4] lunar_surface_temperatures_BUTLER1997([S::Type{<:AbstractFloat}], 
                                               X::Vector{Union{Tuple{Real, Real, Real}, 
                                                               AbstractVector{<:Real}, 
                                                               GlobalSphericalPosition}})
    [5] lunar_surface_temperatures_BUTLER1997([S::Type{<:AbstractFloat}], 
                                              grid::AbstractGrid)

Calculates the lunar surface temperatures based on the analytic formula given in Butler
1997. The input parameters are in spherical, sub-solar coordinates, with the longitude
`theta` in the range [-π, π] and the latitude `phi` in the range [-π/2, π/2]. Alternatively,
the input can be a vector of `GlobalSphericalPosition` objects or an `AbstractGrid` object.

**Reference(s)**

- Butler (1997), "The migration of volatiles on the surfaces of Mercury and the Moon.", 
  DOI: 10.1029/97JE01347
"""
function lunar_surface_temperatures_BUTLER1997(theta::S, phi::S) where {S<:AbstractFloat}
    @assert -pi <= theta <= pi "Longitude must be in [-pi, pi]!"
    @assert -pi/2 <= phi <= pi/2 "Latitude must be in [-pi/2, pi/2]!"

    return abs(theta) < pi/2 ? S(250*(cos(theta) * cos(phi))^(1/4) + 100) : S(100)
end
function lunar_surface_temperatures_BUTLER1997(theta::Real, phi::Real) 
    return lunar_surface_temperatures_BUTLER1997(promote(theta, phi)...)
end
function lunar_surface_temperatures_BUTLER1997(theta::Integer, phi::Integer) 
    return lunar_surface_temperatures_BUTLER1997(promote(theta, phi, 1.0)[1:2]...)
end
function lunar_surface_temperatures_BUTLER1997(thetas::AbstractVector, phis::AbstractVector)
    return lunar_surface_temperatures_BUTLER1997.(thetas, phis)
end
function lunar_surface_temperatures_BUTLER1997(x::_TVGSP)
    return lunar_surface_temperatures_BUTLER1997(_gettheta(x), _getphi(x))
end
function lunar_surface_temperatures_BUTLER1997(X::Vector{_TVGSP})
    return lunar_surface_temperatures_BUTLER1997.(X)
end
function lunar_surface_temperatures_BUTLER1997(X::Vector{GlobalSphericalPosition{S}}) where {S} 
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
export
    lunar_surface_temperatures_BUTLER1997