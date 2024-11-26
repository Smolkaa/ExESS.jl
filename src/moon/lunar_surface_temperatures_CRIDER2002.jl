############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    lunar_surface_temperatures_CRIDER2002([S], x; kwargs...)
    lunar_surface_temperatures_CRIDER2002([S], lon, lat; kwargs...)
    lunar_surface_temperatures_CRIDER2002([S], grid; kwargs...)

Calculates the lunar surface temperatures based on the analytic formula given in Crider et
al. 2002. The input is given either as the solar zenith angle `x` in radians, or a sub-solar
position `lon` and `lat` in radians. Both approaches can be used with vectorized inputs as
well as a grid of points defined in `grid`. Note that only the surface grids points are
used for the surface temperature calculation (see `surfacecoords`). Additionally, the user
can provide a type `S` to convert the output to the desired type.

# Arguments
- (optional) `S::Type{<:AbstractFloat}`: Output type.
- 'x::Real': Solar zenith angle in the range (-π, π). OR
- `x::GlobalSphericalPosition` or `x::Tuple{Real, Real, Real}` or `x::AbstractVector` or
  an `Abstractvector` with entries of the same: SSE coordinate(s).
- `lon::Real` or `lon::AbstractVector`: Longitude(s) in the range (-π, π).
- `lat::Real` or `lat::AbstractVector`: Latitude(s) in the range (-π/2, π/2).
- `grid::AbstractGrid`: Grid of points to evaluate the temperatures.

# Key-Word Arguments (Based on Crider et al., 2002)
- `T0::Real=100`: Minimum (night time) temperature in Kelvin.
- `T1::Real=280`: Amplitude of cosine temperature variation in Kelvin. In Killen et al.
  2019, this value is assumed to be 250 K.

# Reference(s)
- Butler (1997). The migration of volatiles on the surfaces of Mercury and the Moon.
  DOI: 10.1029/97JE01347
- Crider, D., & Vondrak, R. (2002). Hydrogen migration to the lunar poles by solar wind
  bombardment of the Moon. Advances in Space Research, 30(8), 1869-1874.
- Killen, R. M., Williams, D. R., Park, J., Tucker, O. J., & Kim, S.-J. (2019). The lunar
  neon exosphere seen in LACE data. Icarus, 329, 246-250. DOI: 10.1016/j.icarus.2019.04.018
"""
function lunar_surface_temperatures_CRIDER2002(x::S;T0=100,T1=280) where {S<:AbstractFloat}
    @assert -pi <= x <= pi "Solar zenith angle must be in (-π, π)!"
    return abs(x) < pi/2 ? S(T1*cos(x)^(1/4) + T0) : S(T0)
end
function lunar_surface_temperatures_CRIDER2002(x::Integer; kwargs...)
    return lunar_surface_temperatures_CRIDER2002(float(x); kwargs...)
end
function lunar_surface_temperatures_CRIDER2002(X::AbstractVector; kwargs...)
    return lunar_surface_temperatures_CRIDER2002.(X; kwargs...)inc
end
function lunar_surface_temperatures_CRIDER2002(grid::AbstractGrid; kwargs...)
    return lunar_surface_temperatures_CRIDER2002(surfacecoords(grid); kwargs...)
end
function lunar_surface_temperatures_CRIDER2002(args...; kwargs...)
    return lunar_surface_temperatures_CRIDER2002(solar_angle(args...); kwargs...)
end
function lunar_surface_temperatures_CRIDER2002(S::Type{<:AbstractFloat}, args...; kwargs...)
    return S.(lunar_surface_temperatures_CRIDER2002(args...); kwargs...)
end


############################################################################################
#::. EXPORTS
############################################################################################
export lunar_surface_temperatures_CRIDER2002
