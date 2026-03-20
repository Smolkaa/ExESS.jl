############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    mercury_surface_temperatures([S], lon, lat, TAA; kwargs...)

# Arguments
- (optional) `S::Type{<:AbstractFloat}`: Output type.
- `lon::Real` or `lon::AbstractVector`: Longitude(s) in the range (-π, π).
- `lat::Real` or `lat::AbstractVector`: Latitude(s) in the range (-π/2, π/2).
- `TAA::Real` or `TAA::AbstractVector`: True anomaly angle(s) in the range (-π, π).

# Key-Word Arguments
- `T_day_max::Real=700.0`: Maximum dayside temperature at perihelion (K).
- `T_night_max::Real=140.0`: Maximum nightside temperature (K).
- `T_night_min::Real=90.0`: Minimum nightside temperature (K).
- `Power_Temp::Real=0.2`: Power for the cosine dependence of the temperature on the solar zenith angle.

# References (TODO)
- Butler 1997
- Leblanc 2003
- Verkercke?
"""
function mercury_surface_temperatures(lon::S, lat::S, TAA::S; T_day_max=700.0, 
                                      T_night_max=140.0, T_night_min=90.0, 
                                      Power_Temp=0.2) where {S<:AbstractFloat}
                                      
    @assert -pi/2 <= lat <= pi/2 "Latitude must be in (-π/2, π/2)!"
    lon = pclamp(lon, -pi, pi)

    # nightside temperatures
    if lon < -pi/2; return T_night_min + (T_night_min - T_night_max) * (lon + pi/2) / pi
    elseif lon > pi/2; return T_night_max + (T_night_min - T_night_max) * (lon - pi/2) / pi
    end

    # dayside temperatures 
    d_Sun = MERCURY_SEMI_MAJOR_AXIS * (1-MERCURY_ORBITAL_ECCENTRICITY^2) /
            (1+MERCURY_ORBITAL_ECCENTRICITY * cos(TAA)) # distance to Sun
    T_day = T_day_max * sqrt(MERCURY_PERIHELION / d_Sun)
    T_night_avg = 0.5 * (T_night_max + T_night_min)

    # calculate surface temperature
    return S(T_night_avg + (T_day - T_night_avg) * (cos(lon) * cos(lat))^Power_Temp)
end
function mercury_surface_temperatures(lon::Real, lat::Real, TAA::Real; kwargs...) 
    return mercury_surface_temperatures(promote(lon, lat, TAA)...; kwargs...)
end
function mercury_surface_temperatures(lon::Integer, lat::Integer, TAA::Integer; kwargs...) 
    return mercury_surface_temperatures(float(lon), lat, TAA; kwargs...)
end
function mercury_surface_temperatures(lon::AbstractVector, lat::AbstractVector, 
        TAA; kwargs...) 
    return mercury_surface_temperatures.(lon, lat, TAA; kwargs...)
end
function mercury_surface_temperatures(x, TAA; kwargs...) 
    return mercury_surface_temperatures(_gettheta(x), _getphi(x), TAA; kwargs...)
end
function mercury_surface_temperatures(X::AbstractVector, TAA; kwargs...) 
    return mercury_surface_temperatures.(X, TAA; kwargs...)
end
function mercury_surface_temperatures(grid::AbstractGrid, TAA; kwargs...) 
    return mercury_surface_temperatures(surfacecoords(grid), TAA; kwargs...)
end
function mercury_surface_temperatures(S::Type{<:AbstractFloat}, args...)
    return S.(mercury_surface_temperatures(args...))
end


############################################################################################
#::. EXPORTS
############################################################################################
export mercury_surface_temperatures