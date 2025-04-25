#::. FUNCTIONS
"""
    escape_velocity(r, m)

Calculates the magnitude of the escape velocity for an orbital radius `r` and body mass `m`.

# Arguments
- `r::Real`: Orbital radius in meters.
- `m::Real`: Mass of the planetary body in kilograms.
"""
escape_velocity(r::T, m::T) where {T<:AbstractFloat} = T(sqrt(2*GRAVITATIONAL_CONSTANT*m/r))
escape_velocity(r::Real, m::Real) = escape_velocity(promote(r, m)...)
escape_velocity(r::Integer, m::Integer) = escape_velocity(promote(r, m, 1.0)[1:2]...)


#::. EXPORTS
export escape_velocity