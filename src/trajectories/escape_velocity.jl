#::. FUNCTIONS
"""
    [1] escape_velocity(r::Real, m::Real)

Calculates the magnitude of the escape velocity for an orbital radius `r` and body mass `m`.
"""
escape_velocity(r::T, m::T) where {T<:AbstractFloat} = T(sqrt(2*GRAVITATIONAL_CONSTANT*m/r))
escape_velocity(r::Real, m::Real) = escape_velocity(promote(r, m)...)
escape_velocity(r::Integer, m::Integer) = escape_velocity(promote(r, m, 1.0)[1:2]...)


#::. EXPORTS
export escape_velocity