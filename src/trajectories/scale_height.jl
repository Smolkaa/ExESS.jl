############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    scale_height(T, m; kwargs..)

Calculates the scale height of a particle with mass `m`, given the thermal energy
based on the thermodynamic equilibrium temperature `T`, in the exosphere of a planetary
body with radius `R` and mass `M`.

The default values for the planetary body are corresponding to the Moon.

# Arguments
- `T::Real`: thermodynamic equilibrium temperature (K)
- `m::Real`: mass of the particle (kg)

# Key-Word Arguments
- `R::Real=LUNAR_RADIUS`: radius of the planetary body (m)
- `M::Real=LUNAR_MASS`: mass of the planetary body (kg)
"""
function scale_height(T::S, m::S; R::Real=LUNAR_RADIUS, M::Real=LUNAR_MASS) where {S<:AbstractFloat}
    return S(BOLTZMANN_CONSTANT * T * R^2 / (GRAVITATIONAL_CONSTANT * m * M))
end
scale_height(T::Real, m::Real; kwargs...) = scale_height(promote(T, m)...; kwargs...)
scale_height(T::Integer, m::Integer; kwargs...) = scale_height(float(T),float(m); kwargs...)


############################################################################################
#::. EXPORTS
############################################################################################
export scale_height
