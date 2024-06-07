############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    [1] solar_angle(theta::Real, phi::Real)
    [2] solar_angle(x::AbstractGlobalPosition)
    [3] solar_angle(grid::AbstractGrid)

Calculate the solar incidence angle based in selenocentric coordinates given
as `theta` and `phi` in radians. Alternatively, the function accepts a positional
argument of a subtype of `<:AbstractGlobalPosition`. Providing a `grid` as the input
argument automatically extracts the angular coordinates and returns all solar angles as a
vector.
"""
function solar_angle(theta::S, phi::S) where {S<:AbstractFloat}
    return S(acos(cos(theta) * cos(phi)))
end
solar_angle(theta::Real, phi::Real) = solar_angle(promote(theta, phi)...)
solar_angle(theta::Integer, phi::Integer) = solar_angle(promote(theta, phi, 1.0)[1:2]...)
solar_angle(xs::GlobalSphericalPosition) = solar_angle(xs.theta, xs.phi)
solar_angle(x::AbstractGlobalPosition) = solar_angle(GlobalSphericalPosition(x))
solar_angle(grid::AbstractGrid) = [solar_angle(c.theta, c.phi) for c in coords(grid)]
solar_angle(S::Type{<:AbstractFloat}, args...) = S(solar_angle(args...))


############################################################################################
#::. EXPORTS
############################################################################################
export solar_angle
