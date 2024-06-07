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
solar_angle(theta::S, phi::S) where {S<:AbstractFloat} = acos(cos(theta) * cos(phi))
solar_angle(theta::Integer, phi::Integer) = solar_angle(float(theta), float(phi))
solar_angle(args...) = solar_angle(promote(args...)...)
solar_angle(xs::GlobalSphericalPosition) = solar_angle(xs.theta, xs.phi)
solar_angle(x::AbstractGlobalPosition) = solar_angle(GlobalSphericalPosition(x))
solar_angle(grid::AbstractGrid) = [solar_angle(c.theta, c.phi) for c in coords(grid)]
solar_angle(S::Type{<:AbstractFloat}, args...) = S(solar_angle(args...))


############################################################################################
#::. EXPORTS
############################################################################################
export solar_angle
