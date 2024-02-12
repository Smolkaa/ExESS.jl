#::. functions
"""
    [1] solar_angle(theta::Real, phi::Real)
    [2] solar_angle(xs::GlobalSphericalPosition)
    [3] solar_angle(grid::AbstractGrid)

Calculate the solar incidence angle based in selenocentric coordinates given
as `theta` and `phi` in radians. Providing a `grid` as the input argument automatically
extracts the angular coordinates and returns all solar angles as a vector.
"""
solar_angle(theta::Real, phi::Real) = acos(cos(theta) * cos(phi))
solar_angle(xs::GlobalSphericalPosition) = solar_angle(xs.theta, xs.phi)
solar_angle(x::AbstractPosition) = solar_angle(GlobalSphericalPosition(x))
solar_angle(grid::AbstractGrid) = [solar_angle(c.theta, c.phi) for c in grid.coords]


#::. exports
export solar_angle

