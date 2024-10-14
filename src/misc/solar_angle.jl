############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    solar_angle([S], lon, lat)
    solar_angle([S], x)
    solar_angle([S], grid)

Calculate the solar incidence angle based in selenocentric coordinates given
as `lon` and `lat` in radians. The function can also be called with a position `x` or a
grid `grid`. The optional type `S` can be used to convert the output to the desired type.

# Arguments
- (optional) `S::Type{<:AbstractFloat}`: Output type.
- `lon::Real` or `lon::AbstractVector`: Longitude(s) in the range (-π, π).
- `lat::Real` or `lat::AbstractVector`: Latitude(s) in the range (-π/2, π/2).
- `x::GlobalSphericalPosition` or `x::AbstractGlobalPosition` or `AbstractVector`'s of them.
- `grid::AbstractGrid`: Spherical grid.
"""
solar_angle(lon::S, lat::S) where {S<:Real} = acos(cos(lon) * cos(lat))
solar_angle(lon::Real, lat::Real) = solar_angle(promote(lon, lat)...)
solar_angle(lons::AbstractVector, lats::AbstractVector) = solar_angle.(lons, lats)
solar_angle(xs::GlobalSphericalPosition) = solar_angle(xs.theta, xs.phi)
solar_angle(x::AbstractGlobalPosition) = solar_angle(GlobalSphericalPosition(x))
solar_angle(Xs::AbstractVector) = solar_angle.(Xs)
solar_angle(grid::AbstractGrid) = solar_angle(coords(grid))
solar_angle(S::Type{<:AbstractFloat}, args...) = S(solar_angle(args...))


############################################################################################
#::. EXPORTS
############################################################################################
export solar_angle
