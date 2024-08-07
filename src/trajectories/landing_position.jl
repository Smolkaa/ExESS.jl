############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    [1] landing_position(x0::Tuple, v0::Tuple; kwargs...)
    [2] landing_position(x0::AbstractVector, v0::AbstractVector; kwargs...)
    [3] landing_position(x0::AbstractPosition, v0::AbstractVelocity; kwargs...)

Calculates the landing position of a particle starting at position `x0` (global spherical
coordinates: radius, longitude, latitude), with initial velocity `v0` (local cartesian
coordinates: x (east), y (north), z (up)), flying on a ballistic trajectory, i.e., only
influenced by graviational forces.

Returns the landing position as a tuple of global spherical coordinates (radius, longitude,
latitude). Returns `(NaN, NaN, NaN)` if the particle escapes.


**Key-Word Arguments**

| Field    | Default Value  | Unit       | Description                          |
|:-------- | --------------:|:---------- |:------------------------------------ |
| `m`      | `LUNAR_MASS`   | (kg)       | mass of central object               |


**References**

* N. Schörghofer: "USER GUIDE: Planetary-Code-Collection: Thermal and Ice Evolution Models
  for Planetary Surfaces", https://github.com/nschorgh/Planetary-Code-Collection
* B. J. Butler, 1997, “The migration of volatiles on the surfaces of Mercury and the Moon,”
  Journal of Geophysical Research, vol. 102, no. E8, pp. 19,283--19,291,
  doi: 10.1029/97JE01347.
* Kegerreis et al., 2017, "Evidence for a localized source of the argon in the lunar
  exosphere", Journal of Geophysical Research: Planets, American Geophysical Union (AGU),
  122, 2163-2181
"""
function landing_position(x0::NTuple{3, T}, v0::NTuple{3, T};
                          m::Real=LUNAR_MASS) where {T<:AbstractFloat}

    # extract input, check for escape
    r, lon0, lat0 = x0
    v_esc = escape_velocity(r, T(m))
    if norm(v0) > v_esc * 0.99999; return (NaN, NaN, NaN); end

    # orbital mechanics to determine trajectory
    psi    = zenith(v0)                          # zenith angle at launch
    az     = azimuth(v0)                         # azimuth angle at launch
    epskin = _relative_kinetic_energy(v_esc, v0) # relative kinetic energy of trajectory
    e      = _eccentricity(epskin, psi)          # eccentricity of trajectory
    dtheta = 2pi - 2*_true_anomaly(e, epskin)    # true anomaly of trajectory

    # calculate landing position
    sin_lat1 = sin(lat0) * cos(dtheta) + cos(lat0) * sin(dtheta) * sin(az)
    cos_lat1 = sqrt(1 - sin_lat1^2)
    cos_dlon = (cos(dtheta) - sin(lat0) * sin_lat1) / (cos(lat0) * cos_lat1)

    dlon = limit_acos(cos_dlon)  # correct longitude for numerical errors
    if dtheta >= pi; dlon = 2pi - dlon; end # correct for overflowing cos_dlon

    # update spherical elements, return
    lon1 = lon0 + sign(v0[1]) * dlon
    if lon1 > pi; lon1 -= 2pi; elseif lon1 < -pi; lon1 += 2pi; end
    lat1 = asin(sin_lat1)
    return T.((r, lon1, lat1))
end
function landing_position(x0::Tuple{<:Real, <:Real, <:Real},
                          v0::Tuple{<:Real, <:Real, <:Real}; kwargs...)
    x0_r, x0_theta, x0_phi, v0_x, v0_y, v0_z = promote(x0..., v0...)
    return landing_position((x0_r, x0_theta, x0_phi), (v0_x, v0_y, v0_z); kwargs...)
end
function landing_position(x0::Tuple{<:Integer, <:Integer, <:Integer},
                          v0::Tuple{<:Integer, <:Integer, <:Integer}; kwargs...)
    x0_r, x0_theta, x0_phi, v0_x, v0_y, v0_z = promote(x0..., v0..., 1.0)[1:6]
    return landing_position((x0_r, x0_theta, x0_phi), (v0_x, v0_y, v0_z); kwargs...)
end
function landing_position(x0::AbstractVector, v0::AbstractVector; kwargs...)
    return landing_position(_get(x0), _get(v0); kwargs...)
end
function landing_position(x0::GlobalSphericalPosition, v0::LocalCartesianVelocity; kwargs...)
    return GlobalSphericalPosition(landing_position(_get(x0), _get(v0); kwargs...))
end
function landing_position(x0::AbstractPosition, v0::AbstractVelocity; kwargs...)
    return landing_position(GlobalSphericalPosition(x0), LocalCartesianVelocity(x0, v0); kwargs...)
end


"""
    [1] time_of_flight(x0::Tuple, v0::Tuple; kwargs...)
    [2] time_of_flight(x0::AbstractVector, v0::AbstractVector; kwargs...)
    [3] time_of_flight(x0::AbstractPosition, v0::AbstractVelocity; kwargs...)

Calculates the time of flight of a particle on a ballistic trajectory, starting at position
`x0` (global spherical coordinates: radius, longitude, latitude), with initial velocity `v0`
(local cartesian coordinates: x (east), y (north), z (up)).

Returns the time of flight in seconds. Returns `NaN` if the particle escapes.


**Key-Word Arguments**

| Field    | Default Value  | Unit       | Description                          |
|:-------- | --------------:|:---------- |:------------------------------------ |
| `m`      | `LUNAR_MASS`   | [kg]       | mass of central object               |


**References**

* N. Schörghofer: "USER GUIDE: Planetary-Code-Collection: Thermal and Ice Evolution Models
  for Planetary Surfaces", https://github.com/nschorgh/Planetary-Code-Collection
* B. J. Butler, 1997, “The migration of volatiles on the surfaces of Mercury and the Moon,”
  Journal of Geophysical Research, vol. 102, no. E8, pp. 19,283--19,291,
  doi: 10.1029/97JE01347.
* Kegerreis et al., 2017, "Evidence for a localized source of the argon in the lunar
  exosphere", Journal of Geophysical Research: Planets, American Geophysical Union (AGU),
  122, 2163-2181
"""
function time_of_flight(x0::NTuple{3, T}, v0::NTuple{3, T};
                        m::Real=LUNAR_MASS) where {T<:AbstractFloat}
    # extract input, check for escape
    r = _getr(x0)
    v_esc = escape_velocity(r, T(m))
    if norm(v0) > v_esc; return NaN; end

    # orbital mechanics to determine trajectory time
    psi    = zenith(v0)                          # zenith angle at x0
    epskin = _relative_kinetic_energy(v_esc, v0) # relative kinetic energy of trajectory at x0
    e      = _eccentricity(epskin, psi)          # _eccentricity of trajectory
    theta  = _true_anomaly(e, epskin)            # true anomaly at x0 of trajectory
    a      = _semi_major_axis(r, epskin)         # semi-major axis of trajectory
    P      = _orbit_period(a, T(m))              # orbital period of trajectory
    E      = _eccentric_anomaly(e, theta)        # eccentric anomaly at x0 of trajectory
    M      = _mean_anomaly(e, E)                 # mean anomaly at x0 of trajectory

    # time at x0 & time of flight
    t0 = M / 2pi * P
    return T(P - 2*t0)
end
function time_of_flight(x0::Tuple{<:Real, <:Real, <:Real},
                        v0::Tuple{<:Real, <:Real, <:Real}; kwargs...)
    x0_r, x0_theta, x0_phi, v0_x, v0_y, v0_z = promote(x0..., v0...)
    return time_of_flight((x0_r, x0_theta, x0_phi), (v0_x, v0_y, v0_z); kwargs...)
end
function time_of_flight(x0::Tuple{<:Integer, <:Integer, <:Integer},
                        v0::Tuple{<:Integer, <:Integer, <:Integer}; kwargs...)
    x0_r, x0_theta, x0_phi, v0_x, v0_y, v0_z = promote(x0..., v0..., 1.0)[1:6]
    return time_of_flight((x0_r, x0_theta, x0_phi), (v0_x, v0_y, v0_z); kwargs...)
end
function time_of_flight(x0::AbstractVector, v0::AbstractVector; kwargs...)
    return time_of_flight(_get(x0), _get(v0); kwargs...)
end
function time_of_flight(x0::GlobalSphericalPosition, v0::LocalCartesianVelocity; kwargs...)
    return time_of_flight(_get(x0), _get(v0); kwargs...)
end
function time_of_flight(x0::AbstractPosition, v0::AbstractVelocity; kwargs...)
    return time_of_flight(GlobalSphericalPosition(x0), LocalCartesianVelocity(x0, v0); kwargs...)
end



############################################################################################
#::. INTERNALS
############################################################################################
"""
    [1] _eccentric_anomaly(e::Real, theta::Real)

Calculates the eccentric anomaly given the orbit's eccentricity `e` and the true
anomaly `theta`. The result is given in radians [0, 2π).

Internal function, not exported to user.
"""
function _eccentric_anomaly(e::T, theta::T) where {T<:AbstractFloat}
    if e >= 1; T(NaN); end
    E = 2 * atan( sqrt((1-e)/(1+e)) * tan(theta/2) )
    return E < 0 ? E + pi + pi : E
end
_eccentric_anomaly(e::Real, theta::Real) = _eccentric_anomaly(promote(e, theta)...)


"""
    [1] _eccentricity(epskin::Real, psi::Real)

Calculates the eccentricity of an elliptical trajetory based on the squared escape velocity
fraction `epskin`, and the zenith launch angle `psi`.

Internal function, not exported to user.
"""
function _eccentricity(epskin::T, psi::T) where {T<:AbstractFloat}
    return sqrt(1 - 4 * epskin * (1 - epskin) * sin(psi)^2)
end
_eccentricity(epskin::Real, psi::Real) = _eccentricity(promote(epskin, psi)...)


"""
    [1] _mean_anomaly(e::Real, E::Real)

Calculates the mean anomaly, given the orbit's eccentricity `e` and the eccentric
anomaly `E`.

Internal function, not exported to user.
"""
_mean_anomaly(e::T, E::T) where {T<:AbstractFloat} = e < 1.0 ? E - e * sin(E) : T(NaN)
_mean_anomaly(e::Real, E::Real) = _mean_anomaly(promote(e, E)...)
_mean_anomaly(e::Integer, E::Integer) = _mean_anomaly(promote(e, E, 1.0)[1:2]...)


"""
    [1] _orbit_period(a::Real, m::Real)

Calculates the period of an elliptical orbit, based on the semi-major axis `a`, and the
body's mass `m`.

Internal function, not exported to user.
"""
_orbit_period(a::T, m::T) where {T<:AbstractFloat} = T(2 * pi * sqrt(a^3/(m*GRAVITATIONAL_CONSTANT)))
_orbit_period(a::Real, m::Real) = _orbit_period(promote(a, m)...)
_orbit_period(a::Integer, m::Integer) = _orbit_period(promote(a, m, 1.0)[1:2]...)



"""
    [1] _relative_kinetic_energy(vesc::Real, v::Tuple)
    [2] _relative_kinetic_energy(vesc::Real, v::AbstractVector)
    [3] _relative_kinetic_energy(vesc::Real, v::LocalCartesianVelocity)

Calculates the relative kinetic energy of the velocity vector `v` with respect to the escape
velocity `vesc`.

Internal function, not exported to user.
"""
function _relative_kinetic_energy(vesc::Real, v::Union{Tuple, AbstractVector})
    vesc, v1, v2, v3 = promote(vesc, v[1:3]...)
    return (norm((v1, v2, v3)) / vesc)^2
end
function _relative_kinetic_energy(vesc::Real, v::LocalCartesianVelocity)
    return _relative_kinetic_energy(vesc, _get(v))
end


"""
    [1] _semi_major_axis(r::Real, epskin::Real)

Calculates the semi-major axis of an elliptical orbit with the current orbital
radius `r` and the trajectory's squared escape velocity fraction `epskin`.

Internal function, not exported to user.
"""
_semi_major_axis(r::Real, epskin::Real) = r / (2 * (1 - epskin))


"""
    [1] _true_anomaly(e::Real, epskin::Real)

Calculates the true anomaly based on the eccentricity of the orbit `e` and the current
squared escape velocity fraction `epskin`.

Internal function, not exported to user.
"""
_true_anomaly(e::Real, epskin::Real) = limit_acos(1/e * ((1-e^2)/(2-2*epskin) - 1))


############################################################################################
#::. EXPORTS
############################################################################################
export landing_position, time_of_flight
