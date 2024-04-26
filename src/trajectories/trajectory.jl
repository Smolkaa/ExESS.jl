############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    [1] trajectory(x0::AbstractVector, v0::AbstractVector, [ddx::Function]; kwargs...)
    [2] trajectory(x0::AbstractPosition, v0::AbstractVelocity, [ddx::Function]; kwargs...)

Calculate the trajectory of a particle starting at position `x0` (global cartesian
coordinates), with initial velocity `v0` (global cartesian coordinates), given the
acceleration function `ddx` (global cartesian coordinates).

Returns a `ODESolution` object as the trajectory.

**Key-Word Arguments**

| Field    | Value          | Unit       | Description                            |
|:-------- | --------------:|:---------- |:-------------------------------------- |
| `alg`    | `Tsit5()`      |            | numerical solver algorithm             |
| `rmin`   | `LUNAR_RADIUS` | [m]        | minimum radius of computational domain |
| `rmax`   | `1e9`          | [m]        | maximum radius of computational domain |
| `tspan`  | `(0f0,1f10)`   | ([s], [s]) | time span of the integration           |
| `kwargs` |                |            | additional key-word arguments          |

The integration terminates if either the minimum or maximum radius is exceeded or if the
end of the time span is reached.

**Notes**

* 16-bit types are not accurate enough and will be promoted (no type-stability)
* calling the function with `<:Integer` arguments will promote them to `Float64`
"""
function trajectory(x0::NTuple{3, T}, v0::NTuple{3, T}, ddx::Function=ddx_gravity;
        alg=Tsit5(), rmin::Real=LUNAR_RADIUS, rmax::Real=1e9, tspan::Tuple=(0f0,1f10),
        kwargs...) where {T<:AbstractFloat}

    x0 = SA{T}[collect(x0 ./ norm(x0) .* (rmin + 0.1))...] # normalize starting position to surface
    v0 = SA[collect(v0)...]

    function _f(du, u, p, t); return SA[ddx(u, du)...]; end
    prob = SecondOrderODEProblem(_f, v0, x0, tspan)

    function condition(out,u,t,integrator)
        out[1] = norm(u[4:6]) - rmin # terminate if ground is intersected
        out[2] = rmax - norm(u[4:6]) # terminate if particle escapes
    end
    affect!(integrator, idx) = terminate!(integrator)
    cb = VectorContinuousCallback(condition, affect!, 2; save_positions=(false,true))

    return solve(prob, alg; callback=cb, dtmin=1e-2, reltol=1e-4, kwargs...)
end
function trajectory(x0::Tuple{<:Real, <:Real, <:Real},
                    v0::Tuple{<:Real, <:Real, <:Real}, args...; kwargs...)
    x0_r, x0_theta, x0_phi, v0_x, v0_y, v0_z = promote(x0..., v0...)
    return trajectory((x0_r, x0_theta, x0_phi), (v0_x, v0_y, v0_z), args...; kwargs...)
end
function trajectory(x0::Tuple{<:Integer, <:Integer, <:Integer},
                    v0::Tuple{<:Integer, <:Integer, <:Integer}, args...; kwargs...)
    return trajectory(float.(x0), float.(v0), args...; kwargs...)
end
function trajectory(x0::AbstractVector, v0::AbstractVector, args...; kwargs...)
    return trajectory(_get(x0), _get(v0), args...; kwargs...)
end
function trajectory(x0::GlobalCartesianPosition, v0::GlobalCartesianVelocity, args...;
                    kwargs...)
    return trajectory(_get(x0), _get(v0), args...; kwargs...)
end
function trajectory(x0::AbstractPosition, v0::AbstractVelocity, args...; kwargs...)
    return trajectory(GlobalCartesianPosition(x0), GlobalCartesianVelocity(x0, v0), args...;
                      kwargs...)
end
SA

############################################################################################
#::. UTILITY FUNCTIONS
############################################################################################
"""
    [1] ddx_gravity(x::NTuple{3}, [args...]; kwargs...)
    [1] ddx_gravity(x::AbstractVector, [args...]; kwargs...)


Acceleration function for gravity.  Assumes a global cartesian coordinate system, which is
centered at the center of the central object.

**Key-Word Arguments**

| Field    | Default Value  | Unit       | Description                          |
|:-------- | --------------:|:---------- |:------------------------------------ |
| `m`      | `LUNAR_MASS`   | [kg]       | mass of central object               |
"""
function ddx_gravity(x::NTuple{3, T}; m::Real=LUNAR_MASS) where {T<:AbstractFloat}
    return T.( .- x .* GRAVITATIONAL_CONSTANT .* m ./ norm(x)^3)
end
function ddx_gravity(x::Tuple{<:Real, <:Real, <:Real}; kwargs...)
    return ddx_gravity((promote(x...)); kwargs...)
end
function ddx_gravity(x::Tuple{<:Integer, <:Integer, <:Integer})
    return ddx_gravity((promote(1.0, x...)[2:4]); kwargs...)
end
function ddx_gravity(x::AbstractVector{T}; m::Real=LUNAR_MASS) where {T<:AbstractFloat}
    return T.(-x .* GRAVITATIONAL_CONSTANT * m / norm(x)^3)
end
function ddx_gravity(x::AbstractVector; m::Real=LUNAR_MASS)
    return -x .* GRAVITATIONAL_CONSTANT * m / norm(x)^3
end
ddx_gravity(x, args...; kwargs...) = ddx_gravity(x; kwargs...)


############################################################################################
#::. EXTENSIONS
############################################################################################
"""
    [1] GlobalCartesianPosition(sol::ODESolution, t::Real)
    [2] GlobalCartesianPosition(sol::ODESolution, t::AbstractVector)

Extracts `GlobalCartesianPosition` from `ODESolution` at time or times `t`. Note that this
constructor assumes the solution object to be from a `SecondOrderODEProblem`, where the
components of the solution vector are ordered as `[vx, vy, vz, x, y, z]`.
"""
function GlobalCartesianPosition(sol::ODESolution, t::Real)
    solt = sol(float(t))
    return GlobalCartesianPosition(solt[4], solt[5], solt[6])
end
function GlobalCartesianPosition(sol::ODESolution, t::AbstractVector)
    return [GlobalCartesianPosition(sol, tt) for tt in t]
end


"""
    [1] GlobalSphericalPosition(sol::ODESolution, t::Real)
    [2] GlobalSphericalPosition(sol::ODESolution, t::AbstractVector)

Extracts `GlobalSphericalPosition` from `ODESolution` at time or times `t`. Note that this
constructor assumes the solution object to be from a `SecondOrderODEProblem`, where the
components of the solution vector are ordered as `[vx, vy, vz, x, y, z]`.
"""
function GlobalSphericalPosition(sol::ODESolution, t::Real)
    return GlobalSphericalPosition(GlobalCartesianPosition(sol, t))
end
function GlobalSphericalPosition(sol::ODESolution, t::AbstractVector)
    return [GlobalSphericalPosition(sol, tt) for tt in t]
end


"""
    [1] GlobalCartesianVelocity(sol::ODESolution, t::Real)
    [2] GlobalCartesianVelocity(sol::ODESolution, t::AbstractVector)

Extracts `GlobalCartesianVelocity` from `ODESolution` at time or times 't'. Note that this
constructor assumes the solution object to be from a `SecondOrderODEProblem`, where the
components of the solution vector are ordered as `[vx, vy, vz, x, y, z]`.
"""
function GlobalCartesianVelocity(sol::ODESolution, t::Real)
    solt = sol(float(t))
    return GlobalCartesianVelocity(solt[1], solt[2], solt[3])
end
function GlobalCartesianVelocity(sol::ODESolution, t::AbstractVector)
    return [GlobalCartesianVelocity(sol, tt) for tt in t]
end


"""
    [1] LocalCartesianVelocity(sol::ODESolution, t::Real)
    [2] LocalCartesianVelocity(sol::ODESolution, t::AbstractVector)

Extracts `LocalCartesianVelocity` from `ODESolution` at time or times 't'. Note that this
constructor assumes the solution object to be from a `SecondOrderODEProblem`, where the
components of the solution vector are ordered as `[vx, vy, vz, x, y, z]`.
"""
function LocalCartesianVelocity(sol::ODESolution, t::Real)
    return LocalCartesianVelocity(GlobalCartesianPosition(sol, t),
                                  GlobalCartesianVelocity(sol, t))
end
function LocalCartesianVelocity(sol::ODESolution, t::AbstractVector)
    return [LocalCartesianVelocity(sol, tt) for tt in t]
end


############################################################################################
#::. EXPORTS
############################################################################################
export trajectory, ddx_gravity
