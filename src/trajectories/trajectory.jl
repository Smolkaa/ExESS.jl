############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    trajectory(x0, v0; kwargs...)

Calculate the trajectory of a particle starting at position `x0` (global cartesian
coordinates), with initial velocity `v0` (global cartesian coordinates), given the
acceleration function `ddx` (global cartesian coordinates). Returns a `ODESolution` object 
as the trajectory.

# Arguments
- `x0::Tuple` or `x0::AbstractVector` or `x0::AbstractPosition`: Initial position of the 
  particle in global cartesian coordinates.
- `v0::Tuple` or `v0::AbstractVector` or `v0::AbstractVelocity`: Initial velocity of the
  particle in global cartesian coordinates.

# Key-Word Arguments
- `alg=Tsit5()`: Numerical solver algorithm.
- `ddx::Function=ddx_gravity`: Acceleration function of the particle in glob. cart. coords.
- `rmin::Real=LUNAR_RADIUS`: Minimum radius of computational domain (m).
- `rmax::Real=1e9`: Maximum radius of computational domain (m).
- `tspan::Tuple=(0f0,1f10)`: Time span of the integration ((s), (s)).

# Notes
- 16-bit types are not accurate enough and will be promoted (no type-stability).
- The integration terminates if either the minimum or maximum radius is exceeded or if the
end of the time span is reached.
"""
function trajectory(x0::NTuple{3, T}, v0::NTuple{3, T}; ddx::Function=ddx_gravity,
        alg=Tsit5(), rmin::Real=LUNAR_RADIUS, rmax::Real=1e9, tspan::Tuple=(0f0,1f10),
        kwargs...) where {T<:AbstractFloat}

    # input checks and conversions
    @assert norm(x0) > rmin "Starting position `x0` must be greater than `rmin`!"
    x0, v0 = SA[collect(x0)...], SA[collect(v0)...]

    # set up of the 2nd-order ODE problem
    function _f(du, u, p, t); return SA[ddx(u, du)...]; end
    prob = SecondOrderODEProblem(_f, v0, x0, tspan)

    # set up of event handling callbacks
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
    return trajectory(Tuple(x0), Tuple(v0), args...; kwargs...)
end
function trajectory(x0::GlobalCartesianPosition, v0::GlobalCartesianVelocity, args...; kwargs...)
    return trajectory(Tuple(x0), Tuple(v0), args...; kwargs...)
end
function trajectory(x0::AbstractPosition, v0::AbstractVelocity, args...; kwargs...)
    return trajectory(GlobalCartesianPosition(x0), GlobalCartesianVelocity(x0, v0), args...; kwargs...)
end


############################################################################################
#::. UTILITY FUNCTIONS
############################################################################################
"""
    ddx_gravity(x; kwargs...)

Acceleration function for gravity. Assumes a global cartesian coordinate system, which is
centered at the center of the central object.

# Arguments
- `x::Tuple` or `x::AbstractVector` or `x::AbstractPosition`: Position of the particle in
  global cartesian coordinates (also accepts `GlobalSphericalPosition` type).

# Key-Word Arguments
- `m::Real=LUNAR_MASS`: mass of central object (kg).
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
function ddx_gravity(x::AbstractVector; kwargs...)
    return ddx_gravity(Tuple(x); kwargs...)
end
function ddx_gravity(x::GlobalCartesianPosition; kwargs...)
    return ddx_gravity(Tuple(x); kwargs...)
end
function ddx_gravity(x::GlobalSphericalPosition; kwargs...)
    return ddx_gravity(GlobalCartesianPosition(x); kwargs...)
end
ddx_gravity(x, args...; kwargs...) = ddx_gravity(x; kwargs...)


"""
    ddx0(args...; kwargs...) = (0, 0, 0)

Accelerator function for zero acceleration. Used to simulate ray tracing. 
"""
ddx0(x::NTuple{3, T}) where {T<:AbstractFloat} = T.((0,0,0))
ddx0(x::Tuple{<:Real, <:Real, <:Real}) = ddx0((promote(x...)); kwargs...)
ddx0(x::Tuple{<:Integer, <:Integer, <:Integer}) = ddx0(float.(x); kwargs...)
ddx0(x::AbstractVector) = ddx0(Tuple(x))
ddx0(x::GlobalCartesianPosition) = ddx0(Tuple(x))
ddx0(x::GlobalSphericalPosition) = ddx0(GlobalCartesianPosition(x))
ddx0(x, args...; kwargs...) = ddx0(x)


############################################################################################
#::. EXTENSIONS
############################################################################################
"""
    GlobalCartesianPosition(sol, t)

Extracts `GlobalCartesianPosition` from the trajectory `sol` at time or times `t`. 
"""
function GlobalCartesianPosition(sol::ODESolution, t::Real)
    solt = sol(float(t))
    return GlobalCartesianPosition(solt[4], solt[5], solt[6])
end
function GlobalCartesianPosition(sol::ODESolution, t::AbstractVector)
    return [GlobalCartesianPosition(sol, tt) for tt in t]
end


"""
    GlobalSphericalPosition(sol, t)

Extracts `GlobalSphericalPosition` from the trajectory `sol` at time or times `t`. 
"""
function GlobalSphericalPosition(sol::ODESolution, t::Real)
    return GlobalSphericalPosition(GlobalCartesianPosition(sol, t))
end
function GlobalSphericalPosition(sol::ODESolution, t::AbstractVector)
    return [GlobalSphericalPosition(sol, tt) for tt in t]
end


"""
    GlobalCartesianVelocity(sol, t)

Extracts `GlobalCartesianVelocity` from the trajectory `sol` at time or times `t`. 
"""
function GlobalCartesianVelocity(sol::ODESolution, t::Real)
    solt = sol(float(t))
    return GlobalCartesianVelocity(solt[1], solt[2], solt[3])
end
function GlobalCartesianVelocity(sol::ODESolution, t::AbstractVector)
    return [GlobalCartesianVelocity(sol, tt) for tt in t]
end


"""
    LocalCartesianVelocity(sol, t)

Extracts `LocalCartesianVelocity` from the trajectory `sol` at time or times `t`. 
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
