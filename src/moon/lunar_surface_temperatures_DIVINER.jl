############################################################################################
#::. FUNCTIONS

# internal union for simplified type handling
_TVGSP = Union{Tuple{Real, Real, Real}, AbstractVector{<:Real}, GlobalSphericalPosition}
############################################################################################
"""
    [1] lunar_surface_temperatures_DIVINER([S::Type{<:AbstractFloat}], theta0::Real)
    [2] lunar_surface_temperatures_DIVINER([S::Type{<:AbstractFloat}], theta0::Real,
                                           theta::Real, phi::Real)
    [3] lunar_surface_temperatures_DIVINER([S::Type{<:AbstractFloat}], theta0::Real,
                                           thetas::AbstractVector, phis::AbstractVector)
    [4] lunar_surface_temperatures_DIVINER([S::Type{<:AbstractFloat}], theta0::Real,
                                           x::Union{Tuple{Real, Real, Real},
                                                    AbstractVector{<:Real},
                                                    GlobalSphericalPosition})
    [5] lunar_surface_temperatures_DIVINER([S::Type{<:AbstractFloat}], theta0::Real,
                                           XS::Vector{Union{Tuple{Real, Real, Real},
                                                            AbstractVector{<:Real},
                                                            GlobalSphericalPosition}})
    [6] lunar_surface_temperatures_DIVINER([S::Type{<:AbstractFloat}], theta0::Real,
                                           grid::AbstractGrid)

Returns the Diviner measurements based lunar surface temperatures, with the subsolar point
shifted by `theta0` (in radians) in the range of [0, 2π] from the center. The input
parameters are in spherical,  sub-solar coordinates, with the longitude `theta` in the range
[-π, π] and the latitude `phi` in the range [-π/2, π/2]. Alternatively, the input can be a
vector of `GlobalSphericalPosition` objects or an `AbstractGrid` object.
"""
function lunar_surface_temperatures_DIVINER(theta0::Real)
    @assert 0 <= theta0 <= 2pi "Longitudinal shift must be in [0, 2pi]!"
    dtheta = rad2deg(theta0)
    S = typeof(dtheta)

    # find correct index & interpolation factor
    idx = findlast(x -> x<=dtheta, 0:15:360)
    l = (idx-1)*15;
    h, f = l + 15, (dtheta - l) / 15

    # setup correct files
    if idx == 1;         high = "diviner_tbol_snapshot_015E.xyz";
                         low  = "diviner_tbol_snapshot_000E.xyz";
    elseif idx == 2;     high = "diviner_tbol_snapshot_030E.xyz";
                         low  = "diviner_tbol_snapshot_015E.xyz";
    elseif 2 < idx < 7;  high = "diviner_tbol_snapshot_0$(h)E.xyz";
                         low  = "diviner_tbol_snapshot_0$(l)E.xyz";
    elseif idx == 7;     high = "diviner_tbol_snapshot_105E.xyz";
                         low  = "diviner_tbol_snapshot_090E.xyz";
    elseif 7 < idx < 24; high = "diviner_tbol_snapshot_$(h)E.xyz";
                         low  = "diviner_tbol_snapshot_$(l)E.xyz";
    else;                high = "diviner_tbol_snapshot_000E.xyz";
                         low  = "diviner_tbol_snapshot_345E.xyz";
    end

    # load lower and higher file
    T_low = readdlm(joinpath(@__DIR__, "..", "..", "data", "diviner_snapshots", low), S)
    T_high = readdlm(joinpath(@__DIR__, "..", "..", "data", "diviner_snapshots", high), S)

    # return temperatures
    return S.(T_low[:,3] .+ f * (T_high[:,3] .- T_low[:,3]))[:]
end
function lunar_surface_temperatures_DIVINER(theta0::S, theta::S, phi::S) where {S<:AbstractFloat}
    @assert -pi <= theta <= pi "Longitude must be in [-pi, pi]!"
    @assert -pi/2 <= phi <= pi/2 "Latitude must be in [-pi/2, pi/2]!"

    T = lunar_surface_temperatures_DIVINER(theta0)
    T = reshape(T, 360, :)' |> Matrix{Float64}

    # setup interpolate object for the temperature map
    T = vcat(T, T[1,:]'); T = hcat(T, T[:, end]); T = hcat(T[:,1], T);
    itp = interpolate(T, BSpline(Linear()))

    # interpolate and return at given SSE coordinates
    idx_theta = (rad2deg(theta) + 179.75) / 359.5 * 719 + 1
    idx_theta += idx_theta < 1 ? 720 : 0
    idx_phi = (rad2deg(phi) +  89.75) / 179.5 * 359 + 2
    return S(itp(idx_theta, idx_phi))
end
function lunar_surface_temperatures_DIVINER(theta0::Real, theta::Real, phi::Real)
    return lunar_surface_temperatures_DIVINER(promote(theta0, theta, phi)...)
end
function lunar_surface_temperatures_DIVINER(theta0::Integer, theta::Integer, phi::Integer)
    return lunar_surface_temperatures_DIVINER(promote(theta0, theta, phi, 1.0)[1:3]...)
end
function lunar_surface_temperatures_DIVINER(theta0::S, thetas::AbstractVector{S},
                                            phis::AbstractVector{S}) where {S<:AbstractFloat}
    T = lunar_surface_temperatures_DIVINER(theta0)
    T = reshape(T, 360, :)' |> Matrix{S}

    # setup interpolate object for the temperature map
    T = vcat(T, T[1,:]'); T = hcat(T, T[:, end]); T = hcat(T[:,1], T);
    itp = interpolate(T, BSpline(Linear()))

    # interpolate at given SSE coordinates
    TT = zeros(S, length(thetas))
    for i in eachindex(thetas)
        @assert -pi <= thetas[i] <= pi "Longitude must be in [-pi, pi]!"
        @assert -pi/2 <= phis[i] <= pi/2 "Latitude must be in [-pi/2, pi/2]!"

        idx_theta = (rad2deg(thetas[i]) + 179.75) / 359.5 * 719 + 1
        idx_theta += idx_theta < 1 ? 720 : 0
        idx_phi = (rad2deg(phis[i]) +  89.75) / 179.5 * 359 + 2
        TT[i] = itp(idx_theta, idx_phi)
    end
    return S.(TT)
end
function lunar_surface_temperatures_DIVINER(theta0::Real, thetas::AbstractVector,
                                            phis::AbstractVector)
    S = typeof(promote(theta0, thetas[1], phis[1])[1])
    return lunar_surface_temperatures_DIVINER(S(theta0), S.(thetas), S.(phis))
end
function lunar_surface_temperatures_DIVINER(theta0::Integer, thetas::AbstractVector{<:Integer},
                                            phis::AbstractVector{<:Integer})
    S = promote_type(typeof(theta0), typeof(thetas[1]), typeof(phis[1]), Float64)
    return lunar_surface_temperatures_DIVINER(S(theta0), S.(thetas), S.(phis))
end
function lunar_surface_temperatures_DIVINER(theta0::Real, x::_TVGSP)
    return lunar_surface_temperatures_DIVINER(theta0, _gettheta(x), _getphi(x))
end
function lunar_surface_temperatures_DIVINER(theta0::Real, X::Vector{_TVGSP})
    thetas, phis = [_gettheta(x) for x in X], [_getphi(x) for x in X]
    return lunar_surface_temperatures_DIVINER(theta0, thetas, phis)
end
function lunar_surface_temperatures_DIVINER(theta0::Real,
                                            X::Vector{GlobalSphericalPosition{S}}) where {S}
    thetas, phis = [_gettheta(x) for x in X], [_getphi(x) for x in X]
    return lunar_surface_temperatures_DIVINER(theta0, thetas, phis)
end
function lunar_surface_temperatures_DIVINER(theta0::Real, grid::AbstractGrid)
    return lunar_surface_temperatures_DIVINER(theta0, surfacecoords(grid))
end
function lunar_surface_temperatures_DIVINER(S::Type{<:AbstractFloat}, args...)
    return S.(lunar_surface_temperatures_DIVINER(args...))
end



"""
    [1] lunar_surface_temperatures_DIVINER_avg([S::Type{<:AbstractFloat}])
    [2] lunar_surface_temperatures_DIVINER_avg([S::Type{<:AbstractFloat}], theta::Real,
                                               phi::Real)
    [3] lunar_surface_temperatures_DIVINER_avg([S::Type{<:AbstractFloat}],
                                               thetas::AbstractVector, phis::AbstractVector)
    [4] lunar_surface_temperatures_DIVINER_avg([S::Type{<:AbstractFloat}],
                                               x::Union{Tuple{Real, Real, Real},
                                                        AbstractVector{<:Real},
                                                        GlobalSphericalPosition})
    [5] lunar_surface_temperatures_DIVINER_avg([S::Type{<:AbstractFloat}],
                                               XS::Vector{Union{Tuple{Real, Real, Real},
                                                                AbstractVector{<:Real},
                                                                GlobalSphericalPosition}})
    [6] lunar_surface_temperatures_DIVINER_avg([S::Type{<:AbstractFloat}],
                                               grid::AbstractGrid)

Returns the Diviner measurements based averaged lunar surface temperatures. The input
parameters are in spherical, sub-solar coordinates, with the longitude `theta` in the range
[-π, π] and the latitude `phi` in the range [-π/2, π/2]. Alternatively, the input can be a
vector of `GlobalSphericalPosition` objects or an `AbstractGrid` object.
"""
function lunar_surface_temperatures_DIVINER_avg()
    T = readdlm(joinpath(@__DIR__,"..","..","data","lunar_surface_temperatures_DIVINER.csv"), ',')
    return T
    return T[:,3]
end
function lunar_surface_temperatures_DIVINER_avg(theta::S, phi::S) where {S<:AbstractFloat}
    @assert -pi <= theta <= pi "Longitude must be in [-pi, pi]!"
    @assert -pi/2 <= phi <= pi/2 "Latitude must be in [-pi/2, pi/2]!"

    T = lunar_surface_temperatures_DIVINER_avg(S)
    T = reshape(T, 360, :)' |> Matrix{S}

    # setup interpophiion object for the temperature map
    T = vcat(T, T[1,:]'); T = hcat(T, T[:, end]); T = hcat(T[:,1], T);
    itp = interpolate(T, BSpline(Linear()))

    # interpolate at given SSE coordinates
    idx_theta = (rad2deg(theta) + 179.75) / 359.5 * 719 + 1
    idx_theta += idx_theta < 1 ? 720 : 0
    idx_phi = (rad2deg(phi) +  89.75) / 179.5 * 359 + 2
    return S(itp(idx_theta, idx_phi))
end
function lunar_surface_temperatures_DIVINER_avg(theta::Real, phi::Real)
    return lunar_surface_temperatures_DIVINER_avg(promote(theta, phi)...)
end
function lunar_surface_temperatures_DIVINER_avg(theta::Integer, phi::Integer)
    return lunar_surface_temperatures_DIVINER_avg(promote(theta, phi, 1.0)[1:2]...)
end
function lunar_surface_temperatures_DIVINER_avg(thetas::AbstractVector{S},
                                                phis::AbstractVector{S}) where {S<:AbstractFloat}
    T = lunar_surface_temperatures_DIVINER_avg(S)
    T = reshape(T, 360, :)' |> Matrix{S}

    # setup interpophiion object for the temperature map
    T = vcat(T, T[1,:]'); T = hcat(T, T[:, end]); T = hcat(T[:,1], T);
    itp = interpolate(T, BSpline(Linear()))

    # interpolate at given SSE coordinates
    TT = zeros(S, length(thetas))
    for i in eachindex(thetas)
        @assert -pi <= thetas[i] <= pi "Longitude must be in [-pi, pi]!"
        @assert -pi/2 <= phis[i] <= pi/2 "Latitude must be in [-pi/2, pi/2]!"

        idx_theta = (rad2deg(thetas[i]) + 179.75) / 359.5 * 719 + 1
        idx_theta += idx_theta < 1 ? 720 : 0
        idx_phi = (rad2deg(phis[i]) +  89.75) / 179.5 * 359 + 2
        TT[i] = itp(idx_theta, idx_phi)
    end
    return TT
end
function lunar_surface_temperatures_DIVINER_avg(thetas::AbstractVector, phis::AbstractVector)
    S = typeof(promote(thetas[1], phis[1])[1])
    return lunar_surface_temperatures_DIVINER_avg(S.(thetas), S.(phis))
end
function lunar_surface_temperatures_DIVINER_avg(thetas::AbstractVector{<:Integer},
                                                phis::AbstractVector{<:Integer})
    S = promote_type(typeof(thetas[1]), typeof(phis[1]), Float64)
    return lunar_surface_temperatures_DIVINER_avg(S.(thetas), S.(phis))
end
function lunar_surface_temperatures_DIVINER_avg(x::_TVGSP)
    return lunar_surface_temperatures_DIVINER_avg(_gettheta(x), _getphi(x))
end
function lunar_surface_temperatures_DIVINER_avg(X::Vector{_TVGSP})
    thetas, phis = [_gettheta(x) for x in X], [_getphi(x) for x in X]
    return lunar_surface_temperatures_DIVINER_avg(thetas, phis)
end
function lunar_surface_temperatures_DIVINER_avg(X::Vector{GlobalSphericalPosition{S}}) where {S}
    thetas, phis = [_gettheta(x) for x in X], [_getphi(x) for x in X]
    return lunar_surface_temperatures_DIVINER_avg(thetas, phis)
end
function lunar_surface_temperatures_DIVINER_avg(grid::AbstractGrid)
    return lunar_surface_temperatures_DIVINER_avg(surfacecoords(grid))
end
function lunar_surface_temperatures_DIVINER_avg(S::Type{<:AbstractFloat}, args...)
    return S.(lunar_surface_temperatures_DIVINER_avg(args...))
end


############################################################################################
#::. EXPORTS
############################################################################################
export
    lunar_surface_temperatures_DIVINER,
    lunar_surface_temperatures_DIVINER_avg
