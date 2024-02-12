#::. functions
"""
    [1] lunar_surface_temperatures_BUTLER1997(theta::Real, phi::Real)
    [2] lunar_surface_temperatures_BUTLER1997(thetas::AbstractVector, phis::AbstractVector; matrix=true)
    [3] lunar_surface_temperatures_BUTLER1997(xs::GlobalSphericalPosition)
    [4] lunar_surface_temperatures_BUTLER1997(XS::Vector{GlobalSphericalPosition})
    [5] lunar_surface_temperatures_BUTLER1997(grid::AbstractGrid)

Calculates the surface temperature based on the approximation given in Butler, 1997. Takes
the SSE coordinates `theta`/`thetas` (subsolar longitudes) and `phi`/`phis` (subsolar latitudes)
and returns the temperature vector of equal size at the given pair. Alternatively, a 
`GlobalSphericalPosition`/`Vector{GlobalSphericalPosition}` or a `grid` can be passed,
from which the subsolar coordinates will be extracted.

Note that the range for the sunsolar longitude is `[-pi, pi]`, and for the subsolar latitude
is `[-pi/2, pi/2]`. 

`T = 250*cos(Z)^1/4 + 100` on the Sun-side, and `T = 100` on the night-side, in kelvin.
"""
function lunar_surface_temperatures_BUTLER1997(theta::S, phi::S) where {S<:AbstractFloat}
    @assert -pi <= theta <= pi "Longitude must be in [-pi, pi]!"
    @assert -pi/2 <= phi <= pi/2 "Latitude must be in [-pi/2, pi/2]!"

    return abs(theta) < pi/2 ? S(250*(cos(theta) * cos(phi))^(1/4) + 100) : S(100)
end
function lunar_surface_temperatures_BUTLER1997(theta::Real, phi::Real) 
    return lunar_surface_temperatures_BUTLER1997(promote(theta, phi)...)
end
function lunar_surface_temperatures_BUTLER1997(theta::Integer, phi::Integer) 
    return lunar_surface_temperatures_BUTLER1997(promote(theta, phi, 1.0)[1:2]...)
end
function lunar_surface_temperatures_BUTLER1997(thetas::AbstractVector, phis::AbstractVector; matrix=true)
    return matrix ? 
        [lunar_surface_temperatures_BUTLER1997(theta, phi) for theta in thetas, phi in phis] :
        [lunar_surface_temperatures_BUTLER1997(thetas[i], phis[i]) for i in eachindex(thetas)]
end
function lunar_surface_temperatures_BUTLER1997(xs::GlobalSphericalPosition)
    return lunar_surface_temperatures_BUTLER1997(xs.theta, xs.phi)
end
function lunar_surface_temperatures_BUTLER1997(XS::Vector{GlobalSphericalPosition{S}}) where {S} 
    return lunar_surface_temperatures_BUTLER1997.(XS)
end
function lunar_surface_temperatures_BUTLER1997(grid::AbstractGrid)
    return lunar_surface_temperatures_BUTLER1997(surfacecoords(grid))
end


"""
    [1] lunar_surface_temperatures_diviner(theta0::Real)
    [2] lunar_surface_temperatures_diviner(theta0::Real, theta::Real, phi::Real)
    [3] lunar_surface_temperatures_diviner(theta0::Real, thetas::AbstractVector, phis::AbstractVector; matrix=true)
    [4] lunar_surface_temperatures_diviner(theta0::Real, xs::GlobalSphericalPosition) 
    [5] lunar_surface_temperatures_diviner(theta0::Real, XS::Vector{GlobalSphericalPosition}) 
    [6] lunar_surface_temperatures_diviner(theta0::Real, grid::Abstract2DGrid)

Returns the Diviner measurements based lunar surface temperatures, with the subsolar point
shifted by `theta` (in radians) from the center.
"""
function lunar_surface_temperatures_diviner(theta0::Real)
    @assert 0 <= theta0 <= 2pi "Longitude shift must be in [0, 2pi]!"
    dtheta = rad2deg(theta0)
    S = typeof(dtheta)

    # find correct index
    idx = findlast(x -> x<=dtheta, 0:15:360)
    l = (idx-1)*15; 
    h, f = l + 15, (dtheta - l) / 15 

    # setup correct files
    if idx == 1;         high =  "diviner_tbol_snapshot_015E.xyz";   low = "diviner_tbol_snapshot_000E.xyz";
    elseif idx == 2;     high =  "diviner_tbol_snapshot_030E.xyz";   low = "diviner_tbol_snapshot_015E.xyz";
    elseif 2 < idx < 7;  high =  "diviner_tbol_snapshot_0$(h)E.xyz"; low = "diviner_tbol_snapshot_0$(l)E.xyz";
    elseif idx == 7;     high =  "diviner_tbol_snapshot_105E.xyz";   low = "diviner_tbol_snapshot_090E.xyz";
    elseif 7 < idx < 24; high =  "diviner_tbol_snapshot_$(h)E.xyz";  low = "diviner_tbol_snapshot_$(l)E.xyz";
    else;                high =  "diviner_tbol_snapshot_000E.xyz";   low = "diviner_tbol_snapshot_345E.xyz";
    end

    # load lower and higher file
    T_low = readdlm(joinpath(@__DIR__, "..", "..", "data", "diviner_snapshots", low), S)
    T_high = readdlm(joinpath(@__DIR__, "..", "..", "data", "diviner_snapshots", high), S)
    
    return deg2rad.(T_low[:,1]), deg2rad.(T_low[:,2]), T_low[:,3] .+ f * (T_high[:,3] .- T_low[:,3])
end
function lunar_surface_temperatures_diviner(theta0::S, theta::S, phi::S) where {S<:AbstractFloat}
    @assert -pi <= theta <= pi "Longitude must be in [-pi, pi]!"
    @assert -pi/2 <= phi <= pi/2 "Latitude must be in [-pi/2, pi/2]!"

    _, _, T = lunar_surface_temperatures_diviner(theta0)
    T = reshape(T, 360, :)' |> Matrix{Float64}
    
    # setup interpolate object for the temperature map
    T = vcat(T, T[1,:]'); T = hcat(T, T[:, end]); T = hcat(T[:,1], T);
    itp = interpolate(T, BSpline(Linear()))

    # interpolate at given SSE coordinates
    idx_theta = (rad2deg(theta) + 179.75) / 359.5 * 719 + 1
    idx_theta += idx_theta < 1 ? 720 : 0
    idx_phi = (rad2deg(phi) +  89.75) / 179.5 * 359 + 2
    return S(itp(idx_theta, idx_phi))
end 
lunar_surface_temperatures_diviner(theta0::Real, theta::Real, phi::Real) = lunar_surface_temperatures_diviner(promote(theta0, theta, phi)...)
lunar_surface_temperatures_diviner(theta0::Integer, theta::Integer, phi::Integer) = lunar_surface_temperatures_diviner(promote(theta0, theta, phi, 1.0)[1:3]...)
function lunar_surface_temperatures_diviner(theta0::S, thetas::AbstractVector{S}, phis::AbstractVector{S}; matrix=true) where {S<:AbstractFloat}
    _, _, T = lunar_surface_temperatures_diviner(theta0)
    T = reshape(T, 360, :)' |> Matrix{Float64}

    # setup interpolate object for the temperature map
    T = vcat(T, T[1,:]'); T = hcat(T, T[:, end]); T = hcat(T[:,1], T);
    itp = interpolate(T, BSpline(Linear()))

    # interpolate at given SSE coordinates
    TT = matrix ? zeros(S, length(thetas), length(phis)) : zeros(S, length(thetas))
    for i in eachindex(thetas), j in eachindex(phis)*matrix
        @assert -pi <= thetas[i] <= pi "Longitude must be in [-pi, pi]!"
        @assert -pi/2 <= phis[matrix ? j : i] <= pi/2 "Latitude must be in [-pi/2, pi/2]!"

        idx_theta = (rad2deg(thetas[i]) + 179.75) / 359.5 * 719 + 1
        idx_theta += idx_theta < 1 ? 720 : 0
        idx_phi = (rad2deg(phis[matrix ? j : i]) +  89.75) / 179.5 * 359 + 2
        TT[i,matrix ? j : 1] = itp(idx_theta, idx_phi)
    end
    return TT
end 
function lunar_surface_temperatures_diviner(theta0::Real, thetas::AbstractVector, phis::AbstractVector; kwargs...) 
    S = typeof(promote(theta0, thetas..., phis...)[1])
    if S <: Integer; S = typeof(promote(one(S), 1.0)[1]); end
    return lunar_surface_temperatures_diviner(S(theta0), S.(thetas), S.(phis); kwargs...)
end
function lunar_surface_temperatures_diviner(theta0::Real, xs::GlobalSphericalPosition) 
    return lunar_surface_temperatures_diviner(theta0, xs.theta, xs.phi)
end
function lunar_surface_temperatures_diviner(theta0::Real, XS::Vector{GlobalSphericalPosition{S}}) where {S}
    thetas, phis = [xs.theta for xs in XS], [xs.phi for xs in XS]
    return lunar_surface_temperatures_diviner(theta0, thetas, phis; matrix=false)
end
lunar_surface_temperatures_diviner(theta0::Real, grid::AbstractGrid) = lunar_surface_temperatures_diviner(theta0, surfacecoords(grid))


"""
    [1] lunar_surface_temperatures_diviner_avg()
    [2] lunar_surface_temperatures_diviner_avg(theta::Real, phi::Real)
    [3] lunar_surface_temperatures_diviner_avg(thetas::AbstractVector, phis::AbstractVector; matrix=true)
    [4] lunar_surface_temperatures_diviner_avg(xs::GlobalSphericalPosition)
    [5] lunar_surface_temperatures_diviner_avg(XS::Vector{GlobalSphericalPosition})
    [6] lunar_surface_temperatures_diviner_avg(grid::Abstract2DGrid)

Returns the Diviner measurements based averaged lunar surface temperatures.
"""
function lunar_surface_temperatures_diviner_avg()
    T = readdlm(joinpath(@__DIR__, "..", "..", "data", "lunar_surface_temperatures_diviner.csv"), '\t')
    return deg2rad.(T[:,1]), deg2rad.(T[:,2]), T[:,3]
end
function lunar_surface_temperatures_diviner_avg(theta::S, phi::S) where {S<:AbstractFloat}
    @assert -pi <= theta <= pi "Longitude must be in [-pi, pi]!"
    @assert -pi/2 <= phi <= pi/2 "Latitude must be in [-pi/2, pi/2]!"

    _, _, T = lunar_surface_temperatures_diviner_avg()
    T = reshape(T, 360, :)' |> Matrix{Float64}

    # setup interpophiion object for the temperature map
    T = vcat(T, T[1,:]'); T = hcat(T, T[:, end]); T = hcat(T[:,1], T);
    itp = interpolate(T, BSpline(Linear()))

    # interpolate at given SSE coordinates
    idx_theta = (rad2deg(theta) + 179.75) / 359.5 * 719 + 1
    idx_theta += idx_theta < 1 ? 720 : 0
    idx_phi = (rad2deg(phi) +  89.75) / 179.5 * 359 + 2
    return S(itp(idx_theta, idx_phi))
end
lunar_surface_temperatures_diviner_avg(theta::Real, phi::Real) = lunar_surface_temperatures_diviner_avg(promote(theta, phi)...)
lunar_surface_temperatures_diviner_avg(theta::Integer, phi::Integer) = lunar_surface_temperatures_diviner_avg(promote(theta, phi, 1.0)[1:2]...)
function lunar_surface_temperatures_diviner_avg(thetas::AbstractVector{S}, phis::AbstractVector{S}; matrix=true) where {S<:AbstractFloat}
    _, _, T = lunar_surface_temperatures_diviner_avg()
    T = reshape(T, 360, :)' |> Matrix{Float64}

    # setup interpophiion object for the temperature map
    T = vcat(T, T[1,:]'); T = hcat(T, T[:, end]); T = hcat(T[:,1], T);
    itp = interpolate(T, BSpline(Linear()))

    # interpolate at given SSE coordinates
    TT = matrix ? zeros(S, length(thetas), length(phis)) : zeros(S, length(thetas))
    for i in eachindex(thetas), j in eachindex(phis)*matrix
        @assert -pi <= thetas[i] <= pi "Longitude must be in [-pi, pi]!"
        @assert -pi/2 <= phis[matrix ? j : i] <= pi/2 "Latitude must be in [-pi/2, pi/2]!"

        idx_theta = (rad2deg(thetas[i]) + 179.75) / 359.5 * 719 + 1
        idx_theta += idx_theta < 1 ? 720 : 0
        idx_phi = (rad2deg(phis[matrix ? j : i]) +  89.75) / 179.5 * 359 + 2
        TT[i,matrix ? j : 1] = itp(idx_theta, idx_phi)
    end
    return TT
end
function lunar_surface_temperatures_diviner_avg(thetas::AbstractVector, phis::AbstractVector; kwargs...)
    S = typeof(promote(thetas..., phis...)[1])
    if S <: Integer; S = typeof(promote(one(S), 1.0)[1]); end
    return lunar_surface_temperatures_diviner_avg(S.(thetas), S.(phis); kwargs...)
end
function lunar_surface_temperatures_diviner_avg(xs::GlobalSphericalPosition) 
    return lunar_surface_temperatures_diviner_avg(xs.theta, xs.phi)
end
function lunar_surface_temperatures_diviner_avg(XS::Vector{GlobalSphericalPosition{S}}) where {S}
    thetas, phis = [xs.theta for xs in XS], [xs.phi for xs in XS]
    return lunar_surface_temperatures_diviner_avg(thetas, phis; matrix=false)
end
lunar_surface_temperatures_diviner_avg(grid::AbstractGrid) = lunar_surface_temperatures_diviner_avg(surfacecoords(grid))


"""
    [1] lunar_surface_temperatures_HURLEY2015(theta::Real, phi::Real)
    [2] lunar_surface_temperatures_HURLEY2015(thetas::AbstractVector, phis::AbstractVector; matrix=true)
    [3] lunar_surface_temperatures_HURLEY2015(xs::GlobalSphericalPosition)
    [4] lunar_surface_temperatures_HURLEY2015(XS::Vector{GlobalSphericalPosition})
    [5] lunar_surface_temperatures_HURLEY2015(grid::Abstract2DGrid)

Calcuphies the lunar surface temperatures based on the analytic formula given in Hurley et
al. 2015. All angular arguments are in radians.
"""
function lunar_surface_temperatures_HURLEY2015(theta::S, phi::S) where {S<:AbstractFloat}
    @assert -pi <= theta <= pi "Longitude must be in [-pi, pi]!"
    @assert -pi/2 <= phi <= pi/2 "Latitude must be in [-pi/2, pi/2]!"   

    if abs(theta) >= pi/2
        a = [444.738, -448.937, 239.668, -63.8844, 8.34064, -0.423502]
        if theta < 0; theta += 2pi; end
        cophi = -(phi - pi/2)
        return S(sum([a[i] * theta^(i-1) for i in 1:6]) + 35 * (sin(cophi)-1))
    end
    return S(262*(cos(theta) * cos(phi))^(1/2) + 130)
end
lunar_surface_temperatures_HURLEY2015(theta::Real, phi::Real) = lunar_surface_temperatures_HURLEY2015(promote(theta, phi)...)
lunar_surface_temperatures_HURLEY2015(theta::Integer, phi::Integer) = lunar_surface_temperatures_HURLEY2015(promote(theta, phi, 1.0)[1:2]...)
function lunar_surface_temperatures_HURLEY2015(thetas::AbstractVector, phis::AbstractVector; matrix=true)
    return matrix ? 
        [lunar_surface_temperatures_HURLEY2015(theta, phi) for theta in thetas, phi in phis] :
        [lunar_surface_temperatures_HURLEY2015(thetas[i], phis[i]) for i in eachindex(thetas)]
end
function lunar_surface_temperatures_HURLEY2015(xs::GlobalSphericalPosition)
    return lunar_surface_temperatures_HURLEY2015(xs.theta, xs.phi)
end
function lunar_surface_temperatures_HURLEY2015(XS::Vector{GlobalSphericalPosition{S}}) where {S}
    return lunar_surface_temperatures_HURLEY2015.(XS)
end
function lunar_surface_temperatures_HURLEY2015(grid::AbstractGrid)
    return lunar_surface_temperatures_HURLEY2015(surfacecoords(grid))
end


#::. exports
export 
    lunar_surface_temperatures_diviner, 
    lunar_surface_temperatures_diviner_avg, 
    lunar_surface_temperatures_BUTLER1997, 
    lunar_surface_temperatures_HURLEY2015