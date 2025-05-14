############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    lunar_surface_temperatures_DIVINER([S], lon0)
    lunar_surface_temperatures_DIVINER([S], lon0, lon, lat)
    lunar_surface_temperatures_DIVINER([S], lon0, x)
    lunar_surface_temperatures_DIVINER([S], lon0, grid)

Returns the Diviner measurements based lunar surface temperatures, with the subsolar point
shifted by `lon0` (in radians) in the range of (0, 2π) from the center. The input
parameters are in spherical,  sub-solar coordinates, with the longitude `lon` in the range
(-π, π) and the latitude `lat` in the range (-π/2, π/2). Additionally, the user can provide
a type `S` to convert the output to the desired type.

# Arguments
- (optional) `S::Type{<:AbstractFloat}`: Output type.
- `lon0::Real`: Subsolar longitude shift.
- `lon::Real` or `lon::AbstractVector`: Longitude(s) in the range (-π, π).
- `lat::Real` or `lat::AbstractVector`: Latitude(s) in the range (-π/2, π/2).
- `x::GlobalSphericalPosition` or `x::Tuple{Real, Real, Real}` or `x::AbstractVector` or
  an `Abstractvector` with entries of the same: SSE coordinate(s).
- `grid::AbstractGrid`: Grid of points to evaluate the temperatures.

# References
- Williams, J.-P., Paige, D. A., Greenhagen, B. T., & Sefton-Nash, E. (2017). The global
  surface temperatures of the Moon as measured by the Diviner Lunar Radiometer Experiment.
  Icarus, 283, 300--325. DOI: 10.1016/j.icarus.2016.08.012
"""
function lunar_surface_temperatures_DIVINER(lon0::Real)
    dlon = rad2deg(pclamp(lon0, 0, 2pi))
    S = typeof(dlon)

    # find correct index & interpolation factor
    idx = findlast(x -> x<=dlon, 0:15:360)
    l = (idx-1)*15;
    h, f = l + 15, (dlon - l) / 15

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
function lunar_surface_temperatures_DIVINER(lon0::S, lon::S, lat::S) where {S<:AbstractFloat}
    @assert -pi/2 <= lat <= pi/2 "Latitude must be in (-π/2, π/2)!"
    lon = pclamp(lon, -pi, pi)

    T = lunar_surface_temperatures_DIVINER(lon0)
    T = reshape(T, 360, :)' |> Matrix{Float64}

    # setup interpolate object for the temperature map
    T = vcat(T, T[1,:]'); T = hcat(T, T[:, end]); T = hcat(T[:,1], T);
    itp = interpolate(T, BSpline(Linear()))

    # interpolate and return at given SSE coordinates
    idx_lon = (rad2deg(lon) + 179.75) / 359.5 * 719 + 1
    idx_lon += idx_lon < 1 ? 720 : 0
    idx_lat = (rad2deg(lat) +  89.75) / 179.5 * 359 + 2
    return S(itp(idx_lon, idx_lat))
end
function lunar_surface_temperatures_DIVINER(lon0::Real, lon::Real, lat::Real)
    return lunar_surface_temperatures_DIVINER(promote(lon0, lon, lat)...)
end
function lunar_surface_temperatures_DIVINER(lon0::Integer, lon::Integer, lat::Integer)
    return lunar_surface_temperatures_DIVINER(promote(lon0, lon, lat, 1.0)[1:3]...)
end
function lunar_surface_temperatures_DIVINER(lon0::S, lons::AbstractVector{S},
                                            lats::AbstractVector{S}) where {S<:AbstractFloat}
    T = lunar_surface_temperatures_DIVINER(lon0)
    T = reshape(T, 360, :)' |> Matrix{S}

    # setup interpolate object for the temperature map
    T = vcat(T, T[1,:]'); T = hcat(T, T[:, end]); T = hcat(T[:,1], T);
    itp = interpolate(T, BSpline(Linear()))

    # interpolate at given SSE coordinates
    TT = zeros(S, length(lons))
    for i in eachindex(lons)
        @assert -pi <= lons[i] <= pi "Longitude must be in [-pi, pi]!"
        @assert -pi/2 <= lats[i] <= pi/2 "Latitude must be in (-π/2, π/2)!"

        idx_lon = (rad2deg(lons[i]) + 179.75) / 359.5 * 719 + 1
        idx_lon += idx_lon < 1 ? 720 : 0
        idx_lat = (rad2deg(lats[i]) +  89.75) / 179.5 * 359 + 2
        TT[i] = itp(idx_lon, idx_lat)
    end
    return S.(TT)
end
function lunar_surface_temperatures_DIVINER(lon0::Real, lons::AbstractVector,
                                            lats::AbstractVector)
    S = typeof(promote(lon0, lons[1], lats[1])[1])
    return lunar_surface_temperatures_DIVINER(S(lon0), S.(lons), S.(lats))
end
function lunar_surface_temperatures_DIVINER(lon0::Integer, lons::AbstractVector{<:Integer},
                                            lats::AbstractVector{<:Integer})
    S = promote_type(typeof(lon0), typeof(lons[1]), typeof(lats[1]), Float64)
    return lunar_surface_temperatures_DIVINER(S(lon0), S.(lons), S.(lats))
end
function lunar_surface_temperatures_DIVINER(lon0::Real, x)
    return lunar_surface_temperatures_DIVINER(lon0, _gettheta(x), _getphi(x))
end
function lunar_surface_temperatures_DIVINER(lon0::Real, X::AbstractVector)
    lons, lats = [_gettheta(x) for x in X], [_getphi(x) for x in X]
    return lunar_surface_temperatures_DIVINER(lon0, lons, lats)
end
function lunar_surface_temperatures_DIVINER(lon0::Real, grid::AbstractGrid)
    return lunar_surface_temperatures_DIVINER(lon0, surfacecoords(grid))
end
function lunar_surface_temperatures_DIVINER(S::Type{<:AbstractFloat}, args...)
    return S.(lunar_surface_temperatures_DIVINER(args...))
end


"""
    lunar_surface_temperatures_DIVINER_avg([S])
    lunar_surface_temperatures_DIVINER_avg([S], lon, lat)
    lunar_surface_temperatures_DIVINER_avg([S], x)
    lunar_surface_temperatures_DIVINER_avg([S], grid)

Returns the Diviner measurements based averaged lunar surface temperatures. The input
parameters are in spherical, sub-solar coordinates, with the longitude `lon` in the range
(-π, π) and the latitude `lat` in the range (-π/2, π/2). Additionally, the user can provide
a type `S` to convert the output to the desired type.

# Arguments
- (optional) `S::Type{<:AbstractFloat}`: Output type.
- `lon::Real` or `lon::AbstractVector`: Longitude(s) in the range (-π, π).
- `lat::Real` or `lat::AbstractVector`: Latitude(s) in the range (-π/2, π/2).
- `x::GlobalSphericalPosition` or `x::Tuple{Real, Real, Real}` or `x::AbstractVector` or
  an `Abstractvector` with entries of the same: SSE coordinate(s).
- `grid::AbstractGrid`: Grid of points to evaluate the temperatures.

# References
- Williams, J.-P., Paige, D. A., Greenhagen, B. T., & Sefton-Nash, E. (2017). The global
  surface temperatures of the Moon as measured by the Diviner Lunar Radiometer Experiment.
  Icarus, 283, 300--325. DOI: 10.1016/j.icarus.2016.08.012
"""
function lunar_surface_temperatures_DIVINER_avg()
    T = readdlm(joinpath(@__DIR__,"..","..","data","lunar_surface_temperatures_DIVINER.dat"),',')
    return T[:,3]
end
function lunar_surface_temperatures_DIVINER_avg(lon::S, lat::S) where {S<:AbstractFloat}
    @assert -pi/2 <= lat <= pi/2 "Latitude must be in (-π/2, π/2)!"
    lon = pclamp(lon, -pi, pi)

    T = lunar_surface_temperatures_DIVINER_avg(S)
    T = reshape(T, 360, :)' |> Matrix{S}

    # setup interpolation object for the temperature map
    T = vcat(T, T[1,:]'); T = hcat(T, T[:, end]); T = hcat(T[:,1], T);
    itp = interpolate(T, BSpline(Linear()))

    # interpolate at given SSE coordinates
    idx_lon = (rad2deg(lon) + 179.75) / 359.5 * 719 + 1
    idx_lon += idx_lon < 1 ? 720 : 0
    idx_lat = (rad2deg(lat) +  89.75) / 179.5 * 359 + 2
    return S(itp(idx_lon, idx_lat))
end
function lunar_surface_temperatures_DIVINER_avg(lon::Real, lat::Real)
    return lunar_surface_temperatures_DIVINER_avg(promote(lon, lat)...)
end
function lunar_surface_temperatures_DIVINER_avg(lon::Integer, lat::Integer)
    return lunar_surface_temperatures_DIVINER_avg(promote(lon, lat, 1.0)[1:2]...)
end
function lunar_surface_temperatures_DIVINER_avg(lons::AbstractVector{S},
                                                lats::AbstractVector{S}) where {S<:AbstractFloat}
    T = lunar_surface_temperatures_DIVINER_avg(S)
    T = reshape(T, 360, :)' |> Matrix{S}

    # setup interpolation object for the temperature map
    T = vcat(T, T[1,:]'); T = hcat(T, T[:, end]); T = hcat(T[:,1], T);
    itp = interpolate(T, BSpline(Linear()))

    # interpolate at given SSE coordinates
    TT = zeros(S, length(lons))
    for i in eachindex(lons)
        @assert -pi/2 <= lats[i] <= pi/2 "Latitude must be in (-pi/2, pi/2)!"
        lons[i] = pclamp(lons[i], -pi, pi)

        idx_lon = (rad2deg(lons[i]) + 179.75) / 359.5 * 719 + 1
        idx_lon += idx_lon < 1 ? 720 : 0
        idx_lat = (rad2deg(lats[i]) +  89.75) / 179.5 * 359 + 2
        TT[i] = itp(idx_lon, idx_lat)
    end
    return TT
end
function lunar_surface_temperatures_DIVINER_avg(lons::AbstractVector, lats::AbstractVector)
    S = typeof(promote(lons[1], lats[1])[1])
    return lunar_surface_temperatures_DIVINER_avg(S.(lons), S.(lats))
end
function lunar_surface_temperatures_DIVINER_avg(lons::AbstractVector{<:Integer},
                                                lats::AbstractVector{<:Integer})
    S = promote_type(typeof(lons[1]), typeof(lats[1]), Float64)
    return lunar_surface_temperatures_DIVINER_avg(S.(lons), S.(lats))
end
function lunar_surface_temperatures_DIVINER_avg(x)
    return lunar_surface_temperatures_DIVINER_avg(_gettheta(x), _getphi(x))
end
function lunar_surface_temperatures_DIVINER_avg(X::AbstractVector)
    lons, lats = [_gettheta(x) for x in X], [_getphi(x) for x in X]
    return lunar_surface_temperatures_DIVINER_avg(lons, lats)
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
