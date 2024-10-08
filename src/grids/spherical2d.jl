############################################################################################
#::. TYPES
############################################################################################
abstract type AbstractSpherical2DGrid <: AbstractSphericalGrid end


############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    [1] struct Spherical2DGrid{T<:AbstractFloat} <: AbstractSpherical2DGrid end
    [2] Spherical2DGrid([T::Type,] r::Real, N_lon::Integer, N_lat::Integer; kwargs...)

Global structured grid of surface coordinates (2D) of type `GlobalSphericalPosition{T}` over
a sphere of radius `r`.

| Field      | Type; with `T<:AbstractFloat`        | Description                       |
|:---------- |:------------------------------------ |:--------------------------------- |
| `r`        | `T`                                  | radius of global body (sphere)    |
| `N_lon`    | `Int64`                              | # elements in longitude direction |
| `N_lat`    | `Int64`                              | # elements in latitude direction  |
| `coords`   | `Vector{GlobalSphericalPosition{T}}` | coordinates                       |
| `areas`    | `Vector{T}`                          | surface area                      |
| `lonrange` | `Tuple{T, T}`                        | longitude range                   |
| `latrange` | `Tuple{T, T}`                        | latitude range                    |

To ensure expected behavior, the grid object should generally be created with the outer
constructor [2].
"""
struct Spherical2DGrid{T<:AbstractFloat} <: AbstractSpherical2DGrid
    r::T
    N_lon::Integer
    N_lat::Integer
    coords::Vector{GlobalSphericalPosition{T}}
    areas::Vector{T}
    lonrange::Tuple{T, T}
    latrange::Tuple{T, T}
end
function Spherical2DGrid(T::Type{<:AbstractFloat}, r::Real, N_lon::Integer,
                         N_lat::Integer; lonrange=(-pi, pi), latrange=(-pi/2, pi/2))

    @assert lonrange[1] <= lonrange[2] "Invalid longitude range"
    @assert latrange[1] <= latrange[2] "Invalid latitude range"

    lon = T.(repeat(range(1/2/N_lon, 1 - 1/2/N_lon, length=N_lon) *
        (lonrange[2]-lonrange[1]) .+ lonrange[1]; inner=N_lat))
    lat = T.(repeat(range(1/2/N_lat, 1 - 1/2/N_lat, length=N_lat) *
        (latrange[2]-latrange[1]) .+ latrange[1]; outer=N_lon))

    dlon, dlat = lon[N_lat+1] - lon[1], lat[2] - lat[1]
    areas = [r^2 * dlon * (sin(p+dlat/2) - sin(p-dlat/2)) for p in lat]

    return Spherical2DGrid{T}(T(r), N_lon, N_lat,
                              GlobalSphericalPosition.(T(r), lon, lat),
                              areas, T.(lonrange), T.(latrange))
end
function Spherical2DGrid(r::Real, N_lon::Integer, N_lat::Integer; kwargs...)
    return Spherical2DGrid(typeof(r), r, N_lon, N_lat; kwargs...)
end
function Spherical2DGrid(r::Integer, N_lon::Integer, N_lat::Integer; kwargs...)
    return Spherical2DGrid(promote_type(typeof(r), Float64), r, N_lon, N_lat; kwargs...)
end


"""
    [1] struct Spherical2DGrid_EqSim{T<:AbstractFloat} <: AbstractSpherical2DGrid end
    [2] Spherical2DGrid_EqSim([T::Type,] r::Real, N_lon::Integer, N_lat::Integer; kwargs...)

Global structured grid of surface coordinates (2D) of type `GlobalSphericalPosition{T}` over
the upper hemisphere with radius `r`, assuming equatorial symmetry.

| Field      | Type; with `T<:AbstractFloat`        | Description                       |
|:---------- |:------------------------------------ |:--------------------------------- |
| `r`        | `T`                                  | radius of global body (sphere)    |
| `N_lon`    | `Int64`                              | # elements in longitude direction |
| `N_lat`    | `Int64`                              | # elements in latitude direction  |
| `coords`   | `Vector{GlobalSphericalPosition{T}}` | coordinates                       |
| `areas`    | `Vector{T}`                          | surface area                      |
| `lonrange` | `Tuple{T, T}`                        | longitude range                   |
| `latrange` | `Tuple{T, T}`                        | latitude range                    |

To ensure expected behavior, the grid object should generally be created with the outer
constructor [2]. The function comes with a unique key-word argument `latmax` to set the
maximum latitude of the grid. Due to its equatorial symmetry, `latrange[1]=0` and `latmax>0`.
"""
struct Spherical2DGrid_EqSim{T<:AbstractFloat} <: AbstractSpherical2DGrid
    r::T
    N_lon::Integer
    N_lat::Integer
    coords::Vector{GlobalSphericalPosition{T}}
    areas::Vector{T}
    lonrange::Tuple{T, T}
    latrange::Tuple{T, T}
end
function Spherical2DGrid_EqSim(T::Type{<:AbstractFloat}, r::Real, N_lon::Integer,
                               N_lat::Integer; lonrange=(-pi, pi), latmax=pi/2)

    @assert lonrange[1] <= lonrange[2] "Invalid longitude range"
    @assert latmax >= 0 "Invalid `latmax`"

    lon = T.(repeat(range(1/2/N_lon, 1 - 1/2/N_lon, length=N_lon) *
        (lonrange[2]-lonrange[1]) .+ lonrange[1]; inner=N_lat))
    lat = T.(repeat(range(1/2/N_lat, 1 - 1/2/N_lat, length=N_lat) * latmax; outer=N_lon))

    dlon, dlat = lon[N_lat+1] - lon[1], lat[2] - lat[1]
    areas = [r^2 * dlon * (sin(p+dlat/2) - sin(p-dlat/2)) for p in lat]

    return Spherical2DGrid_EqSim{T}(T(r), N_lon, N_lat,
                                    GlobalSphericalPosition.(T(r), lon, lat),
                                    areas, T.(lonrange), T.((0, latmax)))
end
function Spherical2DGrid_EqSim(r::Real, N_lon::Integer, N_lat::Integer; kwargs...)
    return Spherical2DGrid_EqSim(typeof(r), r, N_lon, N_lat; kwargs...)
end
function Spherical2DGrid_EqSim(r::Integer, N_lon::Integer, N_lat::Integer; kwargs...)
    return Spherical2DGrid_EqSim(promote_type(typeof(r), Float64), r, N_lon, N_lat; kwargs...)
end


"""
    [1] struct Spherical2DGrid_Reduced{T<:AbstractFloat} <: AbstractSpherical2DGrid end
    [2] Spherical2DGrid_Reduced([T::Type,] r::Real, N_lat::Integer; kwargs...)

Global structured grid of surface coordinates (2D) of type `GlobalSphericalPosition{T}` over
the sphere radius `r`. The grid is reduced in the longitude direction to have approximately
equal `(lonrange[2]-lonrange[1])*r*cos(lat)/N_lon` grid element lengths.

| Field      | Type; with `T<:AbstractFloat`        | Description                       |
|:---------- |:------------------------------------ |:--------------------------------- |
| `r`        | `T`                                  | radius of global body (sphere)    |
| `N_lon`    | `Vector{Int64}`                      | # elements in longitude direction |
| `N_lat`    | `Int64`                              | # elements in latitude direction  |
| `coords`   | `Vector{GlobalSphericalPosition{T}}` | coordinates                       |
| `areas`    | `Vector{T}`                          | surface area                      |
| `lonrange` | `Tuple{T, T}`                        | longitude range                   |
| `latrange` | `Tuple{T, T}`                        | latitude range                    |

To ensure expected behavior, the grid object should generally be created with the outer
constructor [2].
"""
struct Spherical2DGrid_Reduced{T<:AbstractFloat} <: AbstractSpherical2DGrid
    r::T
    N_lon::Vector{Int64}
    N_lat::Integer
    coords::Vector{GlobalSphericalPosition{T}}
    areas::Vector{T}
    lonrange::Tuple{T, T}
    latrange::Tuple{T, T}
end
function Spherical2DGrid_Reduced(T::Type{<:AbstractFloat}, r::Real, N_lat::Integer;
                                 lonrange=(-pi, pi), latrange=(-pi/2, pi/2))

    @assert lonrange[1] <= lonrange[2] "Invalid longitude range"
    @assert latrange[1] <= latrange[2] "Invalid latitude range"

    # calculate maximum longitude size and latitude grid
    dlon_max = T((latrange[2] - latrange[1])/N_lat)
    lat0 = T.(range(1/2/N_lat, 1 - 1/2/N_lat, length=N_lat) *
        (latrange[2]-latrange[1]) .+ latrange[1])

    # claculate coordinates for longitudes and latitudes
    lon, lat, N_lon = T[], T[], Int64[]
    for p0 in lat0
        n_lon = ceil(Int64, (lonrange[2] - lonrange[1])*cos(p0)/dlon_max + eps(T))
        push!(lon, range(1/2/n_lon, 1 - 1/2/n_lon, length=n_lon) *
            (lonrange[2]-lonrange[1]) .+ lonrange[1]...)
        push!(lat, repeat([p0], n_lon)...)
        push!(N_lon, n_lon)
    end

    # areas calculation
    areas = T[]
    for i in eachindex(N_lon)
        dlon, dlat = (lonrange[2] - lonrange[1]) / N_lon[i], lat0[2] - lat0[1]
        push!(areas, repeat(
            [r^2 * dlon * (sin(lat0[i]+dlat/2) - sin(lat0[i]-dlat/2))], N_lon[i])...)
    end

    return Spherical2DGrid_Reduced(T(r), N_lon, N_lat,
                                   GlobalSphericalPosition.(T(r), lon, lat),
                                   areas, T.(lonrange), T.(latrange))
end
function Spherical2DGrid_Reduced(r::Real, N_lat::Integer; kwargs...)
    return Spherical2DGrid_Reduced(typeof(r), r, N_lat; kwargs...)
end
function Spherical2DGrid_Reduced(r::Integer, N_lat::Integer; kwargs...)
    return Spherical2DGrid_Reduced(promote_type(typeof(r), Float64), r, N_lat; kwargs...)
end


"""
    [1] struct Spherical2DGrid_Reduced_EqSim{T<:AbstractFloat} <: AbstractSpherical2DGrid end
    [2] Spherical2DGrid_Reduced_EqSim([T::Type,] r::Real, N_lat::Integer; kwargs...)

Global structured grid of surface coordinates (2D) of type `GlobalSphericalPosition{T}` over
the upper hemisphere with radius `r`, assuming equatorial symmetry. The grid is reduced in
the longitude direction to have approximately equal
`(lonrange[2]-lonrange[1])*r*cos(lat)/N_lon` grid element lengths.

| Field     | Type; with `T<:AbstractFloat`        | Description                       |
|:--------- |:------------------------------------ |:--------------------------------- |
| `r`       | `T`                                  | radius of global body (sphere)    |
| `N_lon`   | `Vector{Int64}`                      | # elements in longitude direction |
| `N_lat`   | `Int64`                              | # elements in latitude direction  |
| `coords`  | `Vector{GlobalSphericalPosition{T}}` | coordinates                       |
| `areas`   | `Vector{T}`                          | surface area                      |

To ensure expected behavior, the grid object should generally be created with the outer
constructor [2]. The function comes with a unique key-word argument `latmax` to set the
maximum latitude of the grid. Due to its equatorial symmetry, `latrange[1]=0` and `latmax>0`.
"""
struct Spherical2DGrid_Reduced_EqSim{T<:AbstractFloat} <: AbstractSpherical2DGrid
    r::T
    N_lon::Vector{Int64}
    N_lat::Integer
    coords::Vector{GlobalSphericalPosition{T}}
    areas::Vector{T}
    lonrange::Tuple{T, T}
    latrange::Tuple{T, T}
end
function Spherical2DGrid_Reduced_EqSim(T::Type{<:AbstractFloat}, r::Real, N_lat::Integer;
                                       lonrange=(-pi, pi), latmax=pi/2)

    @assert lonrange[1] <= lonrange[2] "Invalid longitude range"
    @assert latmax >= 0 "Invalid `latmax`"

    dlon_max = T(latmax/N_lat)
    lat0 = T.(range(1/2/N_lat, 1 - 1/2/N_lat, length=N_lat) * latmax)

    lon, lat, N_lon = T[], T[], Int64[]
    for p0 in lat0
        n_lon = ceil(Int64, (lonrange[2] - lonrange[1])*cos(p0)/dlon_max + eps(T))
        push!(lon, range(1/2/n_lon, 1 - 1/2/n_lon, length=n_lon) *
            (lonrange[2]-lonrange[1]) .+ lonrange[1]...)
        push!(lat, repeat([p0], n_lon)...)
        push!(N_lon, n_lon)
    end

    areas = T[]
    for i in eachindex(N_lon)
        # idx = accumulate(+, N_lon[1:i])[end]
        # dlon, dlat = lon[idx] - lon[idx-1], lat0[2] - lat0[1]
        dlon, dlat = (lonrange[2] - lonrange[1]) / N_lon[i], lat0[2] - lat0[1]
        push!(areas, repeat([r^2 * dlon * (sin(lat0[i]+dlat/2) - sin(lat0[i]-dlat/2))], N_lon[i])...)
    end

    return Spherical2DGrid_Reduced_EqSim(T(r), N_lon, N_lat,
                                         GlobalSphericalPosition.(T(r), lon, lat),
                                         areas, T.(lonrange), T.((0, latmax)))
end
function Spherical2DGrid_Reduced_EqSim(r::Real, N_lat::Integer; kwargs...)
    return Spherical2DGrid_Reduced_EqSim(typeof(r), r, N_lat; kwargs...)
end
function Spherical2DGrid_Reduced_EqSim(r::Integer, N_lat::Integer; kwargs...)
    return Spherical2DGrid_Reduced_EqSim(promote_type(typeof(r), Float64), r, N_lat; kwargs...)
end


############################################################################################
#::. UTILITY FUNCTIONS
############################################################################################
function coord2idx(grid::Spherical2DGrid, lon::T, lat::T)::Int64 where {T<:AbstractFloat}
    PI  = T(pi) # to prevent numerical cutoff/rounding issues
    lon = pclamp(lon, -PI, PI)
    lat = clamp(lat, -PI/2, PI/2)

    lonrange, latrange = grid.lonrange, grid.latrange
    if !(lonrange[1] <= lon <= lonrange[2]); return 0; end
    if !(latrange[1] <= lat <= latrange[2]); return 0; end

    idxlon = max(1,ceil(Int64, (lon-lonrange[1])/(lonrange[2]-lonrange[1])*grid.N_lon))
    idxlat = max(1,ceil(Int64, (lat-latrange[1])/(latrange[2]-latrange[1])*grid.N_lat))
    return (idxlon-1) * grid.N_lat + idxlat
end
function coord2idx(grid::Spherical2DGrid_EqSim, lon::T, lat::T)::Int64 where {T<:AbstractFloat}
    PI  = T(pi) # to prevent numerical cutoff/rounding issues
    lon = pclamp(lon, -PI, PI)
    lat = clamp(lat, -PI/2, PI/2)

    lonrange, latrange = grid.lonrange, grid.latrange
    if !(lonrange[1] <= lon <= lonrange[2]); return 0; end
    if !(-latrange[2] <= lat <= latrange[2]); return 0; end

    idxlon = max(1,ceil(Int64, (lon-lonrange[1])/(lonrange[2]-lonrange[1])*grid.N_lon))
    idxlat = max(1,ceil(Int64, abs(lat)/latrange[2]*grid.N_lat))
    return (idxlon-1) * grid.N_lat + idxlat
end
function coord2idx(grid::Spherical2DGrid_Reduced, lon::T, lat::T)::Int64 where {T<:AbstractFloat}
    PI  = T(pi) # to prevent numerical cutoff/rounding issues
    lon = pclamp(lon, -PI, PI)
    lat = clamp(lat, -PI/2, PI/2)

    lonrange, latrange = grid.lonrange, grid.latrange
    if !(lonrange[1] <= lon <= lonrange[2]); return 0; end
    if !(latrange[1] <= lat <= latrange[2]); return 0; end

    idxlat = max(1,ceil(Int64, (lat-latrange[1])/(latrange[2]-latrange[1])*grid.N_lat))
    idxlon = max(1,ceil(Int64, (lon-lonrange[1])/(lonrange[2]-lonrange[1])*grid.N_lon[idxlat]))
    if idxlat == 1; return idxlon; end
    return idxlon + accumulate(+, grid.N_lon[1:idxlat])[end-1]
end
function coord2idx(grid::Spherical2DGrid_Reduced_EqSim, lon::T, lat::T)::Int64 where {T<:AbstractFloat}
    PI  = T(pi) # to prevent numerical cutoff/rounding issues
    lon = pclamp(lon, -PI, PI)
    lat = clamp(lat, -PI/2, PI/2)

    lonrange, latrange = grid.lonrange, grid.latrange
    if !(lonrange[1] <= lon <= lonrange[2]); return 0; end
    if !(-latrange[2] <= lat <= latrange[2]); return 0; end

    idxlat = max(1,ceil(Int64, abs(lat)/latrange[2]*grid.N_lat))
    idxlon = max(1,ceil(Int64, (lon-lonrange[1])/(lonrange[2]-lonrange[1])*grid.N_lon[idxlat]))
    if idxlat == 1; return idxlon; end
    return idxlon + accumulate(+, grid.N_lon[1:idxlat])[end-1]
end
function coord2idx(grid::AbstractSpherical2DGrid, lon::Real, lat::Real)::Int64
    return coord2idx(grid, promote(lon, lat)...)
end
function coord2idx(grid::AbstractSpherical2DGrid, lon::Integer, lat::Integer)::Int64
    return coord2idx(grid, float(lon), float(lat))
end
function coord2idx(grid::AbstractSpherical2DGrid, r::Real, lon::Real, lat::Real)::Int64
    return coord2idx(grid, lon, lat)
end


surfacecoords(grid::AbstractSpherical2DGrid) = coords(grid)


volumes(grid::AbstractSpherical2DGrid) = zeros(typeof(grid.r), length(grid.coords))


############################################################################################
#::. EXTENSIONS
############################################################################################
Base.length(grid::AbstractSpherical2DGrid) = length(grid.coords)


Base.size(grid::AbstractSpherical2DGrid) = (grid.N_lon, grid.N_lat)


Base.show(io::IO, ::MIME"text/plain", grid::AbstractSpherical2DGrid) =
    print(io, "$(typeof(grid)):\n"*
            " r:        $(grid.r)\n"*
            " lonrange: $(grid.lonrange)\n"*
            " latrange: $(grid.latrange)\n"*
            " N_lon:    $(grid.N_lon)\n"*
            " N_lat:    $(grid.N_lat)\n"*
            " coords:   $(length(grid.coords))-element $(typeof(coords(grid)))\n"*
            " areas:    $(length(grid.areas))-element $(typeof(grid.areas))\n"*
            "             min: $(min(grid.areas...))\n"*
            "             max: $(max(grid.areas...))\n")

Base.show(io::IO, ::MIME"text/plain", grid::Union{Spherical2DGrid_Reduced, Spherical2DGrid_Reduced_EqSim}) =
    print(io, "$(typeof(grid)):\n"*
            " r:        $(grid.r)\n"*
            " lonrange: $(grid.lonrange)\n"*
            " latrange: $(grid.latrange)\n"*
            " N_lon:    $(length(grid.N_lon))-element Vector{Int64}\n"*
            "             [$(grid.N_lon[1]) .. $(grid.N_lon[end])]\n"*
            " N_lat:    $(grid.N_lat)\n"*
            " coords:   $(length(grid.coords))-element $(typeof(coords(grid)))\n"*
            " areas:    $(length(grid.areas))-element $(typeof(grid.areas))\n"*
            "             min: $(min(grid.areas...))\n"*
            "             max: $(max(grid.areas...))\n")


############################################################################################
#::. EXPORTS
############################################################################################
export
    Spherical2DGrid,
    Spherical2DGrid_EqSim,
    Spherical2DGrid_Reduced,
    Spherical2DGrid_Reduced_EqSim
