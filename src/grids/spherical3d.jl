############################################################################################
#::. TYPES
############################################################################################
abstract type AbstractSpherical3DGrid <: AbstractSphericalGrid end


############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    struct Spherical3DGrid{T} <: AbstractSpherical3DGrid
    Spherical3DGrid([T], r0, h, N_lon, N_lat; kwargs...)
    Spherical3DGrid([T], r, N_lon, N_lat; kwargs...)

Global structured volume grid (3D) of type `GlobalSphericalPosition{T}` over a sphere of
radius `r0` with heights of the individual radial layers `h`. Alternatively, the radial
distances can directly be specified by `r` (`r = r0 .+ h`). Note that all heights have
to be positive. Note that the radial position of each individual grid element is located
at the bottom (base/surface) of each cell.

To ensure expected behavior, the grid object should generally be created with the outer
constructors. Those functions take the key-word `lonrange` and `latrange` to specify the
range of the longitude and latitude. The default values are `(-pi, pi)` and `(-pi/2, pi/2)`,
respectively.

The optional argument `T can be used to control the type of the structs fields.

# Fields & Arguments

| Field      | Type; with `T<:AbstractFloat`        | Description                       |
|:---------- |:------------------------------------ |:--------------------------------- |
| `r0`       | `T`                                  | radius of global body (sphere)    |
| `h`        | `Vector{T}`                          | heights above radial base `r0`    |
| `N_r`      | `Int64`                              | # elements in radial direction    |
| `N_lon`    | `Int64`                              | # elements in longitude direction |
| `N_lat`    | `Int64`                              | # elements in laditude direction  |
| `coords`   | `Vector{GlobalSphericalPosition{T}}` | coordinates                       |
| `areas`    | `Vector{T}`                          | surface area                      |
| `volumes`  | `Vector{T}`                          | volumes                           |
| `lonrange` | `Tuple{T, T}`                        | longitude range                   |
| `latrange` | `Tuple{T, T}`                        | latitude range                    |
"""
struct Spherical3DGrid{T<:AbstractFloat} <: AbstractSpherical3DGrid
    r0::T
    h::Vector{T}
    N_r::Int64
    N_lon::Int64
    N_lat::Int64
    coords::Vector{GlobalSphericalPosition{T}}
    areas::Vector{T}
    volumes::Vector{T}
    lonrange::Tuple{T, T}
    latrange::Tuple{T, T}
end
function Spherical3DGrid(T::Type{<:AbstractFloat}, r0::Real, h::AbstractVector,
            N_lon::Integer, N_lat::Integer; lonrange=(-pi, pi), latrange=(-pi/2, pi/2))

    # input assertions
    @assert sort(h)[1] > 0 "Invalid heights"
    @assert lonrange[1] <= lonrange[2] "Invalid longitude range"
    @assert latrange[1] <= latrange[2] "Invalid latitude range"

    # coordinates
    h, N_r = sort(h), length(h)
    r = T.(repeat(r0 .+ vcat(0, h[1:end-1]), inner=N_lon*N_lat))
    lon = T.(repeat(range(1/2/N_lon, 1 - 1/2/N_lon, length=N_lon) *
        (lonrange[2]-lonrange[1]) .+ lonrange[1]; inner=N_lat, outer=N_r))
    lat = T.(repeat(range(1/2/N_lat, 1 - 1/2/N_lat, length=N_lat) *
        (latrange[2]-latrange[1]) .+ latrange[1]; outer=N_lon*N_r))

    # area and volume calculation
    dR = vcat(h[1], diff(h))
    dlon, dlat = (lonrange[2] - lonrange[1]) / N_lon, (latrange[2] - latrange[1]) / N_lat
    areas, volumes = T[], T[]
    for i in 1:N_r*N_lon*N_lat
        ri, dr = r[i], dR[Int64(1 + floor((i-1) / (N_lon*N_lat)))]
        push!(areas, ri^2 * dlon * (sin(lat[i]+dlat/2) - sin(lat[i]-dlat/2)))
        push!(volumes, ((ri + dr)^3 - ri^3)/3 * dlon * (sin(lat[i]+dlat/2) - sin(lat[i]-dlat/2)))
    end

    return Spherical3DGrid(T(r0), T.(h), N_r, N_lon, N_lat,
                           GlobalSphericalPosition.(r, lon, lat),
                           T.(areas), T.(volumes), T.(lonrange), T.(latrange))
end
function Spherical3DGrid(r0::Real, h::AbstractVector, N_lon::Integer, N_lat::Integer;
                         kwargs...)
    return Spherical3DGrid(typeof(r0), r0, h, N_lon, N_lat; kwargs...)
end
function Spherical3DGrid(r0::Integer, h::AbstractVector, N_lon::Integer, N_lat::Integer;
                         kwargs...)
    return Spherical3DGrid(promote_type(typeof(r0),Float64), r0, h, N_lon, N_lat; kwargs...)
end
function Spherical3DGrid(T::Type{<:AbstractFloat}, r::AbstractVector, N_lon::Integer,
                         N_lat::Integer; kwargs...)
    r0, h = r[1], accumulate(+,diff(r))
    return Spherical3DGrid(T, r0, h, N_lon, N_lat; kwargs...)
end
function Spherical3DGrid(r::AbstractVector, N_lon::Integer, N_lat::Integer; kwargs...)
    return Spherical3DGrid(typeof(r[1]), r, N_lon, N_lat; kwargs...)
end
function Spherical3DGrid(r::AbstractVector{<:Integer}, N_lon::Integer, N_lat::Integer;
                         kwargs...)
    return Spherical3DGrid(promote_type(typeof(r[1]), Float64), r, N_lon, N_lat; kwargs...)
end

"""
    struct Spherical3DGrid_EqSim{T} <: AbstractSpherical3DGrid
    Spherical3DGrid_EqSim([T], r0, h, N_lon, N_lat; kwargs...)
    Spherical3DGrid_EqSim([T], r, N_lon, N_lat; kwargs...)

Global structured volume grid (3D) of type `GlobalSphericalPosition{T}` over a hemisphere of
radius `r0` with heights of the individual radial layers `h`. Alternatively, the radial
distances can directly be specified by `r` (`r = r0 .+ h`). Note that all heights have
to be positive. Note that the radial position of each individual grid element is located
at the bottom (base/surface) of each cell.

To ensure expected behavior, the grid object should generally be created with the outer
constructors. Those functions take the key-word `lonrange` and `latmax` to
specify the range of the longitude and latitude. The default values are `(-pi, pi)` and
`pi/2`, respectively (the lower boundary of latitudes for equatorially symmetric grids is
always zero).

The optional argument `T can be used to control the type of the structs fields.

# Fields & Arguments

| Field      | Type; with `T<:AbstractFloat`        | Description                       |
|:---------- |:------------------------------------ |:--------------------------------- |
| `r0`       | `T`                                  | radius of global body (sphere)    |
| `h`        | `Vector{T}`                          | heights above radial base `r0`    |
| `N_r`      | `Int64`                              | # elements in radial direction    |
| `N_lon`    | `Int64`                              | # elements in azimuth direction   |
| `N_lat`    | `Int64`                              | # elements in elevation direction |
| `coords`   | `Vector{GlobalSphericalPosition{T}}` | coordinates                       |
| `areas`    | `Vector{T}`                          | surface area                      |
| `volumes`  | `Vector{T}`                          | volumes                           |
| `lonrange` | `Tuple{T, T}`                        | longitude range                   |
| `latrange` | `Tuple{T, T}`                        | latitude range                    |
"""
struct Spherical3DGrid_EqSim{T<:AbstractFloat} <: AbstractSpherical3DGrid
    r0::T
    h::Vector{T}
    N_r::Int64
    N_lon::Int64
    N_lat::Int64
    coords::Vector{GlobalSphericalPosition{T}}
    areas::Vector{T}
    volumes::Vector{T}
    lonrange::Tuple{T, T}
    latrange::Tuple{T, T}
end
function Spherical3DGrid_EqSim(T::Type{<:AbstractFloat}, r0::Real, h::AbstractVector,
            N_lon::Integer, N_lat::Integer; lonrange=(-pi, pi), latmax=pi/2)

    # input assertions
    @assert sort(h)[1] > 0 "Invalid heights"
    @assert lonrange[1] <= lonrange[2] "Invalid longitude range"
    @assert latmax >= 0 "Invalid `latmax`"

    # coordinates
    h, N_r = sort(h), length(h)
    r = T.(repeat(r0 .+ vcat(0, h[1:end-1]), inner=N_lon*N_lat))
    lon = T.(repeat(range(1/2/N_lon, 1-1/2/N_lon, length=N_lon) *
        (lonrange[2]-lonrange[1]) .+ lonrange[1]; inner=N_lat, outer=N_r))
    lat = T.(repeat(range(1/2/N_lat, 1-1/2/N_lat, length=N_lat) * latmax; outer=N_lon*N_r))

    # area and volume calculation
    dR, dlon, dlat = vcat(h[1], diff(h)), (lonrange[2] - lonrange[1]) / N_lon, latmax / N_lat
    areas, volumes = T[], T[]
    for i in 1:N_r*N_lon*N_lat
        ri, dr = r[i], dR[Int64(1 + floor((i-1) / (N_lon*N_lat)))]
        push!(areas, ri^2 * dlon * (sin(lat[i]+dlat/2) - sin(lat[i]-dlat/2)))
        push!(volumes, ((ri + dr)^3 - ri^3)/3 * dlon * (sin(lat[i]+dlat/2) - sin(lat[i]-dlat/2)))
    end

    return Spherical3DGrid_EqSim(T(r0), T.(h), N_r, N_lon, N_lat,
                                 GlobalSphericalPosition.(r, lon, lat),
                                 T.(areas), T.(volumes), T.(lonrange), T.((0, latmax)))
end
function Spherical3DGrid_EqSim(r0::Real, h::AbstractVector, N_lon::Integer, N_lat::Integer;
                               kwargs...)
    return Spherical3DGrid_EqSim(typeof(r0), r0, h, N_lon, N_lat; kwargs...)
end
function Spherical3DGrid_EqSim(r0::Integer, h::AbstractVector, N_lon::Integer,
                               N_lat::Integer; kwargs...)
    return Spherical3DGrid_EqSim(promote_type(typeof(r0),Float64),r0,h,N_lon,N_lat;kwargs...)
end
function Spherical3DGrid_EqSim(T::Type{<:AbstractFloat}, r::AbstractVector, N_lon::Integer,
                               N_lat::Integer; kwargs...)
    r0, h = r[1], accumulate(+,diff(r))
    return Spherical3DGrid_EqSim(T, r0, h, N_lon, N_lat; kwargs...)
end
function Spherical3DGrid_EqSim(r::AbstractVector, N_lon::Integer, N_lat::Integer; kwargs...)
    return Spherical3DGrid_EqSim(typeof(r[1]), r, N_lon, N_lat; kwargs...)
end
function Spherical3DGrid_EqSim(r::AbstractVector{<:Integer}, N_lon::Integer, N_lat::Integer;
                               kwargs...)
    return Spherical3DGrid_EqSim(promote_type(typeof(r[1]),Float64),r,N_lon,N_lat;kwargs...)
end


"""
    struct Spherical3DGrid_Reduced{T} <: AbstractSpherical3DGrid
    Spherical3DGrid_Reduced([T], r0, h, N_lat; kwargs...)
    Spherical3DGrid_Reduced([T], r, N_lat; kwargs...)

Global structured volume grid (3D) of type `GlobalSphericalPosition{T}` over a sphere of
radius `r0` with heights of the individual radial layers `h`. Alternatively, the radial
distances can directly be specified by `r` (`r = r0 .+ h`). Note that all heights have
to be positive. The grid is reduced in the azimuth direction to have approximately equal
`(lonrange[2]-lonrange[1])*r*cos(lat)/N_lon` grid element lengths.

To ensure expected behavior, the grid object should generally be created with the outer
constructors. Those functions take the key-word `lonrange` and `latrange` to specify the
range of the longitude and latitude. The default values are `(-pi, pi)` and `(-pi/2, pi/2)`,
respectively.

The optional argument `T can be used to control the type of the structs fields.

# Fields & Arguments

| Field      | Type; with `T<:AbstractFloat`        | Description                       |
|:---------- |:------------------------------------ |:--------------------------------- |
| `r0`       | `T`                                  | radius of global body (sphere)    |
| `h`        | `Vector{T}`                          | heights above radial base `r0`    |
| `N_r`      | `Int64`                              | # elements in radial direction    |
| `N_lon`    | `Int64`                              | # elements in azimuth direction   |
| `N_lat`    | `Int64`                              | # elements in elevation direction |
| `coords`   | `Vector{GlobalSphericalPosition{T}}` | coordinates                       |
| `areas`    | `Vector{T}`                          | surface area                      |
| `volumes`  | `Vector{T}`                          | volumes                           |
| `lonrange` | `Tuple{T, T}`                        | longitude range                   |
| `latrange` | `Tuple{T, T}`                        | latitude range                    |
"""
struct Spherical3DGrid_Reduced{T<:AbstractFloat} <: AbstractSpherical3DGrid
    r0::T
    h::Vector{T}
    N_r::Int64
    N_lon::Vector{Int64}
    N_lat::Int64
    coords::Vector{GlobalSphericalPosition{T}}
    areas::Vector{T}
    volumes::Vector{T}
    lonrange::Tuple{T, T}
    latrange::Tuple{T, T}
end
function Spherical3DGrid_Reduced(T::Type{<:AbstractFloat}, r0::Real, h::AbstractVector,
                                 N_lat::Integer; lonrange=(-pi, pi), latrange=(-pi/2, pi/2))

    # input assertions
    @assert sort(h)[1] > 0 "Invalid heights"
    @assert lonrange[1] <= lonrange[2] "Invalid longitude range"
    @assert latrange[1] <= latrange[2] "Invalid latitude range"

    # prepare coordiante calculation
    h, N_r = sort(h), length(h)
    dlat = dlon_max = T((latrange[2] - latrange[1])/N_lat)
    lat0 = T.(range(1/2/N_lat, 1 - 1/2/N_lat, length=N_lat) *
        (latrange[2]-latrange[1]) .+ latrange[1])

    # coordinates
    r, lon, lat, N_lon = T[], T[], T[], Int64[]
    for p0 in lat0
        n_lon = ceil(Int64, (lonrange[2] - lonrange[1])*cos(p0)/dlon_max + eps(T))
        push!(lon, range(1/2/n_lon, 1 - 1/2/n_lon, length=n_lon) *
            (lonrange[2]-lonrange[1]) .+ lonrange[1]...)
        push!(lat, repeat([p0], n_lon)...)
        push!(N_lon, n_lon)
    end
    r = T.(repeat(r0 .+ vcat(0, h[1:end-1]), inner=sum(N_lon)))
    lon = repeat(lon, outer=N_r)
    lat = repeat(lat, outer=N_r)

    # area and volume calculation
    areas, volumes = T[], T[]
    for i in 1:N_r, j in 1:N_lat, _ in 1:N_lon[j]
        ri, dr = [r0 .+ vcat(0, h[1:end-1])...][i], vcat(h[1], diff(h))[i]
        dlon = (lonrange[2] - lonrange[1]) / N_lon[j]
        push!(areas, ri^2 * dlon * (sin(lat0[j]+dlat/2) - sin(lat0[j]-dlat/2)))
        push!(volumes, ((ri + dr)^3 - ri^3)/3 * dlon * (sin(lat0[j]+dlat/2) - sin(lat0[j]-dlat/2)))
    end

    return Spherical3DGrid_Reduced(T(r0), T.(h), N_r, N_lon, N_lat,
                                   GlobalSphericalPosition.(r, lon, lat),
                                   T.(areas), T.(volumes), T.(lonrange), T.(latrange))
end
function Spherical3DGrid_Reduced(r0::Real, h::AbstractVector, N_lat::Integer; kwargs...)
    return Spherical3DGrid_Reduced(typeof(r0), r0, h, N_lat; kwargs...)
end
function Spherical3DGrid_Reduced(r0::Integer, h::AbstractVector, N_lat::Integer; kwargs...)
    return Spherical3DGrid_Reduced(promote_type(typeof(r0), Float64), r0, h, N_lat; kwargs...)
end
function Spherical3DGrid_Reduced(T::Type{<:AbstractFloat}, r::AbstractVector, N_lat::Integer; kwargs...)
    r0, h = r[1], accumulate(+,diff(r))
    Spherical3DGrid_Reduced(T, r0, h, N_lat; kwargs...)
end
function Spherical3DGrid_Reduced(r::AbstractVector, N_lat::Integer; kwargs...)
    return Spherical3DGrid_Reduced(typeof(r[1]), r, N_lat; kwargs...)
end
function Spherical3DGrid_Reduced(r::AbstractVector{<:Integer}, N_lat::Integer; kwargs...)
    return Spherical3DGrid_Reduced(promote_type(typeof(r[1]), Float64), r, N_lat; kwargs...)
end


"""
    struct Spherical3DGrid_Reduced_EqSim{T} <: AbstractSpherical3DGrid
    Spherical3DGrid_Reduced_EqSim([T], r0, h, N_lat)
    Spherical3DGrid_Reduced_EqSim([T], r, N_lat)

Global structured volume grid (3D) of type `GlobalSphericalPosition{T}` over a sphere of
radius `r0` with heights of the individual radial layers `h`. Alternatively, the radial
distances can directly be specified by `r` (`r = r0 .+ h`). Note that all heights have
to be positive. The grid is reduced in the azimuth direction to have approximately equal
`2*pi*r*cos(lat)/N_lon` grid element lengths.

| Field      | Type; with `T<:AbstractFloat`        | Description                       |
|:---------- |:------------------------------------ |:--------------------------------- |
| `r0`       | `T`                                  | radius of global body (sphere)    |
| `h`        | `Vector{T}`                          | heights above radial base `r0`    |
| `N_r`      | `Int64`                              | # elements in radial direction    |
| `N_lon`    | `Int64`                              | # elements in azimuth direction   |
| `N_lat`    | `Int64`                              | # elements in elevation direction |
| `coords`   | `Vector{GlobalSphericalPosition{T}}` | coordinates                       |
| `areas`    | `Vector{T}`                          | surface area                      |
| `volumes`  | `Vector{T}`                          | volumes                           |
| `lonrange` | `Tuple{T, T}`                        | longitude range                   |
| `latrange` | `Tuple{T, T}`                        | latitude range                    |

To ensure expected behavior, the grid object should generally be created with the outer
constructors [2] or [3].
"""
struct Spherical3DGrid_Reduced_EqSim{T<:AbstractFloat} <: AbstractSpherical3DGrid
    r0::T
    h::Vector{T}
    N_r::Int64
    N_lon::Vector{Int64}
    N_lat::Int64
    coords::Vector{GlobalSphericalPosition{T}}
    areas::Vector{T}
    volumes::Vector{T}
    lonrange::Tuple{T, T}
    latrange::Tuple{T, T}
end
function Spherical3DGrid_Reduced_EqSim(T::Type{<:AbstractFloat}, r0::Real,
            h::AbstractVector, N_lat::Integer; lonrange=(-pi, pi), latmax=pi/2)

    # input assertions
    @assert sort(h)[1] > 0 "Invalid heights"
    @assert lonrange[1] <= lonrange[2] "Invalid longitude range"
    @assert latmax >= 0 "Invalid `latmax`"

    # prepare coordiante calculation
    h, N_r = sort(h), length(h)
    dlat = dlon_max = T(latmax/N_lat)
    lat0 = T.(range(1/2/N_lat, 1 - 1/2/N_lat, length=N_lat) * latmax)

    # coordinate calculation
    r, lon, lat, N_lon = T[], T[], T[], Int64[]
    for p0 in lat0
        n_lon = ceil(Int64, (lonrange[2] - lonrange[1])*cos(p0)/dlon_max + eps(T))
        push!(lon, range(1/2/n_lon, 1 - 1/2/n_lon, length=n_lon) *
            (lonrange[2]-lonrange[1]) .+ lonrange[1]...)
        push!(lat, repeat([p0], n_lon)...)
        push!(N_lon, n_lon)
    end
    r = T.(repeat(r0 .+ vcat(0, h[1:end-1]), inner=sum(N_lon)))
    lon = repeat(lon, outer=N_r)
    lat = repeat(lat, outer=N_r)

    # area and volume calculation
    areas, volumes = T[], T[]
    for i in 1:N_r, j in 1:N_lat, _ in 1:N_lon[j]
        ri, dr = [r0 .+ vcat(0, h[1:end-1])...][i], vcat(h[1], diff(h))[i]
        dlon = (lonrange[2] - lonrange[1]) / N_lon[j]
        push!(areas, ri^2 * dlon * (sin(lat0[j]+dlat/2) - sin(lat0[j]-dlat/2)))
        push!(volumes, ((ri + dr)^3 - ri^3)/3 * dlon * (sin(lat0[j]+dlat/2) - sin(lat0[j]-dlat/2)))
    end

    return Spherical3DGrid_Reduced_EqSim(T(r0), T.(h), N_r, N_lon, N_lat,
                                         GlobalSphericalPosition.(r, lon, lat),
                                         T.(areas), T.(volumes), T.(lonrange), T.((0,latmax)))
end
function Spherical3DGrid_Reduced_EqSim(r0::Real, h::AbstractVector, N_lat::Integer; kwargs...)
    return Spherical3DGrid_Reduced_EqSim(typeof(r0), r0, h, N_lat; kwargs...)
end
function Spherical3DGrid_Reduced_EqSim(r0::Integer, h::AbstractVector, N_lat::Integer; kwargs...)
    return Spherical3DGrid_Reduced_EqSim(promote_type(typeof(r0), Float64), r0, h, N_lat; kwargs...)
end
function Spherical3DGrid_Reduced_EqSim(T::Type{<:AbstractFloat}, r::AbstractVector,
            N_lat::Integer; kwargs...)
    r0, h = r[1], accumulate(+,diff(r))
    Spherical3DGrid_Reduced_EqSim(T, r0, h, N_lat; kwargs...)
end
function Spherical3DGrid_Reduced_EqSim(r::AbstractVector, N_lat::Integer; kwargs...)
    return Spherical3DGrid_Reduced_EqSim(typeof(r[1]), r, N_lat; kwargs...)
end
function Spherical3DGrid_Reduced_EqSim(r::AbstractVector{<:Integer}, N_lat::Integer; kwargs...)
    return Spherical3DGrid_Reduced_EqSim(promote_type(typeof(r[1]), Float64), r, N_lat; kwargs...)
end


############################################################################################
#::. UTILITY FUNCTIONS
############################################################################################
function coord2idx(grid::Spherical3DGrid, r::T, lon::T, lat::T)::Int64 where {T<:AbstractFloat}
    PI     = T(pi) # to prevent numerical cutoff/rounding mistakes
    lon    = pclamp(lon, -PI, PI)
    lat    = clamp(lat, -PI/2, PI/2)

    lonrange, latrange = grid.lonrange, grid.latrange
    if !(grid.r0 <= r < grid.r0 + grid.h[end]); return 0; end
    if !(lonrange[1] <= lon <= lonrange[2]); return 0; end
    if !(latrange[1] <= lat <= latrange[2]); return 0; end

    idxr   = findfirst(grid.h .+ grid.r0 .> r)
    idxlon = max(1,ceil(Int64, (lon-lonrange[1])/(lonrange[2]-lonrange[1])*grid.N_lon)) # max(...) to prevent index < 1 at boundary
    idxlat = max(1,ceil(Int64, (lat-latrange[1])/(latrange[2]-latrange[1])*grid.N_lat))
    return (idxr-1) * grid.N_lat*grid.N_lon + (idxlon-1) * grid.N_lat + idxlat
end
function coord2idx(grid::Spherical3DGrid_EqSim, r::T, lon::T, lat::T)::Int64 where {T<:AbstractFloat}
    PI     = T(pi) # to prevent numerical cutoff/rounding mistakes
    lon    = pclamp(lon, -PI, PI)
    lat    = clamp(lat, -PI/2, PI/2)

    lonrange, latrange = grid.lonrange, grid.latrange
    if !(grid.r0 <= r < grid.r0 + grid.h[end]); return 0; end
    if !(lonrange[1] <= lon <= lonrange[2]); return 0; end
    if !(-latrange[2] <= lat <= latrange[2]); return 0; end

    idxr   = findfirst(grid.h .+ grid.r0 .> r)
    idxlon = max(1,ceil(Int64, (lon-lonrange[1])/(lonrange[2]-lonrange[1])*grid.N_lon))
    idxlat = max(1,ceil(Int64, abs(lat)/latrange[2]*grid.N_lat))
    return (idxr-1) * grid.N_lat*grid.N_lon + (idxlon-1) * grid.N_lat + idxlat
end
function coord2idx(grid::Spherical3DGrid_Reduced, r::T, lon::T, lat::T)::Int64 where {T<:AbstractFloat}
    PI     = T(pi) # to prevent numerical cutoff/rounding mistakes
    lon    = pclamp(lon, -PI, PI)
    lat    = clamp(lat, -PI/2, PI/2)

    lonrange, latrange = grid.lonrange, grid.latrange
    if !(grid.r0 <= r < grid.r0 + grid.h[end]); return 0; end
    if !(lonrange[1] <= lon <= lonrange[2]); return 0; end
    if !(latrange[1] <= lat <= latrange[2]); return 0; end

    idxr   = findfirst(grid.h .+ grid.r0 .> r)
    idxlat = max(1,ceil(Int64, (lat-latrange[1])/(latrange[2]-latrange[1])*grid.N_lat))
    idxlon = max(1,ceil(Int64, (lon-lonrange[1])/(lonrange[2]-lonrange[1])*grid.N_lon[idxlat]))
    if idxlat == 1; return (idxr-1) * sum(grid.N_lon) + idxlon; end
    return (idxr-1) * sum(grid.N_lon) + idxlon + accumulate(+, grid.N_lon[1:idxlat])[end-1]
end
function coord2idx(grid::Spherical3DGrid_Reduced_EqSim, r::T, lon::T, lat::T)::Int64 where {T<:AbstractFloat}
    PI     = T(pi) # to prevent numerical cutoff/rounding mistakes
    lon    = pclamp(lon, -PI, PI)
    lat    = clamp(lat, -PI/2, PI/2)

    lonrange, latrange = grid.lonrange, grid.latrange
    if !(grid.r0 <= r < grid.r0 + grid.h[end]); return 0; end
    if !(lonrange[1] <= lon <= lonrange[2]); return 0; end
    if !(-latrange[2] <= lat <= latrange[2]); return 0; end

    idxr   = findfirst(grid.h .+ grid.r0 .> r)
    idxlat = max(1,ceil(Int64, abs(lat)/latrange[2]*grid.N_lat))
    idxlon = max(1,ceil(Int64, (lon-lonrange[1])/(lonrange[2]-lonrange[1])*grid.N_lon[idxlat]))
    if idxlat == 1; return (idxr-1) * sum(grid.N_lon) + idxlon; end
    return (idxr-1) * sum(grid.N_lon) + idxlon + accumulate(+, grid.N_lon[1:idxlat])[end-1]

    # if !(grid.r0 <= r <= grid.r0 + grid.h[end]); return 0; end
    # PI       = T(pi) # to prevent numerical cutoff/rounding mistakes
    # lon    = pclamp(lon, -PI, PI - eps(T)) + eps(T)
    # lat      = clamp(lat, -PI/2, PI/2 - eps(T)) + eps(T)
    # idxr     = findfirst(grid.h .+ grid.r0 .> r)
    # idxlat   = ceil(Int64, abs(lat)*grid.N_lat*2/pi)
    # idxlon = ceil(Int64, (lon+pi)/2/pi*grid.N_lon[idxlat])
    # if idxlat == 1; return (idxr-1) * sum(grid.N_lon) + idxlon; end
    # return (idxr-1) * sum(grid.N_lon) + idxlon + accumulate(+, grid.N_lon[1:idxlat])[end-1]
end
function coord2idx(grid::AbstractSpherical3DGrid, r::Real, lon::Real, lat::Real)::Int64
    return coord2idx(grid, promote(r, lon, lat)...)
end
function coord2idx(grid::AbstractSpherical3DGrid, r::Integer, lon::Integer, lat::Integer)::Int64
    return coord2idx(grid, float(r), float(lon), float(lat))
end


surfacecoords(grid::AbstractSpherical3DGrid) = grid.coords[1:grid.N_lon*grid.N_lat]
surfacecoords(grid::Spherical3DGrid_Reduced) = grid.coords[1:sum(grid.N_lon)]
surfacecoords(grid::Spherical3DGrid_Reduced_EqSim) = grid.coords[1:sum(grid.N_lon)]


############################################################################################
#::. EXTENSIONS
############################################################################################
Base.length(grid::AbstractSpherical3DGrid) = length(grid.coords)


Base.size(grid::AbstractSpherical3DGrid) = (grid.N_r, grid.N_lon, grid.N_lat)


Base.show(io::IO, ::MIME"text/plain", grid::AbstractSpherical3DGrid) =
    print(io, "$(typeof(grid)):\n"*
            " r0:      $(grid.r0)\n"*
            " h:       $([h for h in grid.h])\n"*
            " N_r:     $(grid.N_r)\n"*
            " N_lon: $(grid.N_lon)\n"*
            " N_lat:   $(grid.N_lat)\n"*
            " coords:  $(length(grid.coords))-element $(typeof(coords(grid)))\n"*
            " areas:   $(length(grid.areas))-element $(typeof(grid.areas))\n"*
            "            min: $(min(grid.areas...))\n"*
            "            max: $(max(grid.areas...))\n"*
            " volumes: $(length(grid.volumes))-element $(typeof(grid.volumes))\n"*
            "            min: $(min(grid.volumes...))\n"*
            "            max: $(max(grid.volumes...))\n")

Base.show(io::IO, ::MIME"text/plain", grid::Spherical3DGrid_Reduced) =
    print(io, "$(typeof(grid)):\n"*
            " r0:      $(grid.r0)\n"*
            " h:       $([h for h in grid.h])\n"*
            " N_lon: $(length(grid.N_lon))-element Vector{Int64}\n"*
            "            @ quator: $(grid.N_lon[1])\n"*
            "            @ poles:  $(grid.N_lon[end])\n"*
            " N_lat:   $(grid.N_lat)\n"*
            " coords:  $(length(grid.coords))-element $(typeof(coords(grid)))\n"*
            " areas:   $(length(grid.areas))-element $(typeof(grid.areas))\n"*
            "            min: $(min(grid.areas...))\n"*
            "            max: $(max(grid.areas...))\n"*
            " volumes: $(length(grid.volumes))-element $(typeof(grid.volumes))\n"*
            "            min: $(min(grid.volumes...))\n"*
            "            max: $(max(grid.volumes...))\n")


############################################################################################
#::. EXPORTS
############################################################################################
export
    Spherical3DGrid,
    Spherical3DGrid_EqSim,
    Spherical3DGrid_Reduced,
    Spherical3DGrid_Reduced_EqSim
