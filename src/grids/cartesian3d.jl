#::. types
abstract type AbstractCartesian3DGrid <: AbstractCartesianGrid end


#::. FUNCTIONS
"""
    [1] struct Cartesian3DGrid{T<:AbstractFloat} <: AbstractCartesian3DGrid
    [2] Cartesian3DGrid([T::Type,] xrange::Tuple, yrange::Tuple, zrange::Tuple, N_x::Int64, N_y::Int64, N_z::Int64)
    [3] Cartesian3DGrid([T::Type,] xrange::Tuple, yrange::Tuple, z2::Real, N_x::Int64, N_y::Int64, N_z::Int64)

Local structured volume grid (3D) of type `Cartesian3DGrid{T}` over a rectangular
domain. The domain is defined by tuples of the form `xrange=(x1, x2)`, `yrange=(y1, y2)`, 
and `zrange=(z1, z2)`. The number of grid points in each direction is given by `N_x`, 
`N_y`, and `N_z`. The grid points are centered in their respective grid element.

| Field    | Type                                | Description                           |
|:-------- |:----------------------------------- |:------------------------------------- |
| `xrange` | `Tuple{T, T}`                       | x domain limits                       |
| `yrange` | `Tuple{T, T}`                       | y domain limits                       |
| `zrange` | `Tuple{T, T}`                       | z domain limits                       |
| `N_x`    | `Int64`                             | number of grid points in x-direction  |
| `N_y`    | `Int64`                             | number of grid points in y-direction  |
| `N_z`    | `Int64`                             | number of grid points in z-direction  |
| `coords` | `Vector{LocalCartesianPosition{T}}` | coordinates                           |
| `area`   | `T`                                 | grid-element surface areas (xy-plane) |
| `volume` | `T`                                 | grid-element volume                   |

Note that the fields `area` and `volume` are scalar values, aplying to all grid elements. For a 
correct, vectorized output, please use the grid-utility-methods `areas(grid)` and `volumes(grid)`.
 
To ensure expected behavior, the grid object should generally be created with the outer
constructors [2] to [3].
"""
struct Cartesian3DGrid{T<:AbstractFloat} <: AbstractCartesian3DGrid
    xrange::Tuple{T, T}
    yrange::Tuple{T, T}
    zrange::Tuple{T, T}
    N_x::Int64
    N_y::Int64
    N_z::Int64
    coords::Vector{LocalCartesianPosition{T}}
    area::T
    volume::T
end
function Cartesian3DGrid(T::Type, xrange::Tuple, yrange::Tuple, zrange::Tuple, 
                               N_x::Integer, N_y::Integer, N_z::Integer)
    # coordinates
    x1, x2 = sort([xrange[1], xrange[2]]) |> Vector{T}
    y1, y2 = sort([yrange[1], yrange[2]]) |> Vector{T}
    z1, z2 = sort([zrange[1], zrange[2]]) |> Vector{T}
    dx, dy, dz = x2 - x1, y2 - y1, z2 - z1
    x = range(x1 + dx/N_x/2, x2 - dx/N_x/2, N_x)
    y = range(y1 + dy/N_y/2, y2 - dy/N_y/2, N_y)
    z = range(z1 + dz/N_z/2, z2 - dz/N_z/2, N_z)
    # z = range(zrange[1] + dz/N_z, zrange[2], N_z)

    # coordinates
    X = repeat(x, outer=N_y*N_z)
    Y = repeat(y, inner=N_x, outer=N_z)
    Z = repeat(z, inner=N_x*N_y)
    coords = LocalCartesianPosition.(X, Y, Z)

    # area and volume
    A = dx * dy / N_x / N_y
    V = A * dz / N_z
    
    return Cartesian3DGrid((x1,x2), (y1,y2), (z1,z2), N_x, N_y, N_z, coords, T(A), T(V))
end
function Cartesian3DGrid(xrange::Tuple, yrange::Tuple, zrange::Tuple, 
                               N_x::Integer, N_y::Integer, N_z::Integer)
    x1, x2, y1, y2, z1, z2 = promote(xrange..., yrange..., zrange...)
    if typeof(x1) <: Integer; T = Float64; else; T = typeof(x1); end
    return Cartesian3DGrid(T, (x1,x2), (y1,y2), (z1,z2), N_x, N_y, N_z)
end
function Cartesian3DGrid(T::Type, xrange::Tuple, yrange::Tuple, z2::Real,
                               N_x::Integer, N_y::Integer, N_z::Integer)
    return Cartesian3DGrid(T, (x1,x2), (y1,y2), (0,z2), N_x, N_y, N_z)
end
function Cartesian3DGrid(xrange::Tuple, yrange::Tuple, z2::Real, 
                               N_x::Integer, N_y::Integer, N_z::Integer)
    x1, x2, y1, y2, z2 = promote(xrange[1], xrange[2], yrange[1], yrange[2], z2)
    if typeof(x1) <: Integer; T = Float64; else; T = typeof(x1); end
    return Cartesian3DGrid(T, (x1,x2), (y1,y2), (0,z2), N_x, N_y, N_z)
end


"""
    [1] Cartesian3DGrid_Exponential{T<:AbstractFloat} <: AbstractCartesian3DGrid
    [2] Cartesian3DGrid_Exponential([T::Type,] xrange::Tuple, yrange::Tuple, zrange::Tuple, N_x::Int64, N_y::Int64, N_z::Int64; c=1.0)
    [3] Cartesian3DGrid_Exponential([T::Type,] xrange::Tuple, yrange::Tuple, z2::Real, N_x::Int64, N_y::Int64, N_z::Int64; kwargs...)

Local structured volume grid (3D) of type `Cartesian3DGrid_Exponential{T}` over a rectangular
domain. The domain is defined by tuples of the form `xrange=(x1, x2)`, `yrange=(y1, y2)`, 
and `zrange=(z1, z2)`. The number of grid points in each direction is given by `N_x`, 
`N_y`, and `N_z`. In contrast to the `Cartesian3DGrid`, the grid points are distributed
exponentially in the vertical direction, i.e. with exponentially increasing spacing with
higher `z`-coordinates. The grid points are centered in their respective grid element.

The additional keyword argument `c` controls the shape of the exponential distribution of
heights. Defaults to `c=1`, increasing values increase the "curviness", decreasing values
lead to a more linear distribution.


| Field     | Type                                | Description                            |
|:--------- |:----------------------------------- |:---------------------------------------|
| `xrange`  | `Tuple{T, T}`                       | x domain limits                        |
| `yrange`  | `Tuple{T, T}`                       | y domain limits                        |
| `zrange`  | `Tuple{T, T}`                       | z domain limits                        |
| `N_x`     | `Int64`                             | number of grid points in x-direction   |
| `N_y`     | `Int64`                             | number of grid points in y-direction   |
| `N_z`     | `Int64`                             | number of grid points in z-direction   |
| `coords`  | `Vector{LocalCartesianPosition{T}}` | coordinates                            |
| `area`    | `T`                                 | grid-element surface areas (xy-plane)  |
| `volumes` | `Vector{T}`                         | grid-element volumes                   |
| `c`       | `Float64`                           | shape factor of exp. distribution      |


Note that the fields `area` and `volume` are auxiliary values. For a correct, vectorized
output, please use the grid-utility-methods `areas(grid)` and `volumes(grid)`.
 
To ensure expected behavior, the grid object should generally be created with the outer
constructors [2] to [3].
"""
struct Cartesian3DGrid_Exponential{T<:AbstractFloat} <: AbstractCartesian3DGrid
    xrange::Tuple{T, T}
    yrange::Tuple{T, T}
    zrange::Tuple{T, T}
    N_x::Int64
    N_y::Int64
    N_z::Int64
    coords::Vector{LocalCartesianPosition{T}}
    area::T
    volumes::Vector{T}
    c::Float64
end
function Cartesian3DGrid_Exponential(T::Type, xrange::Tuple, yrange::Tuple, zrange::Tuple, N_x::Int64, N_y::Int64, N_z::Int64; c=1.0)
    # coordinates
    x1, x2 = sort([xrange[1], xrange[2]]) |> Vector{T}
    y1, y2 = sort([yrange[1], yrange[2]]) |> Vector{T}
    z1, z2 = sort([zrange[1], zrange[2]]) |> Vector{T}
    dx, dy, dz = x2 - x1, y2 - y1, z2 - z1
    x = range(x1 + dx/N_x/2, x2 - dx/N_x/2, N_x) |> Vector{T}
    y = range(y1 + dy/N_y/2, y2 - dy/N_y/2, N_y) |> Vector{T}
    z = z1 .+ [dz * (exp(c/N_z)^(i-0.5) - 1) / (exp(c) - 1) for i in 1:N_z] |> Vector{T}

    # coordinates
    X = repeat(x, outer=N_y*N_z)
    Y = repeat(y, inner=N_x, outer=N_z)
    Z = repeat(z, inner=N_x*N_y)
    coords = LocalCartesianPosition.(X, Y, Z)

    # area and volumes
    A = dx * dy / N_x / N_y
    V = T[]
    for k in 1:N_z
        dz = k == 1   ? (z[2]   - z[1])     / 2 + (z[1] - z1) : 
             k == N_z ? (z[N_z] - z[N_z-1]) / 2 + (z2   - z[N_z]) : 
                        (z[k+1] - z[k-1])   / 2
        push!(V, A * dz)
    end
    
    return Cartesian3DGrid_Exponential((x1,x2), (y1,y2), (z1,z2), N_x, N_y, N_z, coords, T(A), V, Float64(c))
end
function Cartesian3DGrid_Exponential(xrange::Tuple, yrange::Tuple, zrange::Tuple, N_x::Int64, N_y::Int64, N_z::Int64; kwargs...)
    x1, x2, y1, y2, z1, z2 = promote(xrange[1], xrange[2], yrange[1], yrange[2], zrange[1], zrange[2])
    if typeof(x1) <: Integer; T = Float64; else; T = typeof(x1); end
    return Cartesian3DGrid_Exponential(T, (x1,x2), (y1,y2), (z1,z2), N_x, N_y, N_z; kwargs...)
end
function Cartesian3DGrid_Exponential(T::Type, xrange::Tuple, yrange::Tuple, z2::Real, N_x::Int64, N_y::Int64, N_z::Int64; kwargs...)
    return Cartesian3DGrid_Exponential(T, xrange, yrange, (0,z2), N_x, N_y, N_z; kwargs...)
end
function Cartesian3DGrid_Exponential(xrange::Tuple, yrange::Tuple, z2::Real, N_x::Int64, N_y::Int64, N_z::Int64; kwargs...)
    x1, x2, y1, y2, z2 = promote(xrange[1], xrange[2], yrange[1], yrange[2], z2)
    if typeof(x1) <: Integer; T = Float64; else; T = typeof(x1); end
    return Cartesian3DGrid_Exponential(T, xrange, yrange, (0,z2), N_x, N_y, N_z; kwargs...)
end




#::. utility functions
areas(grid::AbstractCartesian3DGrid) = grid.area .* ones(length(grid))

function coord2idx(grid::Cartesian3DGrid, x::Real, y::Real, z::Real)
    x1, x2, y1, y2, z1, z2 = grid.xrange[1], grid.xrange[2], grid.yrange[1], grid.yrange[2], grid.zrange[1], grid.zrange[2]
    if !(x1 <= x <= x2 && y1 <= y <= y2 && z1 <= z <= z2); return 0; end

    dx, dy, dz = (x2 - x1)/grid.N_x, (y2 - y1)/grid.N_y, (z2 - z1)/grid.N_z

    i = ceil(Int64, (x - x1)/dx)
    j = ceil(Int64, (y - y1)/dy)
    k = ceil(Int64, (z - z1)/dz)
    return (k-1)*grid.N_x*grid.N_y + (j-1)*grid.N_x + i
end
function coord2idx(grid::Cartesian3DGrid_Exponential, x::Real, y::Real, z::Real)
    x1, x2, y1, y2, z1, z2 = grid.xrange[1], grid.xrange[2], grid.yrange[1], grid.yrange[2], grid.zrange[1], grid.zrange[2]
    if !(x1 <= x <= x2 && y1 <= y <= y2 && z1 <= z <= z2); return 0; end

    dx, dy = (x2 - x1)/grid.N_x, (y2 - y1)/grid.N_y
    Z = z1 .+ [z2 - z1 * (exp(grid.c/grid.N_z)^(i-0.5) - 1) / (exp(grid.c) - 1) for i in 1:grid.N_z]

    i = ceil(Int64, (x - x1)/dx)
    j = ceil(Int64, (y - y1)/dy)
    k = findmin(abs.(Z .- z))[2]
    return (k-1)*grid.N_x*grid.N_y + (j-1)*grid.N_x + i
end

surfacecoords(grid::AbstractCartesian3DGrid) = grid.coords[1:grid.N_x*grid.N_y]

volumes(grid::AbstractCartesian3DGrid) = grid.volume .* ones(length(grid))
volumes(grid::Cartesian3DGrid_Exponential) = repeat(grid.volumes, inner=grid.N_x*grid.N_y)


#::. EXTENSIONS
Base.length(grid::AbstractCartesian3DGrid) = grid.N_x * grid.N_y * grid.N_z

Base.size(grid::AbstractCartesian3DGrid) = (grid.N_x, grid.N_y, grid.N_z)

Base.show(io::IO, ::MIME"text/plain", grid::Cartesian3DGrid) = 
    print(io, "$(typeof(grid)):\n"*
            " xrange:  $(grid.xrange)\n"*
            " yrange:  $(grid.yrange)\n"*
            " zrange:  $(grid.zrange)\n"*
            " N_x:     $(grid.N_x)\n"*
            " N_y:     $(grid.N_y)\n"*
            " N_z:     $(grid.N_z)\n"*
            " coords:  $(length(grid.coords))-element $(typeof(coords(grid)))\n"*
            " area:    $(grid.area)\n"*
            " volume:  $(grid.volume)\n")

Base.show(io::IO, ::MIME"text/plain", grid::Cartesian3DGrid_Exponential) = 
    print(io, "$(typeof(grid)):\n"*
            " xrange:  $(grid.xrange)\n"*
            " yrange:  $(grid.yrange)\n"*
            " zrange:  $(grid.zrange)\n"*
            " N_x:     $(grid.N_x)\n"*
            " N_y:     $(grid.N_y)\n"*
            " N_z:     $(grid.N_z)\n"*
            " coords:  $(length(grid.coords))-element $(typeof(coords(grid)))\n"*
            " area:    $(grid.area)\n"*
            " volumes: $(length(grid.volumes))-element $(typeof(grid.volumes))\n"*
            "            min: $(min(grid.volumes...))\n"*
            "            max: $(max(grid.volumes...))\n"*
            " c:       $(grid.c)\n")


#::. EXPORTS
export 
    Cartesian3DGrid, 
    Cartesian3DGrid_Exponential