############################################################################################
#::. TYPES
############################################################################################
abstract type AbstractGrid end

abstract type AbstractCartesianGrid <: AbstractGrid end
abstract type AbstractSphericalGrid <: AbstractGrid end


############################################################################################
#::. FUNCTIONS
############################################################################################
(grid::AbstractGrid)(idx::Integer) = grid.coords[idx]

"""
    [1] areas([T::Type,] grid::AbstractGrid)

For 2D grids, returns the surface area of each grid element. For 3D grids, returns the
base area of each grid element.
"""
areas(T::Type, grid::AbstractGrid) = T.(areas(grid))
areas(grid::AbstractGrid) = grid.areas


"""
    [1] coords([T::Type,] grid::AbstractGrid)

Returns the coordinates of each grid element.
"""
coords(T::Type, grid::AbstractGrid) = T.(coords(grid))
coords(grid::AbstractGrid) = grid.coords


"""
    [1] mapgrid(x::AbstractVector, grid_from::AbstractGrid, grid_to::AbstractGrid)

Maps the values of `x` at the coordinates defined on `grid_from` to the coordinates defined
on `grid_to`. The output is a vector of the same length as `grid_to` with the values of `x`.
The function **does not** interpolate the values of `x` between the grid elements.
"""
function mapgrid(x::AbstractVector, grid_from::AbstractGrid, grid_to::AbstractGrid)
    coords_to = coords(grid_to)
    idx_from = coord2idx(grid_from, coords_to)
    return x[idx_from]
end


"""
    [1] coord2idx(grid::AbstractGrid, coord::AbstractPosition)
    [2] coord2idx(grid::AbstractGrid, coords::Vector{T}) where {T <: AbstractPosition}
    [3] coord2idx(grid::AbstractSphericalGrid, r::AbstractVector, theta::AbstractVector,
                  phi::AbstractVector)

Calculates the index of the grid element containing the given coordinates.
"""
function coord2idx(grid::AbstractGrid, coord::AbstractPosition)::Int64
    return coord2idx(grid, _get(coord)...)
end
function coord2idx(grid::AbstractGrid, coords::Vector{T}) where {T <: AbstractPosition}
    return [coord2idx(grid, coord) for coord in coords]
end
function coord2idx(grid::AbstractSphericalGrid, coord::AbstractGlobalPosition)::Int64
    #TODO: Do not use AbstractGloablPosition union type, but why?
    return coord2idx(grid, _get(GlobalSphericalPosition(coord))...)
end
function coord2idx(grid::AbstractSphericalGrid, r::AbstractVector, theta::AbstractVector,
                   phi::AbstractVector)
    return [coord2idx(grid, r[i], theta[i], phi[i]) for i in eachindex(r)]
end


"""
    [1] surfacecoords([T::Type,] grid::AbstractGrid)

Returns only the coordinates of the surface (i.e. the base) of the discretized geometry.
"""
surfacecoords(T::Type, grid::AbstractGrid) = T.(surfacecoords(grid))


"""
    [1] volumes([T::Type,] grid::AbstractGrid)

For 2D grids, returns a vector of `zeros(T)` for each grid element. For 3D grids, returns the
volume of each grid element.
"""
volumes(T::Type, grid::AbstractGrid) = T.(volumes(grid))
volumes(grid::AbstractGrid) = grid.volumes


############################################################################################
#::. EXPORTS
############################################################################################
export areas, coords, coord2idx, mapgrid, surfacecoords, volumes
