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
    [1] coord2idx(grid::AbstractGrid, coord::AbstractPosition)
    [2] coord2idx(grid::AbstractSphericalGrid, r::AbstractVector, theta::AbstractVector, phi::AbstractVector)
    [2] coord2idx(grid::AbstractSphericalGrid, coords::Vector{AbstractPosition})

Calculates the index of the grid element containing the given coordinates.
"""
coord2idx(grid::AbstractGrid, coord::AbstractPosition) = coord2idx(grid, tup(coord)...)
function coord2idx(grid::AbstractSphericalGrid, coord::AbstractGlobalPosition) #TODO: Do not use AbstractGloablPosition union type
    return coord2idx(grid, tup(GlobalSphericalPosition(coord))...)
end
function coord2idx(grid::AbstractSphericalGrid, coords::Vector{AbstractPosition})
    return [coord2idx(grid, coord) for coord in coords]
end
function coord2idx(grid::AbstractSphericalGrid, r::AbstractVector, theta::AbstractVector, phi::AbstractVector)
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
export areas, coords, coord2idx, surfacecoords, volumes
