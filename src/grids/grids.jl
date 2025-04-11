############################################################################################
#::. TYPES
############################################################################################
abstract type AbstractGrid end

abstract type AbstractCartesianGrid <: AbstractGrid end
abstract type AbstractSphericalGrid <: AbstractGrid end


############################################################################################
#::. FUNCTIONS
############################################################################################
(grid::AbstractGrid)(idx::Integer) = grid.coords[idx] # get grid element by index


"""
    areas([T], grid)

Returns the surface area of each grid element in `grid`. If `grid` is three-dimensional,
the function returns the base area of each grid element. The optional type `T` can be used
to control the output type of the areas.
"""
areas(grid::AbstractGrid) = grid.areas
areas(T::Type{<:AbstractFloat}, grid::AbstractGrid) = T.(areas(grid))


"""
    coords([T], grid)

Returns the coordinates of each grid element in `grid`. The optional type `T` can be used to
control the output type of the coordinates.
"""
coords(grid::AbstractGrid) = grid.coords
coords(T::Type{<:AbstractFloat}, grid::AbstractGrid) = T.(coords(grid))


"""
    mapgrid(x, grid_from, grid_to)

Maps the values of `x` at the coordinates defined on `grid_from` to the coordinates defined
on `grid_to`. The output is a vector of the same length as `grid_to` with the values of `x`.
The function **does not** interpolate the values of `x` between the grid elements.

# Arguments
- `x::AbstractVector`: Values to map.
- `grid_from::AbstractGrid`: Grid of the input values.
- `grid_to::AbstractGrid`: Grid of the output values.

# Notes
- There may be unexpected behaviour if the grids are of different dimensions, i.e. a 2D grid
  mapped to a 3D grid, or vice versa. The function does not assume any mathematical or
  physical meaning to the mapping or reduction of an additional dimension.
"""
function mapgrid(x::AbstractVector, grid_from::AbstractGrid, grid_to::AbstractGrid)
    coords_to = coords(grid_to)
    idx_from = coord2idx(grid_from, coords_to)
    return x[idx_from]
end


"""
    coord2idx(grid, coord)
    coord2idx(grid, [r], lon, lat)

Calculates the index of the grid element in `grid` containing the given coordinates defined
in `coord` or `r`, `lon`, and `lat`.

# Arguments
- `grid::AbstractGrid`: Grid to calculate the index.
- `coord::AbstractPosition` or `coord::Tuple`: Coordinates to calculate the index.
- `r::Real`: Radial coordinate. Can be omitted for 2D grids.
- `lon::Real`: Longitude coordinate.
- `lat::Real`: Latitude coordinate.

The coordinate `coord` or the individual coordinates `r`, `lon`, and `lat` can also be
vectors of the same length.
"""
coord2idx(grid::AbstractGrid, coord::Tuple)::Int64 = coord2idx(grid, coord...)
coord2idx(grid::AbstractGrid, coord::AbstractPosition)::Int64 = coord2idx(grid, Tuple(coord))
function coord2idx(grid::AbstractGrid, coords::Vector{T}) where {T <: AbstractPosition}
    return [coord2idx(grid, coord) for coord in coords]
end
function coord2idx(grid::AbstractSphericalGrid, coord::AbstractGlobalPosition)::Int64
    #TODO: Do not use AbstractGloablPosition union type, but why?
    return coord2idx(grid, Tuple(GlobalSphericalPosition(coord))...)
end
function coord2idx(grid::AbstractSphericalGrid, r::AbstractVector, theta::AbstractVector,
                   phi::AbstractVector)
    return [coord2idx(grid, r[i], theta[i], phi[i]) for i in eachindex(r)]
end


"""
    surfacecoords([T], grid)

Returns the surface coordinates of each grid element in `grid`. The optional type `T` can be
used to control the output type of the coordinates.
"""
surfacecoords(T::Type{<:AbstractFloat}, grid::AbstractGrid) = T.(surfacecoords(grid))


"""
    volumes([T], grid)

Returns the volume of each grid element in `grid`. The optional type `T` can be used to
control the output type of the volumes.

# Notes
- returns `zeros(T)` for 2D grids.
"""
volumes(grid::AbstractGrid) = grid.volumes
volumes(T::Type{<:AbstractFloat}, grid::AbstractGrid) = T.(volumes(grid))


############################################################################################
#::. EXPORTS
############################################################################################
export areas, coords, coord2idx, mapgrid, surfacecoords, volumes
