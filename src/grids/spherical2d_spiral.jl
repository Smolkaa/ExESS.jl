############################################################################################
#::. TYPES
############################################################################################
abstract type AbstractSpherical2DGrid_Spiral <: AbstractSphericalGrid end


############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    [1] struct Spherical2DGrid_Spiral{T<:AbstractFloat} <: AbstractSpherical2DGrid_Spiral end
    [2] Spherical2DGrid_Spiral([T::Type,] r::Real, N::Integer)

Global spiral grid of surface coordinates (2D) of type `GlobalSphericalPosition{T}` over
a sphere of radius `r`. The spiral is created based on the Fibonnaci sequence and the
Lambert cylindrical equal area projection. The surface area `area` is a scalar value and
assumes a prefectly equal distribution of points.

| Field    | Type; with `T<:AbstractFloat`        | Description                    |
|:-------- |:------------------------------------ |:------------------------------ |
| `r`      | `T`                                  | radius of global body (sphere) |
| `N`      | `Int64`                              | # elements                     |
| `coords` | `Vector{GlobalSphericalPosition{T}}` | coordinates                    |
| `area`   | `T`                                  | surface area                   |
| `tree`   | `KDTree`                             | nearest neighbors tree         |

To ensure expected behavior, the grid object should generally be created with the outer
constructor [2].
"""
struct Spherical2DGrid_Spiral{T<:AbstractFloat} <: AbstractSpherical2DGrid_Spiral
    r::T
    N::Integer
    coords::Vector{GlobalSphericalPosition{T}}
    area::T
    tree::KDTree
end
function Spherical2DGrid_Spiral(T::Type{<:AbstractFloat}, r::Real, N::Integer)
    GOLDEN_RATIO = (1 + sqrt(5)) / 2

    # calculate coordinates
    theta = 2pi*[mod(i/GOLDEN_RATIO, 1) for i in 0:N-1] .- pi
    phi = acos.([2*i/(N-1) for i in 0:N-1] .- 1) .- pi/2 # Lambert cylindrical equal area projection
    area = 4*pi*r^2 / N

    # setup KDTree
    x, y, z =  cos.(theta) .* cos.(phi), sin.(theta) .* cos.(phi), sin.(phi)
    tree = KDTree([x'; y'; z'])

    return Spherical2DGrid_Spiral(T(r), N,
                                  GlobalSphericalPosition.(T(r), T.(theta), T.(phi)),
                                  T(area), tree)
end
function Spherical2DGrid_Spiral(r::Real, N::Integer)
    return Spherical2DGrid_Spiral(typeof(r), r, N)
end
function Spherical2DGrid_Spiral(r::Integer, N::Integer)
    return Spherical2DGrid_Spiral(promote_type(typeof(r), Float64), r, N)
end


############################################################################################
#::. UTILITY FUNCTIONS
############################################################################################
areas(grid::AbstractSpherical2DGrid_Spiral) = repeat([grid.area], grid.N)


function coord2idx(grid::Spherical2DGrid_Spiral, theta::Real, phi::Real)::Int64
    x, y, z = _get(GlobalCartesianPosition(GlobalSphericalPosition(1.0, theta, phi)))
    idxs, _ = knn(grid.tree, [x, y, z], 1, true)
    return idxs[1]
end
function coord2idx(grid::AbstractSpherical2DGrid_Spiral, r::Real, theta::Real, phi::Real)::Int64
    return coord2idx(grid, theta, phi)
end


surfacecoords(grid::AbstractSpherical2DGrid_Spiral) = coords(grid)


volumes(grid::AbstractSpherical2DGrid_Spiral) = zeros(typeof(grid.r), length(grid.coords))


############################################################################################
#::. EXTENSIONS
############################################################################################
Base.length(grid::AbstractSpherical2DGrid_Spiral) = grid.N


Base.size(grid::AbstractSpherical2DGrid_Spiral) = (grid.N, )


Base.show(io::IO, ::MIME"text/plain", grid::AbstractSpherical2DGrid_Spiral) =
    print(io, "$(typeof(grid)):\n"*
            " r:       $(grid.r)\n"*
            " N:       $(grid.N)\n"*
            " coords:  $(length(grid.coords))-element $(typeof(coords(grid)))\n"*
            " area:    $(grid.area)")


############################################################################################
#::. EXPORTS
############################################################################################
export Spherical2DGrid_Spiral
