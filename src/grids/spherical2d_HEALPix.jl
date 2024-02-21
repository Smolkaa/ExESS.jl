#::. types
abstract type AbstractSpherical2DGrid_HEALPix <: AbstractSphericalGrid end


#::. FUNCTIONS
"""
    [1] struct Spherical2DGrid_HEALPix{T<:AbstractFloat} <: AbstractSpherical2DGrid_HEALPix
    [2] Spherical2DGrid_HEALPix([T::Type,] r::Real, k::Int64)

Global grid of surface coordinates (2D) of type `Spherical2DGrid_HEALPix{T}` over
a sphere of radius `r`. The spiral is created based on the HEALPix (Hierarchical Equal Area
isoLatitude Pixelisation) sequence and the Lambert cylindrical equal area projection. The 
surface area `area` is a scalar value and assumes a prefectly equal distribution of points. 

| Field    | Type                                 | Description                    |
|:-------- |:------------------------------------ |:------------------------------ |
| `r`      | `T`                                  | radius of global body (sphere) |
| `k`      | `Int64`                              | number of rings                |
| `coords` | `Vector{GlobalSphericalPosition{T}}` | coordinates                    |
| `area`   | `T`                                  | surface area                   |
| `tree`   | `KDTree`                             | nearest neighbors tree         |

To ensure expected behavior, the grid object should generally be created with the outer
constructor [2].


## Number of Grid_HEALPix Elements

The number of elements in the `Spherical2DGrid_HEALPix` is not explicitly specified, but is
instead determined by the number of rings `k`: `N = 12*k^2`.
"""
struct Spherical2DGrid_HEALPix{T<:AbstractFloat} <: AbstractSpherical2DGrid_HEALPix
    r::T
    k::Int64
    coords::Vector{GlobalSphericalPosition{T}}
    area::T
    tree::KDTree
end
function Spherical2DGrid_HEALPix(T::Type, r::Real, k::Int64)
    theta, phi = T[], T[]
    
    # polar region (k rings)
    for i in 1:k, j in 1:4*i
        push!(theta, (pi/(2*i)) * (j-0.5) - pi)
        push!(phi, acos(1 - i^2/(3*k^2)) - pi/2)
    end

    # equator region (2k rings)
    for i in k+1:2*k, j in 1:4*k
        s = mod((i-k+1), 2)
        push!(theta, (pi/(2*k)*(j-s/2)) - pi)
        push!(phi, acos(4/3 - (2*i)/(3*k)) - pi/2)
    end

    # full grid
    theta = vcat(theta, reverse(theta[1:end-4*k]))
    phi = vcat(phi, reverse(-phi[1:end-4*k]))

    # check theta boundaries
    for i in eachindex(theta)
        if abs(theta[i]) > pi; theta[i] = sign(theta[i]) * (pi - eps(T)); end
    end

    # area calculation (equal area)
    area = 4*pi*r^2 / (12*k^2)

    # setup KDTree
    x, y, z =  cos.(theta) .* cos.(phi), sin.(theta) .* cos.(phi), sin.(phi)
    tree = KDTree([x'; y'; z'])

    return Spherical2DGrid_HEALPix(T(r), k, GlobalSphericalPosition.(T(r), T.(theta), T.(phi)), T(area), tree)
end
function Spherical2DGrid_HEALPix(r::Real, k::Int64)
    if typeof(r) <: Integer; r = Float64(r); end
    return Spherical2DGrid_HEALPix(typeof(r), r, k)
end


#::. utility functions
areas(grid::Spherical2DGrid_HEALPix) = length(grid) .* grid.area

function coord2idx(grid::Spherical2DGrid_HEALPix, theta::Real, phi::Real)
    x, y, z = _get(GlobalCartesianPosition(GlobalSphericalPosition(1.0, theta, phi)))
    idxs, _ = knn(grid.tree, [x, y, z], 1, true)
    return idxs[1]
end
function coord2idx(grid::AbstractSpherical2DGrid_HEALPix, r::Real, theta::Real, phi::Real) 
    return coord2idx(grid, theta, phi)
end


#::. EXTENSIONS
Base.length(grid::AbstractSpherical2DGrid_HEALPix) = length(grid.coords)
Base.size(grid::AbstractSpherical2DGrid_HEALPix) = (length(grid), )

Base.show(io::IO, ::MIME"text/plain", grid::AbstractSpherical2DGrid_HEALPix) = 
    print(io, "$(typeof(grid)):\n"*
            " r:       $(grid.r)\n"*
            " N:       $(length(grid))\t(k=$(grid.k))\n"*
            " coords:  $(length(grid))-element $(typeof(coords(grid)))\n"*
            " area:    $(grid.area)")


#::. EXPORTS
export Spherical2DGrid_HEALPix