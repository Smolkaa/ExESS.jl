############################################################################################
#::. TYPES
############################################################################################
abstract type AbstractSpherical2DGrid <: AbstractSphericalGrid end


############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    [1] struct Spherical2DGrid{T<:AbstractFloat} <: AbstractSpherical2DGrid end
    [2] Spherical2DGrid([T::Type,] r::Real, N_theta::Integer, N_phi::Integer)

Global structured grid of surface coordinates (2D) of type `GlobalSphericalPosition{T}` over
a sphere of radius `r`.

| Field     | Type; with `T<:AbstractFloat`        | Description                       |
|:--------- |:------------------------------------ |:--------------------------------- |
| `r`       | `T`                                  | radius of global body (sphere)    |
| `N_theta` | `Int64`                              | # elements in azimuth direction   |
| `N_phi`   | `Int64`                              | # elements in elevation direction |
| `coords`  | `Vector{GlobalSphericalPosition{T}}` | coordinates                       |
| `areas`   | `Vector{T}`                          | surface area                      |

To ensure expected behavior, the grid object should generally be created with the outer
constructor [2].
"""
struct Spherical2DGrid{T<:AbstractFloat} <: AbstractSpherical2DGrid
    r::T
    N_theta::Integer
    N_phi::Integer
    coords::Vector{GlobalSphericalPosition{T}}
    areas::Vector{T}
end
function Spherical2DGrid(T::Type{<:AbstractFloat}, r::Real, N_theta::Integer, N_phi::Integer)
    theta = T.(repeat(range(1/N_theta-1, 1-1/N_theta, length=N_theta) * pi, inner=N_phi))
    phi = T.(repeat(range(1/N_phi-1, 1-1/N_phi, length=N_phi) * pi/2, outer=N_theta))
    dtheta, dphi = theta[N_phi+1] - theta[1], phi[2] - phi[1]
    areas = [r^2 * dtheta * (sin(p+dphi/2) - sin(p-dphi/2)) for p in phi]
    return Spherical2DGrid{T}(T(r), N_theta, N_phi, 
                              GlobalSphericalPosition.(T(r), theta, phi), areas)
end
function Spherical2DGrid(r::Real, N_theta::Integer, N_phi::Integer)
    return Spherical2DGrid(typeof(r), r, N_theta, N_phi)
end
function Spherical2DGrid(r::Integer, N_theta::Integer, N_phi::Integer)
    return Spherical2DGrid(promote_type(typeof(r), Float64), r, N_theta, N_phi)
end


"""
    [1] struct Spherical2DGrid_EqSim{T<:AbstractFloat} <: AbstractSpherical2DGrid end
    [2] Spherical2DGrid_EqSim([T::Type,] r::Real, N_theta::Integer, N_phi::Integer)

Global structured grid of surface coordinates (2D) of type `GlobalSphericalPosition{T}` over
the upper hemisphere with radius `r`, assuming equatorial symmetry.

| Field     | Type; with `T<:AbstractFloat`        | Description                       |
|:--------- |:------------------------------------ |:--------------------------------- |
| `r`       | `T`                                  | radius of global body (sphere)    |
| `N_theta` | `Int64`                              | # elements in azimuth direction   |
| `N_phi`   | `Int64`                              | # elements in elevation direction |
| `coords`  | `Vector{GlobalSphericalPosition{T}}` | coordinates                       |
| `areas`   | `Vector{T}`                          | surface area                      |

To ensure expected behavior, the grid object should generally be created with the outer
constructor [2].
"""
struct Spherical2DGrid_EqSim{T<:AbstractFloat} <: AbstractSpherical2DGrid
    r::T
    N_theta::Integer
    N_phi::Integer
    coords::Vector{GlobalSphericalPosition{T}}
    areas::Vector{T}
end
function Spherical2DGrid_EqSim(T::Type{<:AbstractFloat}, r::Real, N_theta::Integer, N_phi::Integer)
    theta = T.(repeat(range(1/N_theta-1, 1-1/N_theta, length=N_theta) * pi, inner=N_phi))
    phi = T.(repeat(range(1/N_phi, 2-1/N_phi, length=N_phi) * pi/4, outer=N_theta))
    dtheta, dphi = theta[N_phi+1] - theta[1], phi[2] - phi[1]
    areas = [r^2 * dtheta * (sin(p+dphi/2) - sin(p-dphi/2)) for p in phi]
    return Spherical2DGrid_EqSim{T}(T(r), N_theta, N_phi, 
                                    GlobalSphericalPosition.(T(r), theta, phi), areas)
end
function Spherical2DGrid_EqSim(r::Real, N_theta::Integer, N_phi::Integer)
    return Spherical2DGrid_EqSim(typeof(r), r, N_theta, N_phi)
end
function Spherical2DGrid_EqSim(r::Integer, N_theta::Integer, N_phi::Integer)
    return Spherical2DGrid_EqSim(promote_type(typeof(r), Float64), r, N_theta, N_phi)
end


"""
    [1] struct Spherical2DGrid_Reduced{T<:AbstractFloat} <: AbstractSpherical2DGrid end
    [2] Spherical2DGrid_Reduced([T::Type,] r::Real, N_phi::Integer)

Global structured grid of surface coordinates (2D) of type `GlobalSphericalPosition{T}` over
the sphere radius `r`. The grid is reduced in the azimuth direction to have approximately 
equal `2*pi*r*cos(phi)/N_theta` grid element lengths.

| Field     | Type; with `T<:AbstractFloat`        | Description                       |
|:--------- |:------------------------------------ |:--------------------------------- |
| `r`       | `T`                                  | radius of global body (sphere)    |
| `N_theta` | `Vector{Int64}`                      | # elements in azimuth direction   |
| `N_phi`   | `Int64`                              | # elements in elevation direction |
| `coords`  | `Vector{GlobalSphericalPosition{T}}` | coordinates                       |
| `areas`   | `Vector{T}`                          | surface area                      |

To ensure expected behavior, the grid object should generally be created with the outer
constructor [2].
"""
struct Spherical2DGrid_Reduced{T<:AbstractFloat} <: AbstractSpherical2DGrid
    r::T
    N_theta::Vector{Int64}
    N_phi::Integer
    coords::Vector{GlobalSphericalPosition{T}}
    areas::Vector{T}
end
function Spherical2DGrid_Reduced(T::Type{<:AbstractFloat}, r::Real, N_phi::Integer)
    dtheta_max = T(pi/N_phi)
    phi0 = T.(range(1/N_phi-1, 1-1/N_phi, length=N_phi) * pi/2)
    
    theta, phi, N_theta = T[], T[], Int64[]
    for p0 in phi0
        n_theta = ceil(Int64, 2*pi*cos(p0)/dtheta_max)
        push!(theta, range(1/n_theta-1, 1-1/n_theta, length=n_theta) .* pi...)
        push!(phi, repeat([p0], n_theta)...)
        push!(N_theta, n_theta)
    end

    areas = T[]
    for i in eachindex(N_theta)
        idx = accumulate(+, N_theta[1:i])[end]
        dtheta, dphi = theta[idx] - theta[idx-1], phi0[2] - phi0[1]
        push!(areas, repeat([r^2 * dtheta * (sin(phi0[i]+dphi/2) - sin(phi0[i]-dphi/2))], N_theta[i])...)
    end

    return Spherical2DGrid_Reduced(T(r), N_theta, N_phi, 
                                   GlobalSphericalPosition.(T(r), theta, phi), areas)
end
function Spherical2DGrid_Reduced(r::Real, N_phi::Integer)
    return Spherical2DGrid_Reduced(typeof(r), r, N_phi)
end
function Spherical2DGrid_Reduced(r::Integer, N_phi::Integer) 
    return Spherical2DGrid_Reduced(promote_type(typeof(r), Float64), r, N_phi)
end


"""
    [1] struct Spherical2DGrid_Reduced_EqSim{T<:AbstractFloat} <: AbstractSpherical2DGrid end
    [2] Spherical2DGrid_Reduced_EqSim([T::Type,] r::Real, N_phi::Integer)

Global structured grid of surface coordinates (2D) of type `GlobalSphericalPosition{T}` over
the upper hemisphere with radius `r`, assuming equatorial symmetry. The grid is reduced in
the azimuth direction to have approximately equal `2*pi*r*cos(phi)/N_theta` grid element 
lengths.

| Field     | Type; with `T<:AbstractFloat`        | Description                       |
|:--------- |:------------------------------------ |:--------------------------------- |
| `r`       | `T`                                  | radius of global body (sphere)    |
| `N_theta` | `Vector{Int64}`                      | # elements in azimuth direction   |
| `N_phi`   | `Int64`                              | # elements in elevation direction |
| `coords`  | `Vector{GlobalSphericalPosition{T}}` | coordinates                       |
| `areas`   | `Vector{T}`                          | surface area                      |

To ensure expected behavior, the grid object should generally be created with the outer
constructor [2].
"""
struct Spherical2DGrid_Reduced_EqSim{T<:AbstractFloat} <: AbstractSpherical2DGrid
    r::T
    N_theta::Vector{Int64}
    N_phi::Integer
    coords::Vector{GlobalSphericalPosition{T}}
    areas::Vector{T}
end
function Spherical2DGrid_Reduced_EqSim(T::Type{<:AbstractFloat}, r::Real, N_phi::Integer)
    dtheta_max = T(pi/2/N_phi)
    phi0 = T.(range(1/N_phi, 2-1/N_phi, length=N_phi) * pi/4)
    
    theta, phi, N_theta = T[], T[], Int64[]
    for p0 in phi0
        n_theta = ceil(Int64, 2*pi*cos(p0)/dtheta_max)
        push!(theta, range(1/n_theta-1, 1-1/n_theta, length=n_theta) .* pi...)
        push!(phi, repeat([p0], n_theta)...)
        push!(N_theta, n_theta)
    end

    areas = T[]
    for i in eachindex(N_theta)
        idx = accumulate(+, N_theta[1:i])[end]
        dtheta, dphi = theta[idx] - theta[idx-1], phi0[2] - phi0[1]
        push!(areas, repeat([r^2 * dtheta * (sin(phi0[i]+dphi/2) - sin(phi0[i]-dphi/2))], N_theta[i])...)
    end

    return Spherical2DGrid_Reduced_EqSim(T(r), N_theta, N_phi, 
                                         GlobalSphericalPosition.(T(r), theta, phi), areas)
end
function Spherical2DGrid_Reduced_EqSim(r::Real, N_phi::Integer) 
    return Spherical2DGrid_Reduced_EqSim(typeof(r), r, N_phi)
end
function Spherical2DGrid_Reduced_EqSim(r::Integer, N_phi::Integer)
    return Spherical2DGrid_Reduced_EqSim(promote_type(typeof(r), Float64), r, N_phi)
end


############################################################################################
#::. UTILITY FUNCTIONS
############################################################################################
function coord2idx(grid::Spherical2DGrid, theta::Real, phi::Real)
    theta = theta == -pi ? pi : theta
    phi = phi == -pi/2 ? -pi/2+eps(Float64) : phi
    idxtheta = ceil(Int64, (theta+pi)*grid.N_theta/2/pi)
    idxphi = ceil(Int64, (phi+pi/2)*grid.N_phi/pi)
    return (idxtheta-1) * grid.N_phi + idxphi
end
function coord2idx(grid::Spherical2DGrid_EqSim, theta::Real, phi::Real)
    theta = theta == -pi ? pi : theta
    phi = phi == 0 ? eps(Float64) : phi
    idxtheta = ceil(Int64, (theta+pi)*grid.N_theta/2/pi)
    idxphi = ceil(Int64, abs(phi)*grid.N_phi*2/pi)
    return (idxtheta-1) * grid.N_phi + idxphi
end
function coord2idx(grid::Spherical2DGrid_Reduced, theta::Real, phi::Real)
    theta = theta == -pi ? pi : theta
    phi = phi == -pi/2 ? -pi/2+eps(Float64) : phi
    idxphi = ceil(Int64, (phi+pi/2)/pi*grid.N_phi)
    idxtheta = ceil(Int64, (theta+pi)/2/pi*grid.N_theta[idxphi])
    if idxphi == 1; return idxtheta; end
    return idxtheta + accumulate(+, grid.N_theta[1:idxphi])[end-1]
end
function coord2idx(grid::Spherical2DGrid_Reduced_EqSim, theta::Real, phi::Real)
    theta = theta == -pi ? pi : theta
    phi = phi == 0 ? eps(Float64) : phi
    idxphi = ceil(Int64, abs(phi)*grid.N_phi*2/pi)
    idxtheta = ceil(Int64, (theta+pi)*grid.N_theta[idxphi]/2/pi)
    if idxphi == 1; return idxtheta; end
    return idxtheta + accumulate(+, grid.N_theta[1:idxphi])[end-1]
end
function coord2idx(grid::AbstractSpherical2DGrid, r::Real, theta::Real, phi::Real)
    return coord2idx(grid, theta, phi)
end


surfacecoords(grid::AbstractSpherical2DGrid) = coords(grid)


volumes(grid::AbstractSpherical2DGrid) = zeros(typeof(grid.r), length(grid.coords))


############################################################################################
#::. EXTENSIONS
############################################################################################
Base.length(grid::AbstractSpherical2DGrid) = length(grid.coords)


Base.size(grid::AbstractSpherical2DGrid) = (grid.N_theta, grid.N_phi)


Base.show(io::IO, ::MIME"text/plain", grid::AbstractSpherical2DGrid) = 
    print(io, "$(typeof(grid)):\n"*
            " r:       $(grid.r)\n"*
            " N_theta: $(grid.N_theta)\n"*
            " N_phi:   $(grid.N_phi)\n"*
            " coords:  $(length(grid.coords))-element $(typeof(coords(grid)))\n"*
            " areas:   $(length(grid.areas))-element $(typeof(grid.areas))\n"*
            "            min: $(min(grid.areas...))\n"*
            "            max: $(max(grid.areas...))\n")

Base.show(io::IO, ::MIME"text/plain", grid::Spherical2DGrid_Reduced) = 
    print(io, "$(typeof(grid)):\n"*
            " r:       $(grid.r)\n"*
            " N_theta: $(length(grid.N_theta))-element Vector{Int64}\n"*
            "            @ quator: $(grid.N_theta[ceil(Int64,grid.N_phi/2)])\n"*
            "            @ poles:  $(grid.N_theta[end])\n"*
            " N_phi:   $(grid.N_phi)\n"*
            " coords:  $(length(grid.coords))-element $(typeof(coords(grid)))\n"*
            " areas:   $(length(grid.areas))-element $(typeof(grid.areas))\n"*
            "            min: $(min(grid.areas...))\n"*
            "            max: $(max(grid.areas...))\n")

Base.show(io::IO, ::MIME"text/plain", grid::Spherical2DGrid_Reduced_EqSim) = 
    print(io, "$(typeof(grid)):\n"*
            " r:       $(grid.r)\n"*
            " N_theta: $(length(grid.N_theta))-element Vector{Int64}\n"*
            "            @ quator: $(grid.N_theta[1])\n"*
            "            @ poles:  $(grid.N_theta[end])\n"*
            " N_phi:   $(grid.N_phi)\n"*
            " coords:  $(length(grid.coords))-element $(typeof(coords(grid)))\n"*
            " areas:   $(length(grid.areas))-element $(typeof(grid.areas))\n"*
            "            min: $(min(grid.areas...))\n"*
            "            max: $(max(grid.areas...))\n") 

            
############################################################################################
#::. EXPORTS
############################################################################################
export 
    Spherical2DGrid,
    Spherical2DGrid_EqSim, 
    Spherical2DGrid_Reduced,
    Spherical2DGrid_Reduced_EqSim