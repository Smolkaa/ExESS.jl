############################################################################################
#::. TYPES
############################################################################################
abstract type AbstractSpherical3DGrid <: AbstractSphericalGrid end


############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    [1] struct Spherical3DGrid{T<:AbstractFloat} <: AbstractSpherical3DGrid
    [2] Spherical3DGrid([T::Type{<:AbstractFloat}, ] r0::Real, h::AbstractVector, 
                        N_theta::Integer, N_phi::Integer)
    [3] Spherical3DGrid([T::Type{<:AbstractFloat}, ] r::AbstractVector, N_theta::Integer, 
                        N_phi::Integer)

Global structured volume grid (3D) of type `GlobalSphericalPosition{T}` over a sphere of
radius `r0` with heights of the individual radial layers `h`. Alternatively, the radial 
distances can directly be specified by `r` (`r = r0 .+ h`). Note that all heights have
to be positive.

| Field     | Type; with `T<:AbstractFloat`        | Description                       |
|:--------- |:------------------------------------ |:--------------------------------- |
| `r0`      | `T`                                  | radius of global body (sphere)    |
| `h`       | `Vector{T}`                          | heights above radial base `r0`    |
| `N_r`     | `Int64`                              | # elements in radial direction    |
| `N_theta` | `Int64`                              | # elements in azimuth direction   |
| `N_phi`   | `Int64`                              | # elements in elevation direction |
| `coords`  | `Vector{GlobalSphericalPosition{T}}` | coordinates                       |
| `areas`   | `Vector{T}`                          | surface area                      |
| `volumes` | `Vector{T}`                          | volumes                           |

To ensure expected behavior, the grid object should generally be created with the outer
constructors [2] or [3].
"""
struct Spherical3DGrid{T<:AbstractFloat} <: AbstractSpherical3DGrid
    r0::T
    h::Vector{T}
    N_r::Int64
    N_theta::Int64
    N_phi::Int64
    coords::Vector{GlobalSphericalPosition{T}}
    areas::Vector{T}
    volumes::Vector{T}
end
function Spherical3DGrid(T::Type{<:AbstractFloat}, r0::Real, h::AbstractVector, 
            N_theta::Integer, N_phi::Integer)
    @assert sort(h)[1] > 0
    h = sort(h)
    N_r = length(h)

    # coordinates
    r = T.(repeat(r0 .+ vcat(0, h[1:end-1]), inner=N_theta*N_phi))
    theta = T.(repeat(range(1/N_theta-1, 1-1/N_theta, length=N_theta) .* pi, inner=N_phi, outer=N_r))
    phi = T.(repeat(range(1/N_phi-1, 1-1/N_phi, length=N_phi) .* pi/2, outer=N_theta*N_r))
    
    # area and volume calculation
    dR, dtheta, dphi = vcat(h[1], diff(h)), theta[N_phi+1] - theta[1], phi[2] - phi[1]
    areas = [r[1]^2 * dtheta * (sin(p+dphi/2) - sin(p-dphi/2)) for p in phi[1:N_theta*N_phi]]
    for rr in [r0 .+ vcat(h[1:end-1])...]
        push!(areas, areas[1:N_theta*N_phi] / r0^2 * rr^2...)
    end
    volumes = [
        ((r[i] + dr)^3 - r[i]^3)/3 * dtheta * (sin(phi[i]+dphi/2) - sin(phi[i]-dphi/2))
        for i in 1:N_r*N_theta*N_phi for dr = dR[Int64(1 + floor((i-1) / (N_theta*N_phi)))]
    ]

    return Spherical3DGrid(T(r0), T.(h), N_r, N_theta, N_phi, 
        GlobalSphericalPosition.(r, theta, phi), T.(areas), T.(volumes))
end
function Spherical3DGrid(r0::Real, h::AbstractVector, N_theta::Integer, N_phi::Integer)
    return Spherical3DGrid(typeof(r0), r0, h, N_theta, N_phi)
end
function Spherical3DGrid(r0::Integer, h::AbstractVector, N_theta::Integer, N_phi::Integer)
    return Spherical3DGrid(promote_type(typeof(r0), Float64), r0, h, N_theta, N_phi)
end
function Spherical3DGrid(T::Type{<:AbstractFloat}, r::AbstractVector, N_theta::Integer, 
                         N_phi::Integer)
    r0, h = r[1], accumulate(+,diff(r))
    return Spherical3DGrid(T, r0, h, N_theta, N_phi)
end
function Spherical3DGrid(r::AbstractVector, N_theta::Integer, N_phi::Integer)
    return Spherical3DGrid(typeof(r[1]), r, N_theta, N_phi)
end
function Spherical3DGrid(r::AbstractVector{<:Integer}, N_theta::Integer, N_phi::Integer)
    return Spherical3DGrid(promote_type(typeof(r[1]), Float64), r, N_theta, N_phi)
end

"""
    [1] struct Spherical3DGrid_EqSim{T<:AbstractFloat} <: AbstractSpherical3DGrid end
    [2] Spherical3DGrid_EqSim([T::Type{<:AbstractFloat}, ] r0::Real, h::AbstractVector, 
                              N_theta::Integer, N_phi::Integer)
    [3] Spherical3DGrid_EqSim([T::Type{<:AbstractFloat}, ] r::AbstractVector, 
                              N_theta::Integer, N_phi::Integer)

Global structured volume grid (3D) of type `GlobalSphericalPosition{T}` over a hemisphere of
radius `r0` with heights of the individual radial layers `h`. Alternatively, the radial 
distances can directly be specified by `r` (`r = r0 .+ h`). Note that all heights have
to be positive.

| Field     | Type; with `T<:AbstractFloat`        | Description                       |
|:--------- |:------------------------------------ |:--------------------------------- |
| `r0`      | `T`                                  | radius of global body (sphere)    |
| `h`       | `Vector{T}`                          | heights above radial base `r0`    |
| `N_r`     | `Int64`                              | # elements in radial direction    |
| `N_theta` | `Int64`                              | # elements in azimuth direction   |
| `N_phi`   | `Int64`                              | # elements in elevation direction |
| `coords`  | `Vector{GlobalSphericalPosition{T}}` | coordinates                       |
| `areas`   | `Vector{T}`                          | surface area                      |
| `volumes` | `Vector{T}`                          | volumes                           |

To ensure expected behavior, the grid object should generally be created with the outer
constructors [2] or [3].
"""
struct Spherical3DGrid_EqSim{T<:AbstractFloat} <: AbstractSpherical3DGrid
    r0::T
    h::Vector{T}
    N_r::Int64
    N_theta::Int64
    N_phi::Int64
    coords::Vector{GlobalSphericalPosition{T}}
    areas::Vector{T}
    volumes::Vector{T}
end
function Spherical3DGrid_EqSim(T::Type{<:AbstractFloat}, r0::Real, h::AbstractVector, 
            N_theta::Integer, N_phi::Integer)

    h = sort(h)
    @assert h[1] > 0
    N_r = length(h)

    # coordinates
    r = T.(repeat(r0 .+ vcat(0, h[1:end-1]), inner=N_theta*N_phi))
    theta = T.(repeat(
        range(1/N_theta-1, 1-1/N_theta, length=N_theta) .* pi, inner=N_phi, outer=N_r))
    phi = T.(repeat(range(1/N_phi, 2-1/N_phi, length=N_phi) .* pi/4, outer=N_theta*N_r))
    
    # area and volume calculation
    dR, dtheta, dphi = vcat(h[1], diff(h)), theta[N_phi+1] - theta[1], phi[2] - phi[1]
    areas = [r[1]^2 * dtheta * (sin(p+dphi/2) - sin(p-dphi/2)) for p in phi[1:N_theta*N_phi]]
    for rr in [r0 .+ vcat(h[1:end-1])...]
        push!(areas, areas[1:N_theta*N_phi] / r0^2 * rr^2...)
    end
    volumes = [
        ((r[i] + dr)^3 - r[i]^3)/3 * dtheta * (sin(phi[i]+dphi/2) - sin(phi[i]-dphi/2))
        for i in 1:N_r*N_theta*N_phi for dr = dR[Int64(1 + floor((i-1) / (N_theta*N_phi)))]
    ]

    return Spherical3DGrid_EqSim(T(r0), T.(h), N_r, N_theta, N_phi, 
                                 GlobalSphericalPosition.(r, theta, phi), 
                                 T.(areas), T.(volumes))
end
function Spherical3DGrid_EqSim(r0::Real, h::AbstractVector, N_theta::Integer, N_phi::Integer)
    return Spherical3DGrid_EqSim(typeof(r0), r0, h, N_theta, N_phi)
end
function Spherical3DGrid_EqSim(r0::Integer, h::AbstractVector, N_theta::Integer, N_phi::Integer)
    return Spherical3DGrid_EqSim(promote_type(typeof(r0), Float64), r0, h, N_theta, N_phi)
end
function Spherical3DGrid_EqSim(T::Type{<:AbstractFloat}, r::AbstractVector, N_theta::Integer, 
                               N_phi::Integer)
    r0, h = r[1], accumulate(+,diff(r))
    return Spherical3DGrid_EqSim(T, r0, h, N_theta, N_phi)
end
function Spherical3DGrid_EqSim(r::AbstractVector, N_theta::Integer, N_phi::Integer)
    return Spherical3DGrid_EqSim(typeof(r[1]), r, N_theta, N_phi)
end
function Spherical3DGrid_EqSim(r::AbstractVector{<:Integer}, N_theta::Integer, N_phi::Integer)
    return Spherical3DGrid_EqSim(promote_type(typeof(r[1]), Float64), r, N_theta, N_phi)
end


"""
    [1] struct Spherical3DGrid_Reduced{T<:AbstractFloat} <: AbstractSpherical3DGrid end
    [2] Spherical3DGrid_Reduced([T::Type{<:AbstractFloat}, ] r0::Real, h::AbstractVector, 
                                N_phi::Integer)
    [3] Spherical3DGrid_Reduced([T::Type{<:AbstractFloat}, ] r::AbstractVector, 
                                N_phi::Integer)

Global structured volume grid (3D) of type `GlobalSphericalPosition{T}` over a sphere of
radius `r0` with heights of the individual radial layers `h`. Alternatively, the radial 
distances can directly be specified by `r` (`r = r0 .+ h`). Note that all heights have
to be positive. The grid is reduced in the azimuth direction to have approximately equal 
`2*pi*r*cos(phi)/N_theta` grid element lengths.

| Field     | Type; with `T<:AbstractFloat`        | Description                       |
|:--------- |:------------------------------------ |:--------------------------------- |
| `r0`      | `T`                                  | radius of global body (sphere)    |
| `h`       | `Vector{T}`                          | heights above radial base `r0`    |
| `N_r`     | `Int64`                              | # elements in radial direction    |
| `N_theta` | `Int64`                              | # elements in azimuth direction   |
| `N_phi`   | `Int64`                              | # elements in elevation direction |
| `coords`  | `Vector{GlobalSphericalPosition{T}}` | coordinates                       |
| `areas`   | `Vector{T}`                          | surface area                      |
| `volumes` | `Vector{T}`                          | volumes                           |

To ensure expected behavior, the grid object should generally be created with the outer
constructors [2] or [3].
"""
struct Spherical3DGrid_Reduced{T<:AbstractFloat} <: AbstractSpherical3DGrid
    r0::T
    h::Vector{T}
    N_r::Int64
    N_theta::Vector{Int64}
    N_phi::Int64
    coords::Vector{GlobalSphericalPosition{T}}
    areas::Vector{T}
    volumes::Vector{T}
end
function Spherical3DGrid_Reduced(T::Type{<:AbstractFloat}, r0::Real, 
            h::AbstractVector, N_phi::Integer)

    h = sort(h)
    @assert h[1] > 0
    N_r = length(h)

    # prepare coordiante calculation
    dtheta_max = T(pi/N_phi)
    phi0 = T.(range(1/N_phi-1, 1-1/N_phi, length=N_phi) * pi/2)

    # coordinates
    r, theta, phi, N_theta = T[], T[], T[], Int64[]
    for p0 in phi0
        n_theta = ceil(Int64, 2*pi*cos(p0)/dtheta_max)
        push!(theta, range(1/n_theta-1, 1-1/n_theta, length=n_theta) .* pi...)
        push!(phi, repeat([p0], n_theta)...)
        push!(N_theta, n_theta)
    end
    r = T.(repeat(r0 .+ vcat(0, h[1:end-1]), inner=sum(N_theta)))
    theta = repeat(theta, outer=N_r)
    phi = repeat(phi, outer=N_r)
    
    # areas calculation
    areas = T[]
    for i in 1:N_phi
        idx = accumulate(+, N_theta[1:i])[end]
        dtheta, dphi = theta[idx] - theta[idx-1], phi0[2] - phi0[1]
        push!(areas, repeat(
            [r[1]^2 * dtheta * (sin(phi0[i]+dphi/2) - sin(phi0[i]-dphi/2))], N_theta[i])...)
    end
    for rr in [r0 .+ vcat(h[1:end-1])...]
        push!(areas, areas[1:sum(N_theta)] / r0^2 * rr^2...)
    end

    # volumes calculation
    volumes = T[]
    for i in 1:N_r, j in 1:N_phi, _ in 1:N_theta[j]
        ri, dr = [r0 .+ vcat(0, h[1:end-1])...][i], vcat(h[1], diff(h))[i]
        dtheta = 2pi / N_theta[j]
        p0, dphi = phi0[j], pi / N_phi
        push!(volumes, ((ri + dr)^3 - ri^3)/3 * dtheta * (sin(p0+dphi/2) - sin(p0-dphi/2)))
    end
    
    return Spherical3DGrid_Reduced(T(r0), T.(h), N_r, N_theta, N_phi, 
                                   GlobalSphericalPosition.(r, theta, phi), 
                                   T.(areas), T.(volumes))
end
function Spherical3DGrid_Reduced(r0::Real, h::AbstractVector, N_phi::Integer)
    return Spherical3DGrid_Reduced(typeof(r0), r0, h, N_phi)
end
function Spherical3DGrid_Reduced(r0::Integer, h::AbstractVector, N_phi::Integer)
    return Spherical3DGrid_Reduced(promote_type(typeof(r0), Float64), r0, h, N_phi)
end
function Spherical3DGrid_Reduced(T::Type{<:AbstractFloat}, r::AbstractVector, N_phi::Integer)
    r0, h = r[1], accumulate(+,diff(r))
    Spherical3DGrid_Reduced(T, r0, h, N_phi)
end
function Spherical3DGrid_Reduced(r::AbstractVector, N_phi::Integer)
    return Spherical3DGrid_Reduced(typeof(r[1]), r, N_phi)
end
function Spherical3DGrid_Reduced(r::AbstractVector{<:Integer}, N_phi::Integer)
    return Spherical3DGrid_Reduced(promote_type(typeof(r[1]), Float64), r, N_phi)
end


"""
    [1] struct Spherical3DGrid_Reduced_EqSim{T<:AbstractFloat} <: AbstractSpherical3DGrid
    [2] Spherical3DGrid_Reduced_EqSim([T::Type{<:AbstractFloat},] r0::Real, 
            h::AbstractVector, N_phi::Integer)
    [3] Spherical3DGrid_Reduced_EqSim([T::Type{<:AbstractFloat},] r::AbstractVector, 
            N_phi::Integer)

Global structured volume grid (3D) of type `GlobalSphericalPosition{T}` over a sphere of
radius `r0` with heights of the individual radial layers `h`. Alternatively, the radial 
distances can directly be specified by `r` (`r = r0 .+ h`). Note that all heights have
to be positive. The grid is reduced in the azimuth direction to have approximately equal 
`2*pi*r*cos(phi)/N_theta` grid element lengths.

| Field     | Type; with `T<:AbstractFloat`        | Description                       |
|:--------- |:------------------------------------ |:--------------------------------- |
| `r0`      | `T`                                  | radius of global body (sphere)    |
| `h`       | `Vector{T}`                          | heights above radial base `r0`    |
| `N_r`     | `Int64`                              | # elements in radial direction    |
| `N_theta` | `Int64`                              | # elements in azimuth direction   |
| `N_phi`   | `Int64`                              | # elements in elevation direction |
| `coords`  | `Vector{GlobalSphericalPosition{T}}` | coordinates                       |
| `areas`   | `Vector{T}`                          | surface area                      |
| `volumes` | `Vector{T}`                          | volumes                           |

To ensure expected behavior, the grid object should generally be created with the outer
constructors [2] or [3].
"""
struct Spherical3DGrid_Reduced_EqSim{T<:AbstractFloat} <: AbstractSpherical3DGrid
    r0::T
    h::Vector{T}
    N_r::Int64
    N_theta::Vector{Int64}
    N_phi::Int64
    coords::Vector{GlobalSphericalPosition{T}}
    areas::Vector{T}
    volumes::Vector{T}
end
function Spherical3DGrid_Reduced_EqSim(T::Type{<:AbstractFloat}, r0::Real, 
            h::AbstractVector, N_phi::Integer)

    h = sort(h)
    @assert h[1] > 0
    N_r = length(h)

    # prepare coordiante calculation
    dtheta_max = T(pi/2/N_phi)
    phi0 = T.(range(1/N_phi, 2-1/N_phi, length=N_phi) * pi/4)

    # coordinates
    r, theta, phi, N_theta = T[], T[], T[], Int64[]
    for p0 in phi0
        n_theta = ceil(Int64, 2*pi*cos(p0)/dtheta_max)
        push!(theta, range(1/n_theta-1, 1-1/n_theta, length=n_theta) .* pi...)
        push!(phi, repeat([p0], n_theta)...)
        push!(N_theta, n_theta)
    end
    r = T.(repeat(r0 .+ vcat(0, h[1:end-1]), inner=sum(N_theta)))
    theta = repeat(theta, outer=N_r)
    phi = repeat(phi, outer=N_r)
    
    # areas calculation
    areas = T[]
    for i in 1:N_phi
        idx = accumulate(+, N_theta[1:i])[end]
        dtheta, dphi = theta[idx] - theta[idx-1], phi0[2] - phi0[1]
        push!(areas, repeat(
            [r[1]^2 * dtheta * (sin(phi0[i]+dphi/2) - sin(phi0[i]-dphi/2))], N_theta[i])...)
    end
    for rr in [r0 .+ vcat(h[1:end-1])...]
        push!(areas, areas[1:sum(N_theta)] / r0^2 * rr^2...)
    end

    # volumes calculation
    volumes = T[]
    for i in 1:N_r, j in 1:N_phi, _ in 1:N_theta[j]
        ri, dr = [r0 .+ vcat(0, h[1:end-1])...][i], vcat(h[1], diff(h))[i]
        dtheta = 2pi / N_theta[j]
        p0, dphi = phi0[j], pi / N_phi / 2
        push!(volumes, ((ri + dr)^3 - ri^3)/3 * dtheta * (sin(p0+dphi/2) - sin(p0-dphi/2)))
    end
    
    return Spherical3DGrid_Reduced_EqSim(T(r0), T.(h), N_r, N_theta, N_phi, 
                                         GlobalSphericalPosition.(r, theta, phi), 
                                         T.(areas), T.(volumes))
end
function Spherical3DGrid_Reduced_EqSim(r0::Real, h::AbstractVector, N_phi::Integer)
    return Spherical3DGrid_Reduced_EqSim(typeof(r0), r0, h, N_phi)
end
function Spherical3DGrid_Reduced_EqSim(r0::Integer, h::AbstractVector, N_phi::Integer)
    return Spherical3DGrid_Reduced_EqSim(promote_type(typeof(r0), Float64), r0, h, N_phi)
end
function Spherical3DGrid_Reduced_EqSim(T::Type{<:AbstractFloat}, r::AbstractVector, 
            N_phi::Integer)
    r0, h = r[1], accumulate(+,diff(r))
    Spherical3DGrid_Reduced_EqSim(T, r0, h, N_phi)
end
function Spherical3DGrid_Reduced_EqSim(r::AbstractVector, N_phi::Integer)
    return Spherical3DGrid_Reduced_EqSim(typeof(r[1]), r, N_phi)
end
function Spherical3DGrid_Reduced_EqSim(r::AbstractVector{<:Integer}, N_phi::Integer)
    return Spherical3DGrid_Reduced_EqSim(promote_type(typeof(r[1]), Float64), r, N_phi)
end


############################################################################################
#::. UTILITY FUNCTIONS
############################################################################################
function coord2idx(grid::Spherical3DGrid, r::Real, theta::Real, phi::Real)
    if !(grid.r0 <= r <= grid.r0 + grid.h[end]); return 0; end

    theta    = theta == -pi ? pi : theta
    phi      = phi == -pi/2 ? -pi/2+eps(Float64) : phi
    idxr     = findfirst(grid.h .+ grid.r0 .> r)
    idxtheta = ceil(Int64, (theta+pi)*grid.N_theta/2/pi)
    idxphi   = ceil(Int64, (phi+pi/2)*grid.N_phi/pi)
    return (idxr-1) * grid.N_phi*grid.N_theta + (idxtheta-1) * grid.N_phi + idxphi
end
function coord2idx(grid::Spherical3DGrid_EqSim, r::Real, theta::Real, phi::Real)
    if !(grid.r0 <= r <= grid.r0 + grid.h[end]); return 0; end

    theta    = theta == -pi ? pi : theta
    phi      = phi == 0 ? eps(Float64) : phi
    idxr     = findfirst(grid.h .+ grid.r0 .> r)
    idxtheta = ceil(Int64, (theta+pi)*grid.N_theta/2/pi)
    idxphi   = ceil(Int64, abs(phi)*grid.N_phi*2/pi)
    return (idxr-1) * grid.N_phi*grid.N_theta + (idxtheta-1) * grid.N_phi + idxphi
end
function coord2idx(grid::Spherical3DGrid_Reduced, r::Real, theta::Real, phi::Real)
    if !(grid.r0 <= r <= grid.r0 + grid.h[end]); return 0; end

    theta    = theta == -pi ? pi : theta
    phi      = phi == -pi/2 ? -pi/2+eps(Float64) : phi
    idxr     = findfirst(grid.h .+ grid.r0 .> r)
    idxphi   = ceil(Int64, (phi+pi/2)/pi*grid.N_phi)
    idxtheta = ceil(Int64, (theta+pi)/2/pi*grid.N_theta[idxphi])
    if idxphi == 1; return (idxr-1) * sum(grid.N_theta) + idxtheta; end
    return (idxr-1) * sum(grid.N_theta) + idxtheta + accumulate(+, grid.N_theta[1:idxphi])[end-1]
end
function coord2idx(grid::Spherical3DGrid_Reduced_EqSim, r::Real, theta::Real, phi::Real)
    if !(grid.r0 <= r <= grid.r0 + grid.h[end]); return 0; end
    
    theta    = theta == -pi ? pi : theta
    phi      = phi == 0 ? eps(Float64) : phi
    idxr     = findfirst(grid.h .+ grid.r0 .> r)
    idxphi   = ceil(Int64, abs(phi)*grid.N_phi*2/pi)
    idxtheta = ceil(Int64, (theta+pi)/2/pi*grid.N_theta[idxphi])
    if idxphi == 1; return (idxr-1) * sum(grid.N_theta) + idxtheta; end
    return (idxr-1) * sum(grid.N_theta) + idxtheta + accumulate(+, grid.N_theta[1:idxphi])[end-1]
end


surfacecoords(grid::AbstractSpherical3DGrid) = grid.coords[1:grid.N_theta*grid.N_phi]
surfacecoords(grid::Spherical3DGrid_Reduced) = grid.coords[1:sum(grid.N_theta)]
surfacecoords(grid::Spherical3DGrid_Reduced_EqSim) = grid.coords[1:sum(grid.N_theta)]


############################################################################################
#::. EXTENSIONS
############################################################################################
Base.length(grid::AbstractSpherical3DGrid) = length(grid.coords)


Base.size(grid::AbstractSpherical3DGrid) = (grid.N_r, grid.N_theta, grid.N_phi)


Base.show(io::IO, ::MIME"text/plain", grid::AbstractSpherical3DGrid) = 
    print(io, "$(typeof(grid)):\n"*
            " r0:      $(grid.r0)\n"*
            " h:       $([h for h in grid.h])\n"*
            " N_r:     $(grid.N_r)\n"*
            " N_theta: $(grid.N_theta)\n"*
            " N_phi:   $(grid.N_phi)\n"*
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
            " N_theta: $(length(grid.N_theta))-element Vector{Int64}\n"*
            "            @ quator: $(grid.N_theta[1])\n"*
            "            @ poles:  $(grid.N_theta[end])\n"*
            " N_phi:   $(grid.N_phi)\n"*
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