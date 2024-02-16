############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    [1] lunar_surface_temperatures_HURLEY2015(theta::Real, phi::Real)
    [2] lunar_surface_temperatures_HURLEY2015(thetas::AbstractVector, phis::AbstractVector)
    [3] lunar_surface_temperatures_HURLEY2015(xs::GlobalSphericalPosition)
    [4] lunar_surface_temperatures_HURLEY2015(XS::Vector{GlobalSphericalPosition})
    [5] lunar_surface_temperatures_HURLEY2015(grid::AbstractGrid)

Calculates the lunar surface temperatures based on the analytic formula given in Hurley et
al. 2015. The input parameters are in spherical, sub-solar coordinates, with the longitude
`theta` in the range [-π, π] and the latitude `phi` in the range [-π/2, π/2]. Alternatively,
the input can be a vector of `GlobalSphericalPosition` objects or an `AbstractGrid` object.

**Reference(s)**

- Hurley et al. (2015), "An analytic function of lunar surface temperature for exospheric
  modeling", DOI: 10.1016/j.icarus.2014.08.043
"""
function lunar_surface_temperatures_HURLEY2015(theta::S, phi::S) where {S<:AbstractFloat}
    @assert -pi <= theta <= pi "Longitude must be in [-pi, pi]!"
    @assert -pi/2 <= phi <= pi/2 "Latitude must be in [-pi/2, pi/2]!"   

    if abs(theta) >= pi/2
        a = [444.738, -448.937, 239.668, -63.8844, 8.34064, -0.423502]
        if theta < 0; theta += 2pi; end
        cophi = -(phi - pi/2)
        return S(sum([a[i] * theta^(i-1) for i in 1:6]) + 35 * (sin(cophi)-1))
    end
    return S(262*(cos(theta) * cos(phi))^(1/2) + 130)
end
function lunar_surface_temperatures_HURLEY2015(theta::Real, phi::Real) 
    return lunar_surface_temperatures_HURLEY2015(promote(theta, phi)...)
end
function lunar_surface_temperatures_HURLEY2015(theta::Integer, phi::Integer) 
    return lunar_surface_temperatures_HURLEY2015(promote(theta, phi, 1.0)[1:2]...)
end
function lunar_surface_temperatures_HURLEY2015(thetas::AbstractVector, phis::AbstractVector)
    return lunar_surface_temperatures_HURLEY2015.(thetas, phis)
end
function lunar_surface_temperatures_HURLEY2015(xs::GlobalSphericalPosition)
    return lunar_surface_temperatures_HURLEY2015(xs.theta, xs.phi)
end
function lunar_surface_temperatures_HURLEY2015(XS::Vector{GlobalSphericalPosition{S}}) where {S}
    return lunar_surface_temperatures_HURLEY2015.(XS)
end
function lunar_surface_temperatures_HURLEY2015(grid::AbstractGrid)
    return lunar_surface_temperatures_HURLEY2015(surfacecoords(grid))
end


############################################################################################
#::. EXPORTS
############################################################################################
export 
    lunar_surface_temperatures_HURLEY2015