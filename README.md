# ExESS.jl

<img src='res/exess_logo.svg' alt="logo" align="right" width = "20%" height="20%">

**Extraterrestrial Exosphere and Surface Simulations** 
</br>
_by A. Smolka ([a.smolka@tum.de](mailto:a.smolka@tum.de))_

Scientific code library for simulating surface-bounded exosphere environments of airless 
bodies in the solar system. Note that the package is under development and subject to 
frequent changes. For questions, please contact the A. Smolka via email.



## Installation

Please refer to the [official guide](https://julialang.org/downloads/platform/) to 
installing Julia on your machine. The `ExESS.jl` package can be installed using Julia's 
package-manager `Pkg`:
```julia
using Pkg
Pkg.add("https://github.com/Smolkaa/ExESS.jl.git")
```
Afterwards, the package can be loaded, as usual, using the `using` command:
```julia
using ExESS
```

## Documentation & Manual

<details>

<summary><b>1. Custom Coordinates</b></summary>

The `ExESS.jl` package provides a set of custom types and functions to work with different types of coordinates. While the library can also be used with Julia types like `Tuple`, the custom types provide a more streamlined way to work with coordinates, and, more importantly, coordinate transformation.

The custom coordinates can be divided into different groups, namely position and velocity types, cartesian and spherical coordinates, and global and local coordinates:

```
ExESS
└─ AbstractCoordinates
    ├─ AbstractPosition
    │  ├─ GlobalCartesianPosition{T<:AbstractFloat}
    │  ├─ LocalCartesianPosition{T<:AbstractFloat}
    │  └─ GlobalSphericalPosition{T<:AbstractFloat}
    │
    └─ AbstractVelocity
       ├─ GlobalCartesianVelocity{T<:AbstractFloat}
       ├─ LocalCartesianVelocity{T<:AbstractFloat}
       └─ GlobalSphericalVelocity{T<:AbstractFloat}
```

Global coordinate systems in this instance have their origin in the center of the respective
central body of the exosphere, for example, the Moon. While the orientation of the coordinate
axes is not generally fixed for all applications, the library assumes that for cartesian
systems, the global x-axis is pointing towards the subsolar point and, thus, the Sun. This 
means, for instance, that for a [`GlobalCartesianPosition`](@ref) the subsolar point on the
Moon is
```julia
gcp_subsol = GlobalCartesianPosition(LUNAR_RADIUS, 0, 0) # subsolar point
```
The global z-axis is pointing towards positive longitude and the global y-axis is pointing 
towards positive latitude, completing the regular three-dimensional system. A spherical
global system is based on the same principle, where the global angular arguments are zero
at the subsolar points, resulting in
```julia
gsp_subsol = GlobalSphericalPosition(LUNAR_RADIUS, 0, 0) # subsolar point
```
The global azimuth angle is rotating around the positive longitude, and the global zenith 
angle is rotating with the positive latitude. In terms of global velocities, their coordinate 
systems are identical to the ones described above. The subsolar reference frame, including
the main axes from the global cartesian system, is shown in the figure below:

![global_cartesian_subsolar_cs.png](img/global_cartesian_subsolar_cs.png)

Local coordinate systems have, as their name suggests, their origin at a specific point in
three-dimensional space. From this origin, the local cartesian system can be constructed as
follows. The x-axis is pointing towards positive longitude (usually east), the y-axis is
pointing towards positive latitude (usually north), and the z-axis is extending the line 
from the local origin to the center of the central body. Thus, at the subsolar point, the
following two velocities would be identical:
```julia
lcv_subsol = LocalCartesianVelocity(0, 0, 100) # local velocity
gcv_subsol = GlobalCartesianVelocity(100, 0, 0) # global velocity
```

Please note that these custom coordinates can be used differently depending on the application, 
though all the integrated functionality of this library assumes, where
necessary, the coordinate system as described above.

---

</details>



<details>

<summary><b>2. Drivers</b></summary>

<details>

<summary><b>2.1 Maxwellian Velocity Distributions</b></summary>

With temperature being one of the main drivers of exospheres, Maxwellian distributions 
are a common way to describe velocity and energy distributions for thermal processes. The
`ExESS.jl` package provides both Maxwell-Boltzmann distribution as well as flux-weighted 
Maxwell-Boltzmann distributions for azimuth, elevation, speed, and velocity.

## Creating a Maxwellian Distribution

All Maxwellian distributions are created based on two parameters: the temperature `T` and 
the mass `m` of the particles. The temperature is given in Kelvin, and the mass is given in
kilograms. 

!!! note
    To convert particle masses from amu to kg, the utility function [`amu2kg`](@ref) can be
    used.

Here are some examples of the creation of Maxwellian distributions:

```julia
# Maxwell-Boltzmann Azimuth distribution @300K for a 4amu particle (He)
MBAzimuthDistribution(300, amu2kg(4))

# Maxwell-Boltzmann Flux Elevation distribution @250K for a 17amu particle (OH)
MBFluxElevationDistribution(250, amu2kg(17))

# Maxwell-Boltzmann Speed distribution @200K for a 2amu particle (H2)
MBSpeedDistribution(200, amu2kg(2))

# Maxwell-Boltzmann Velocity distribution @400K for a 18amu particle (H2O)
MBFluxVelocityDistribution(400, amu2kg(18))
```

## Sampling a distribution

Every distribution can be sampled using the `rand` function from Julia's `Base` module.

!!! note
    Note that it is generally recommended to create the distribution inside of the `rand` 
    call due to performance reasons.

```@repl exess
rand(MBAzimuthDistribution(300, amu2kg(4)))
rand(MBFluxElevationDistribution(250, amu2kg(17)))
rand(MBSpeedDistribution(200, amu2kg(2)))
rand(MBFluxVelocityDistribution(400, amu2kg(18)))
```

Please note that the velocity samples are returned in a local cartesian coordinate system.
For more information on the respective coordinate systems, see the 
[coordinates documentation](../base/coordinates.md). The `rand` call can also take a `Type` 
argument as its first parameter which automatically converts the sample into the desired 
type. Here, we could, for example, directly use the correct coordinate type 
[`LocalCartesianVelocity`](@ref) to get the samples in the correct coordinate system:

```@repl exess
rand(LocalCartesianVelocity, MBFluxVelocityDistribution(400, amu2kg(18)))
```

!!! note
    Note that passing a different velocity type will result in an error, providing an 
    additional safety net for the user:

    ```@repl exess
    rand(GlobalCartesianVelocity, MBFluxVelocityDistribution(400, amu2kg(18)))
    ```

</details>

---

</details>
