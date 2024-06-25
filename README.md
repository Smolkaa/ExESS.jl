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

## Documentation

<details>

<summary><b>1. Custom Coordinates</b></summary>

---

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
