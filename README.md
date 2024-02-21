# ExESS.jl: Extraterrestrial Exosphere and Surfaces Simulations 

<img src='res/exess_logo.svg' alt="logo" align="right" width = "20%" height="20%">

Some information!


## Installation

This simulation tool was created as a custom package for the computer language Julia. Please refer to the [official guide](https://julialang.org/downloads/platform/) to installing Julia on your machine. The entire module can be included directly by specifying the path of the downloaded package in the `include` call.
```julia
include(joinpath(PathToPackage, "src", "ExESS.jl"))
using .ExESS
```


## Required Packages

The `ExESS` package relies on multiple other packages that have to be added and installed by the user. They can be added using Julia's package-manager `Pkg`:
```julia
using Pkg
Pkg.add("Interpolations")
Pkg.add("NearestNeighbors")
Pkg.add("DifferentialEquations")
Pkg.add("SpecialFunctions")
Pkg.add("StaticArrays")
```


## Notes

* all internal functions (not exported) are prefixed with an underscore `_`
* abstract types are prefixed with `Abstract` and are **not** exported

* rules for changing the version number in `Project.toml`
  - if the release is breaking, increment MAJOR
  - if the release adds a new user-visible feature, increment MINOR
  - otherwise (bug fixes, documentation improvements), increment PATCH


## ToDo's until next Version

- add documentation to the repository using `Documenter.jl`
- add references to all applicable functions


## Author Information

Alexander Smolka - [alexander.smolka@tum.de](mailto:a.smolka@tum.de)
