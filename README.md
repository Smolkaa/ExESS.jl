# ExESS.jl

<img src='res/exess_logo.svg' alt="logo" align="right" width = "20%" height="20%">

**Extraterrestrial Exosphere and Surface Simulations** 
</br>
_by A. Peschel ([a.peschel@tum.de](mailto:a.peschel@tum.de))_

Scientific code library for simulating surface-bounded exosphere environments of airless 
bodies in the solar system. Note that the package is under development and subject to 
frequent changes. For questions, please contact the A. Peschel via email.

The documentation and manual can be found in the [GitHub wiki](https://github.com/Smolkaa/ExESS.jl/wiki).


## Installation & Usage

Please refer to the [official guide](https://julialang.org/downloads/platform/) to 
installing Julia on your machine. The `ExESS.jl` package can be installed using Julia's 
package-manager `Pkg`:
```julia
using Pkg
Pkg.add(url="https://github.com/Smolkaa/ExESS.jl.git")
using ExESS
```
If you want to use the latest development version, you can install the package as:
```julia
Pkg.add(url="https://github.com/Smolkaa/ExESS.jl.git", rev="dev")
```
Note that this version may contain bugs and is not recommended for general use.

