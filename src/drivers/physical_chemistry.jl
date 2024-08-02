############################################################################################
#::. FUNCTIONS
############################################################################################
"""
    [1] arrhenius([S::Type{<:AbstractFloat}], [A::Real], T::Real, Ea::Real)

Calculation of the rate constant based on the
[Arrhenius Equation](https://en.wikipedia.org/wiki/Arrhenius_equation), with parameters `A`
as the pre-exponential factor in (s-1), `T` as the temperature in (K), and `Ea` as the
activation energy in (J). The function returns the rate constant in (s-1). Additionally,
the function takes an optional type parameter `S` to specify the output type of the function.

The pre-exponential factor `A` is also an optional input. If not given, the function will
calculate the pre-exponential factor based on `A = k_B * T / h`, with Boltzmann's constant
`k_B`, and Planck's constant `h`.
"""
function arrhenius(A::S, T::S, Ea::S) where {S<:AbstractFloat}
    return S(A * exp(- Ea / (BOLTZMANN_CONSTANT * T)))
end
arrhenius(A::Real, T::Real, Ea::Real) = arrhenius(promote(A, T, Ea)...)
arrhenius(A::Integer, T::Integer, Ea::Integer) = arrhenius(promote(A, T, Ea, 1.0)[1:3]...)
function arrhenius(T::S, Ea::S) where {S<:AbstractFloat}
    return arrhenius(S(BOLTZMANN_CONSTANT * T / PLANCK_CONSTANT), T, Ea)
end
arrhenius(T::Real, Ea::Real) = arrhenius(promote(T, Ea)...)
arrhenius(T::Integer, Ea::Integer) = arrhenius(promote(T, Ea, 1.0)[1:2]...)
arrhenius(S::Type{<:AbstractFloat}, args...) = S(arrhenius(args...))



"""
    [1] diffusion_time([S::Type{<:AbstractFloat}], h::Real, D0::Real, T::Real, E::Real)

Calculates the time scale of diffusive transport of a species over a distance of `h` in (m),
based on a diffusion coefficient `D0` in (m2/s), temperature `T` in (K), and diffusion
energy `E` in (J). The function returns the time scale in (s). The calculation is based on
the Arrhenius law. See the function `arrhenius` for more information.

The optional type parameter `S` is changing the output type of the function.

**References**

* Grumpe, A., Wöhler, C., Berezhnoy, A. A., & Shevchenko, V. V. (2019).
  Time-of-day-dependent behavior of surficial lunar hydroxyl/water: Observations and
  modeling. Icarus, 321, 486-507. doi: 10.1016/j.icarus.2018.11.025
* Morrissey, L. S., Pratt, D., Farrell, W. M., Tucker, O. J., Nakhla, S., & Killen,
  R. M. (2022). Simulating the diffusion of hydrogen in amorphous silicates: a 'jumping'
  migration process and its implications for solar wind implanted lunar volatiles.
  Icarus, 379, 114979. doi: 10.1016/j.icarus.2022.114979
"""
diffusion_time(h::S, D0::S, T::S, E::S) where {S<:AbstractFloat} = arrhenius(h^2 / D0, T, -E)
diffusion_time(h::Real, D0::Real, T::Real, E::Real) = diffusion_time(promote(h, D0, T, E)...)
function diffusion_time(h::Integer, D0::Integer, T::Integer, E::Integer)
    return diffusion_time(promote(h, D0, T, E, 1.0)[1:4]...)
end
diffusion_time(S::Type{<:AbstractFloat}, args...) = S(diffusion_time(args...))


############################################################################################
#::. EXPORTS
############################################################################################
export arrhenius, diffusion_time
