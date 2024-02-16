#::. FUNCTIONS
"""
    [1] projection_CHAMBERLAIN1963(n1::Real, r1::Real, r2::Real, T::Real, m::Real; M::Real=LUNAR_MASS)

Calculates the number density of particles of mass `m` at radial distance `r2` given the 
number density `n1` at radial distance `r1`, assuming constant temperature `T`.

Assumes hydrostatic equilibrium and perfect-gas law with isotropic gas pressure, resulting
in the generalized form of the isothermal barometric law. The calculation is based on the 
equation (14) in Chamberlain (1963).


**Rerences**

- Chamberlain, 1963, "Planetary coronae and atmospheric evaporation"
- Cook et al. 2013, "New upper limits on numerous atmospheric species in the native lunar
  atmosphere"
"""
function projection_CHAMBERLAIN1963(n1::Real, r1::Real, r2::Real, T::Real, m::Real; M::Real=LUNAR_MASS)
    pot_energy_1 = GRAVITATIONAL_CONSTANT * M * m / (BOLTZMANN_CONSTANT * T * r1)
    pot_energy_2 = GRAVITATIONAL_CONSTANT * M * m / (BOLTZMANN_CONSTANT * T * r2)
    
    # pre-calculations
    L = pot_energy_1^2 / (pot_energy_1 + pot_energy_2)
    g =  gamma(1.5)
    gi = gamma_inc(1.5, pot_energy_1)[1] * g
    giL = gamma_inc(1.5, pot_energy_1 - L)[1] * g
    f = sqrt(pot_energy_2^2 - pot_energy_1^2) / pot_energy_2

    # partition functions (ballistic, satellite, escape)
    par_bal = 2/sqrt(pi) * (gi - f*exp(-L)*giL)
    par_sat = 2/sqrt(pi) * (f*exp(-L)*giL)
    par_esc = 1/sqrt(pi) * (g - gi - f*exp(-L)*(g - giL))
    par = par_bal + par_sat + par_esc

    return n1 * exp(-(pot_energy_1 - pot_energy_2)) * par
end



#::. EXPORTS
export projection_CHAMBERLAIN1963