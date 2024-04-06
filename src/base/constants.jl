############################################################################################
#::. CONSTANTS
############################################################################################
"""
    AVOGADRO_CONSTANT::Float64 = 6.02214086e23

Avogadro constant in (mol-1).
Value taken from the [NIST CODATA database](https://physics.nist.gov/cuu/Constants/).
"""
const AVOGADRO_CONSTANT = 6.02214086e23 # (mol-1)


"""
    BOLTZMANN_CONSTANT::Float64 = 1.380649e-23

Boltzmann constant in (J K-1).
Value taken from the [NIST CODATA database](https://physics.nist.gov/cuu/Constants/).
"""
const BOLTZMANN_CONSTANT = 1.380649e-23 # (J K-1)


"""
    ELEMENTARY_CHARGE::Float64 = 1.602176634e-19

Elementary charge in (C).
Value taken from the [NIST CODATA database](https://physics.nist.gov/cuu/Constants/).
"""
const ELEMENTARY_CHARGE = 1.602176634e-19 # (C)


"""
    GRAVITATIONAL_CONSTANT::Float64 = 6.67430e-11

Gravitational constant in (m3 kg-1 s-2).
Value taken from the [NIST CODATA database](https://physics.nist.gov/cuu/Constants/).

**Note**: The constant has a documented uncertainty of `0.000 15e-11` in (m3 kg-1 s-2).
"""
const GRAVITATIONAL_CONSTANT = 6.67430e-11 # (m3 kg-1 s-2)


"""
    PLANCK_CONSTANT::Float64 = 6.62607015e-34

Planck constant in (m2 kg s-1).
Value taken from the [NIST CODATA database](https://physics.nist.gov/cuu/Constants/).
"""
const PLANCK_CONSTANT  = 6.62607015e-34 # (m2 kg s-1)


"""
    STEFAN_BOLTZMANN_CONSTANT::Float64 = 5.670374419e-8

Stefan-Boltzmann constant in (W m-2 K-4).
Value taken from the [NIST CODATA database](https://physics.nist.gov/cuu/Constants/).
"""
const STEFAN_BOLTZMANN_CONSTANT = 5.670374419e-8 # (W m-2 K-4)


"""
    UNIVERSAL_GAS_CONSTANT::Float64 = 8.314462618

Universal gas constant in (J K-1 mol-1).
Value taken from the [NIST CODATA database](https://physics.nist.gov/cuu/Constants/).
"""
const UNIVERSAL_GAS_CONSTANT = 8.314462618 # (J K-1 mol-1)


"""
    AMU_H::Float64 = 1.007975

Average atomic mass of hydrogen in (u). 
Value taken from [Prohaska et al.](https://doi.org/10.1515/pac-2019-0603).
The value used is the average of the reported interval `[1.00784, 1.00811]` in (u).
"""
const AMU_H = (1.00784 +  1.00811) / 2 # (u)


"""
    AMU_He::Float64 = 4.002602

Average atomic mass of helium in (u).
Value taken from [Prohaska et al.](https://doi.org/10.1515/pac-2019-0603).
"""
const AMU_He = 4.002602 # (u)


"""
    AMU_O::Float64 = 15.9994

Average atomic mass of oxygen in (u).
Value taken from [Prohaska et al.](https://doi.org/10.1515/pac-2019-0603).
The value used is the average of the reported interval `[15.99903, 15.99977]` in (u).
"""
const AMU_O = (15.99903 + 15.99977) / 2 # (u)


"""
    AMU_Ne::Float64 = 20.1797

Average atomic mass of neon in (u).
Value taken from [Prohaska et al.](https://doi.org/10.1515/pac-2019-0603).
"""
const AMU_Ne = 20.1797 # (u)


"""
    AMU_Ar::Float64 = 39.8775

Average atomic mass of argon in (u).
Value taken from [Prohaska et al.](https://doi.org/10.1515/pac-2019-0603).
The value used is the average of the reported interval `[39.792, 39.963]` in (u).
"""
const AMU_Ar = (39.792 + 39.963) / 2 # (u)


"""
    AMU_H2::Float64

Average molecular mass of molecular hydrogen in (u).
The value is an alias for `2 * AMU_H`.
"""
const AMU_H2 = 2 * AMU_H # (u)


"""
    AMU_OH::Float64

Average molecular mass of hydroxyl radical in (u). 
The value is an alias for `AMU_H + AMU_O`.
"""
const AMU_OH = AMU_H + AMU_O # (u)


"""
    AMU_H2O::Float64

Average molecular mass of water in (u).
The value is an alias for `2 * AMU_H + AMU_O`.
"""
const AMU_H2O = 2 * AMU_H + AMU_O # (u)


############################################################################################
#::. EXPORTS
############################################################################################
export 
    AVOGADRO_CONSTANT, BOLTZMANN_CONSTANT, ELEMENTARY_CHARGE, GRAVITATIONAL_CONSTANT, 
    PLANCK_CONSTANT, STEFAN_BOLTZMANN_CONSTANT, UNIVERSAL_GAS_CONSTANT, 

    AMU_H, AMU_He, AMU_O, AMU_Ne, AMU_Ar, AMU_H2, AMU_OH, AMU_H2O
