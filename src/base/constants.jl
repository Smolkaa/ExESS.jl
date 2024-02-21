#::. CONSTANTS
const AVOGADRO_CONSTANT         = 6.02214086e23     # [mol-1]
const BOLTZMANN_CONSTANT        = 1.38065e-23       # [J K-1]
const ELEMENTARY_CHARGE         = 1.602176634e-19   # [C]
const GRAVITATIONAL_CONSTANT    = 6.67408e-11       # [m3 kg-1 s-2]
const PLANCK_CONSTANT           = 6.62607004e-34    # [m2 kg s-1]
const STEFAN_BOLTZMANN_CONSTANT = 5.67e-8           # [W m-2 K-4]
const UNIVERSAL_GAS_CONSTANT    = 8.3144626181532e0 # [J K-1 mol-1]


const AMU_H   =  (1.00784 +  1.00811) / 2 # [u] 
const AMU_He  =   4.002602                # [u] 
const AMU_O   = (15.99903 + 15.99977) / 2 # [u]
const AMU_Ne  =  20.1797                  # [u]
const AMU_Ar  = (39.792   + 39.963)   / 2 # [u]

const AMU_H2 = 2 * AMU_H            # [u] 
const AMU_OH = AMU_H + AMU_O        # [u]
const AMU_H2O = 2 * AMU_H + AMU_O   # [u]


#::. EXPORTS
export 
    AVOGADRO_CONSTANT, BOLTZMANN_CONSTANT, ELEMENTARY_CHARGE, GRAVITATIONAL_CONSTANT, 
    PLANCK_CONSTANT, STEFAN_BOLTZMANN_CONSTANT, UNIVERSAL_GAS_CONSTANT, 

    AMU_H, AMU_He, AMU_O, AMU_Ne, AMU_Ar, AMU_H2, AMU_OH, AMU_H2O
