############################################################################################
#::. CONSTANTS
############################################################################################
"""
Mass of the Moon in (kg).

https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
"""
const LUNAR_MASS = 7.346e22 # [kg]


"""
Radius of the Moon in (m) (volumetric mean).

https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
"""
const LUNAR_RADIUS = 1737400e0 # [m]


"""
Length of the revolution period of the Moon in (s).

https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
"""
const LUNAR_REVOLUTION_PERIOD = 3600 * 24 * 27.3217 # [s]


"""
Length of the orbital period of the Moon around Earth in (s).

https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
"""
const LUNAR_ORBITAL_PERIOD = LUNAR_REVOLUTION_PERIOD


"""
Length of a day on the Moon, with respect to the Sun (also: synodic period), in (s).
"""
const LUNAR_DAY =  3600 * 24 * 29.5305903 # [s]


"""
Gravitational surface acceleration on the Moon in (m s-2).

This constant is a precalculated value based on: `G * m_M / r_M^2`, with the gravitational
constant `G`, the lunar mass `m_M`, and the meand lunar radius `r_M`.
"""
const LUNAR_g0 = GRAVITATIONAL_CONSTANT * LUNAR_MASS / LUNAR_RADIUS^2 # [m s-2]



############################################################################################
#::. EXPORTS
############################################################################################
export
    LUNAR_DAY,
    LUNAR_g0,
    LUNAR_MASS,
    LUNAR_RADIUS,
    LUNAR_REVOLUTION_PERIOD,
    LUNAR_ORBITAL_PERIOD
