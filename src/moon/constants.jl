############################################################################################
#::. CONSTANTS
############################################################################################
"""
    LUNAR_MASS::Float64 = 7.346e22

Mass of the Moon in (kg).

https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
"""
const LUNAR_MASS = 7.346e22 # (kg)


"""
    LUNAR_RADIUS::Float64 = 1737400e0

Radius of the Moon in (m) (volumetric mean).

https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
"""
const LUNAR_RADIUS = 1737400e0 # (m)


"""
    LUNAR_ROTATIONA_PERIOD::Float64 = 2.36059488e6

Length of the revolution period of the Moon in (s).

https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
"""
const LUNAR_ROTATION_PERIOD = 3600 * 24 * 27.3217 # (s)


"""
    LUNAR_ORBITAL_PERIOD::Float64 = 2.36059488e6

Length of the orbital period of the Moon around Earth in (s).

https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
"""
const LUNAR_ORBITAL_PERIOD = LUNAR_ROTATION_PERIOD # (s)


"""
    LUNAR_DAY::Float64 = 2.55144300192e6

Length of a day on the Moon, with respect to the Sun (also: synodic period), in (s).

https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
"""
const LUNAR_DAY = 3600 * 24 * 29.5305903 # (s)


"""
    LUNAR_g0::Float64 = 1.6242654756205575

Gravitational surface acceleration on the Moon in (m s-2).

This constant is a precalculated value based on: `G * m_M / r_M^2`, with the gravitational
constant `G`, the lunar mass `m_M`, and the meand lunar radius `r_M`.
"""
const LUNAR_g0 = GRAVITATIONAL_CONSTANT * LUNAR_MASS / LUNAR_RADIUS^2 # (m s-2)



############################################################################################
#::. EXPORTS
############################################################################################
export
    LUNAR_DAY,
    LUNAR_g0,
    LUNAR_MASS,
    LUNAR_RADIUS,
    LUNAR_ROTATION_PERIOD,
    LUNAR_ORBITAL_PERIOD
