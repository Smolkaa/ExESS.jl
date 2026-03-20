############################################################################################
#::. CONSTANTS
############################################################################################
"""
    MERCURY_MASS::Float64 = 3.3010e23

Mass of Mercury in (kg).

https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html
"""
const MERCURY_MASS = 3.3010e23 # (kg)


"""
    MERCURY_ORBITAL_ECCENTRICITY::Float64 = 0.2056

Orbital eccentricity of Mercury.

https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html
"""
const MERCURY_ORBITAL_ECCENTRICITY = 0.2056 # (-)


"""
    MERCURY_ORBITAL_PERIOD::Float64 = 7.77333024e6

Orbital period of Mercury in (s). Equivalent to 87.9691 Earth days.
"""
const MERCURY_ORBITAL_PERIOD = 89.9691 * 24 * 3600 # (s)


"""
    MERCURY_PERIHELION::Float64 = 46e9

Perihelion distance of Mercury in (m).

https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html
"""
const MERCURY_PERIHELION = 46e9 # (m)


"""
    MERCURY_RADIUS::Float64 = 2439.7e3

Radius of Mercury in (m) (volumetric mean).

https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html
"""
const MERCURY_RADIUS = 2439.7e3 # (m)


"""
    MERCURY_ROTATION_PERIOD::Float64 = 5.0670144e6

Length of a rotation period of Mercury in (s). Equivalent to 58.646 Earth days.
"""
const MERCURY_ROTATION_PERIOD = 58.646 * 24 * 3600 # (s)


"""
    MERCURY_SEMI_MAJOR_AXIS::Float64 = 57.909e9

Semi-major axis of Mercury's orbit in (m).

https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html
"""
const MERCURY_SEMI_MAJOR_AXIS = 57.909e9 # (m)


"""
    MERCURY_SYNODIC_ROTATION_PERIOD::Float64 = 1.520136e7

Length of a synodic period (one day) of Mercury in (s). Equivalent to approx. 175.9417 
Earth days.

https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html
"""
const MERCURY_SYNODIC_ROTATION_PERIOD = 175.9416666666667 * 24 * 3600 # (s)


############################################################################################
#::. DERIVED CONSTANTS (for convenience)
############################################################################################
"""
    MERCURY_DAY::Float64 = 1.520136e7

Length of a day on Mercury in (s). Equivalent to approx. 175.9417 Earth days.

https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html
"""
const MERCURY_DAY = MERCURY_SYNODIC_ROTATION_PERIOD # (s)



############################################################################################
#::. EXPORTS
############################################################################################
export 
    MERCURY_DAY, 
    MERCURY_MASS, 
    MERCURY_ORBITAL_ECCENTRICITY,
    MERCURY_ORBITAL_PERIOD,
    MERCURY_PERIHELION,
    MERCURY_RADIUS,
    MERCURY_ROTATION_PERIOD,
    MERCURY_SEMI_MAJOR_AXIS,
    MERCURY_SYNODIC_ROTATION_PERIOD
