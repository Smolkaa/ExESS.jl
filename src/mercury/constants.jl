############################################################################################
#::. CONSTANTS
############################################################################################
"""
    MERCURY_MASS::Float64 = 3.3010e23

Mass of Mercury in kg.

https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html
"""
const MERCURY_MASS = 3.3010e23 # [kg]


"""
    MERCURY_RADIUS::Float64 = 2439.7e3

Radius of Mercury in m (volumetric mean).

https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html
"""
const MERCURY_RADIUS = 2439.7e3 # [m]


"""
    MERCURY_DAY::Float64 = 1.520136e7

Length of a day on Mercury in seconds.

https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html
"""
const MERCURY_DAY = 4222.6 * 3600 # [s]



############################################################################################
#::. EXPORTS
############################################################################################
export MERCURY_DAY, MERCURY_MASS, MERCURY_RADIUS
