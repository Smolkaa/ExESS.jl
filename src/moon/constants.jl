############################################################################################
#::. CONSTANTS
############################################################################################
"""
Mass of the Moon in kg.

https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
"""
const LUNAR_MASS     = 7.24767309e22   # [kg]


"""
Radius of the Moon in m (volumetric mean).

https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
"""
const LUNAR_RADIUS   = 1737400e0       # [m]


const LUNAR_DAY      = 3600*24*295e-1  # [s]


const LUNAR_g0       = GRAVITATIONAL_CONSTANT * LUNAR_MASS / LUNAR_RADIUS^2 # [m s-2]



############################################################################################
#::. EXPORTS
############################################################################################
export LUNAR_DAY, LUNAR_MASS, LUNAR_RADIUS, LUNAR_g0