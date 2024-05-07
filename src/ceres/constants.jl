############################################################################################
#::. CONSTANTS
############################################################################################
"""
Length of a day on Ceres in seconds.

https://nssdc.gsfc.nasa.gov/planetary/factsheet/asteroidfact.html
"""
const CERES_DAY = 9.074 * 3600 # [s]

"""
Mass of Ceres in kg.

https://nssdc.gsfc.nasa.gov/planetary/factsheet/asteroidfact.html
"""
const CERES_MASS = 9.393e20 # [kg]


"""
Mean radius of Ceres in m.

* https://science.nasa.gov/dwarf-planets/ceres/facts/
* https://nssdc.gsfc.nasa.gov/planetary/factsheet/asteroidfact.html (965 x 961 x 891 km)
"""
const CERES_RADIUS = 476e3    # [m]


############################################################################################
#::. EXPORTS
############################################################################################
export CERES_MASS, CERES_RADIUS
