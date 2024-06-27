############################################################################################
#::. CONSTANTS
############################################################################################
"""
    CERES_DAY::Float64 = 32666.40

Length of a (synodic) day on Ceres in (s).

https://nssdc.gsfc.nasa.gov/planetary/factsheet/asteroidfact.html
"""
const CERES_DAY = 9.074 * 3600 # (s)

"""
    CERES_MASS::Float64 = 9.393e20

Mass of Ceres in (kg).

https://nssdc.gsfc.nasa.gov/planetary/factsheet/asteroidfact.html
"""
const CERES_MASS = 9.393e20 # (kg)


"""
    CERES_RADIUS::Float64 = 476e3

Mean radius of Ceres in (m).

* https://science.nasa.gov/dwarf-planets/ceres/facts/
* https://nssdc.gsfc.nasa.gov/planetary/factsheet/asteroidfact.html (965 x 961 x 891 km)
"""
const CERES_RADIUS = 476e3 # (m)


############################################################################################
#::. EXPORTS
############################################################################################
export CERES_DAY, CERES_MASS, CERES_RADIUS
