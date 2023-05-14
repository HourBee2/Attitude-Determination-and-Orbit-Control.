from astropy import units as u
from astropy import constants

# Gravitational constant
# G = 6.67e-11
G = constants.G

# Earth radius
r_e = 6.371e6 * u.m

# Earth mass
# M = 5.9722e24 * u.kg
M = constants.M_earth

# Orbit
r = 400e3 * u.m + constants.R_earth

# Pert
J2 = 1.083e-3

# Radial acc for satellite at 400km due to J2 term
coef = ((3 * G * M) / r**4) * r_e**2 * J2
print(coef)


print(constants.R_earth.value)
