"""
A satellite is in a Keplerian orbit around the Earth with perigee height above 
the Earth of hp = 300 km and an apogee height above the Earth of 
ha = 10000 km. Please calculate 
(a) the semi-major axis a.
(b) the eccentricity e.
(c) the velocity at perigee vp.
(d) the velocity at apogee va.
(e) the orbital period T. 
"""

import math

# Gravitational constant
G = 6.67e-11

# Earth radius
r_e = 6.371e3

# Earth mass
M = 5.9722e24

# GIVEN
hp = 300.0 + r_e
ha = 10000.0 +r_e

# CALCULATIONS
a = (hp + ha) / 2
e = (ha - hp) / (ha + hp)
mu = G * M
vp = math.sqrt(2*(-mu/(2*a*1e3)+mu/(hp*1e3)))
va = math.sqrt(2*(-mu/(2*a*1e3)+mu/(ha*1e3)))
T = 2 * math.pi * math.sqrt((a*1e3)**3 / (G*M))

# DISPLAY RESULTS
print("a = {0:.2f} km".format(a))
print("e = {0:.2f}".format(e))
print("vp = {0:.2f} km/s".format(vp))
print("va = {0:.2f} km/s".format(va))
print("T = {0:.2f} s".format(T))
