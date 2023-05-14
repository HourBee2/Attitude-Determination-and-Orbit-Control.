"""
Consider a satellite at perigee height hp = 800 km altitude above the Earth’s 
surface, and velocity (perpendicular to the radius) of v = 8 km/s.
(a) Compute the semi-major axis a.
(b)Compute the eccentricity e.
(c) Compute the maximum radius of this orbit ra.
(d)Compute the maximum altitude hmax
"""
import math 

# Gravitational constant
G = 6.67e-11

# Earth radius
r_e = 6.371e6

# Earth mass
M = 5.9722e24

# Orbit
hp = 800e3 + r_e
v = 8e3

mu = G * M
a = 1/2 * (-mu / ((v**2)/2-(mu/hp)))

e = 1 - (hp/a)

ra = a * (1 + e)

hmax = ra - r_e

print("a:", round(a/1000,2), "km")
print("e:", round(e,2))
print("ra:", round(ra/1000,2), "km")
print("hmax:", round(hmax/1000,2), "km")