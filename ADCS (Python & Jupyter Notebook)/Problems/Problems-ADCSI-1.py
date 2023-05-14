"""
Make a program in Python to determine the radial acceleration due to gravitational force for the Earth at:
(a) Sea surface.
(b) For an Earth observation satellite at 800 km altitude.
(c) For a GPS satellite at 20200 km altitude.
(d) For a geostationary satellite at 35800 km altitude
"""

# Gravitational constant
G = 6.67430e-11

# Earth radius
r_e = 6.378e6

# Earth mass
m_e = 5.9722e24

# Earth surface
r_s = r_e

# Earth observation satellite
r_os = 6.8e6

# GPS satellite
r_gs = 2.07e7

# Geostationary satellite
r_geo = 3.58e7

# Calculate the radial acceleration due to gravitational force
a_s = G * m_e / r_s ** 2
a_os = G * m_e / r_os ** 2
a_gs = G * m_e / r_gs ** 2
a_geo = G * m_e / r_geo ** 2

print("Radial acceleration due to gravitational force at sea surface:", a_s, "m/s^2")
print("Radial acceleration due to gravitational force at 800 km altitude:", a_os, "m/s^2")
print("Radial acceleration due to gravitational force at 20200 km altitude:", a_gs, "m/s^2")
print("Radial acceleration due to gravitational force at 35800 km altitude:", a_geo, "m/s^2")