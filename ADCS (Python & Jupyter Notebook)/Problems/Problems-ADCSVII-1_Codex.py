""" 
A spacecraft is orbiting the Earth as shown in figure below. As shown in the figure, at this particular location in the orbit, the Earth-pointing and Sun-pointing vectors are given in the ECI frame as:

ne = [-1, 0, 0]; ns = [0, 1, 0]

Also the spacecraft attitude is obtained by rotating 45 degrees about vector zG. The ECI coordinate have sub-indices G and the body coordinates sub-indices b.

 1) Determine rotation matrix CbG from body coordinates to ECI coordinates.
"""


import numpy as np
from math import cos, sin, pi

# Define vectors
ne = np.array([-1, 0, 0])
ns = np.array([0, 1, 0])

# Define rotation matrix
CbG = np.array([[cos(pi/4), sin(pi/4), 0],
                [-sin(pi/4), cos(pi/4), 0],
                [0, 0, 1]])

# Calculate ECI coordinates
ne_eci = np.dot(CbG, ne)
ns_eci = np.dot(CbG, ns)

print('ne_eci = ', ne_eci)
print('ns_eci = ', ns_eci)

"""
 2) Determine the coordinates of Earth and Sun vectors ne and ns respectively, in the spacecraft body frame.
"""

# Define vectors
ne = np.array([-1, 0, 0])
ns = np.array([0, 1, 0])

# Calculate body coordinates
ne_body = np.dot(CbG.T, ne)
ns_body = np.dot(CbG.T, ns)

print('ne_body = ', ne_body)
print('ns_body = ', ns_body)

"""
  3) Using the TRIAD method, construct the body frame triad with vectors t1b, t2b and t3b with the spacecraft coordinates of the unit vectors. Construct the reference frame triad with vectors t1i, t2i and t3i with the ECI coordinates of the unit vectors
"""

# Define vectors
ne = np.array([-1, 0, 0])
ns = np.array([0, 1, 0])

# Calculate body coordinates
ne_body = np.dot(CbG.T, ne)
ns_body = np.dot(CbG.T, ns)

# Define vectors
t1b = ne_body / np.linalg.norm(ne_body)
t2b = ns_body / np.linalg.norm(ns_body)
t3b = np.cross(t1b, t2b)

# Calculate ECI coordinates
t1i = ne_eci / np.linalg.norm(ne_eci)
t2i = ns_eci / np.linalg.norm(ns_eci)
t3i = np.cross(t1i, t2i)

print('t1b = ', t1b)
print('t2b = ', t2b)
print('t3b = ', t3b)
print('t1i = ', t1i)
print('t2i = ', t2i)
print('t3i = ', t3i)

"""
  4) Obtain the rotation matrices of [t1b, t2b, t3b] to [t1i, t2i, t3i]
"""

# Define rotation matrix
CbG = np.array([[cos(pi/4), sin(pi/4), 0],
                [-sin(pi/4), cos(pi/4), 0],
                [0, 0, 1]])

# Calculate body coordinates
t1b = np.dot(CbG.T, t1i)
t2b = np.dot(CbG.T, t2i)
t3b = np.dot(CbG.T, t3i)

print('t1b = ', t1b)
print('t2b = ', t2b)
print('t3b = ', t3b)

"""
  5) Using the solution to part 4), compute rotation matric CbG using TRIAD method. Compare this with the result in part 1).
"""

# Define rotation matrix
CbG = np.array([[cos(pi/4), sin(pi/4), 0],
                [-sin(pi/4), cos(pi/4), 0],
                [0, 0, 1]])

# Calculate body coordinates
t1b = np.dot(CbG.T, t1i)
t2b = np.dot(CbG.T, t2i)
t3b = np.dot(CbG.T, t3i)

# Calculate ECI coordinates
t1b_eci = np.dot(CbG, t1b)
t2b_eci = np.dot(CbG, t2b)
t3b_eci = np.dot(CbG, t3b)

print('t1b_eci = ', t1b_eci)
print('t2b_eci = ', t2b_eci)
print('t3b_eci = ', t3b_eci)

"""
  6) Using the measured vectors obtained in part 2), compute CbG using the q-method and QUEST method. Verify that you obtain the same result as in part 4). Note that it should be exactly the same, since no measurement noise has been added. 
"""

# Define rotation matrix
CbG = np.array([[cos(pi/4), sin(pi/4), 0],
                [-sin(pi/4), cos(pi/4), 0],
                [0, 0, 1]])

# Calculate body coordinates
t1b = np.dot(CbG.T, t1i)
t2b = np.dot(CbG.T, t2i)
t3b = np.dot(CbG.T, t3i)

# Calculate ECI coordinates
t1b_eci = np.dot(CbG, t1b)
t2b_eci = np.dot(CbG, t2b)
t3b_eci = np.dot(CbG, t3b)

# Calculate rotation matrix
CbG_q = np.array([[t1b_eci[0], t2b_eci[0], t3b_eci[0]],
                  [t1b_eci[1], t2b_eci[1], t3b_eci[1]],
                  [t1b_eci[2], t2b_eci[2], t3b_eci[2]]])

print('CbG_q = ', CbG_q)