""" 
A spacecraft is orbiting the Earth as shown in figure below. 
As shown in the figure, at this particular location in the orbit, 
the Earth-pointing and Sun-pointing vectors are given in the ECI frame as:

    ne = [-1, 0, 0]; ns = [0, 1, 0]

Also the spacecraft attitude is obtained by rotating 45 degrees about vector zG. 
The ECI coordinate have sub-indices G and the body coordinates sub-indices b.
"""

#%% Import libraries
import numpy as np
# from math import sin, cos, pi
from sys import path
path.append(
    "c:\\Users\\diego\\Dropbox\\Academic\\MEng Space Systems\\3. DOCA\\ADCS functions")
import ADCS_Functions as adcs
import ADCS_Functions_sym as adcs_sym

#%% Data
# Define vectors
ne_eci = np.array([-1, 0, 0]).T
ns_eci = np.array([0, 1, 0]).T

#%% 1) Determine rotation matrix CbG from body coordinates to ECI coordinates.

# Define rotation matrix
#CbG = np.array([[cos(pi/4), sin(pi/4), 0],
#                [-sin(pi/4), cos(pi/4), 0],
#                [0, 0, 1]])
#CbG = np.array(adcs_sym.C3().subs("theta_3", np.pi/4)).astype(np.float64)
#print('rotation matrix CbG = ', CbG)

# or using the library
CbG = adcs_sym.DCM('num', 3, Eul_ang=np.array([np.pi/4]))
print('rotation matrix CbG = ', CbG)
#%% 2) Determine the coordinates of Earth and Sun vectors ne and ns respectively, in the spacecraft body frame.

# Calculate body coordinates
ne_body = np.dot(CbG, ne_eci)
ns_body = np.dot(CbG, ns_eci)

print('ne_body = ', ne_body)
print('ns_body = ', ns_body)


#%%  3) Using the TRIAD method, construct the body frame triad with vectors t1b, t2b and t3b with the spacecraft coordinates of the unit vectors. Construct the reference frame triad with vectors t1i, t2i and t3i with the ECI coordinates of the unit vectors

rot_triad_b, rot_triad_i, CbG_triad = adcs.triad(ne_body, ns_body, ne_eci, ns_eci)


t1b = rot_triad_b[:,0]
t2b = rot_triad_b[:,1]
t3b = rot_triad_b[:,2]

t1i = rot_triad_i[:,0]
t2i = rot_triad_i[:,1]
t3i = rot_triad_i[:,2]

print('t1b = ', t1b), print('t2b = ', t2b), print('t3b = ', t3b)
print('t1i = ', t1i), print('t2i = ', t2i), print('t3i = ', t3i)
#%%  4) Obtain the rotation matrices of [t1b, t2b, t3b] to [t1i, t2i, t3i]

print('triad_b = ', rot_triad_b)
print('triad_i = ', rot_triad_i)

#%%  5) Using the solution to part 4), compute rotation matric CbG using TRIAD method. Compare this with the result in part 1).
print('\n'); print('CbG = ', CbG_triad)
print('Does the rot matrix coincide with the one from part 1?',
np.allclose(CbG, CbG_triad))

#%%  6) Using the measured vectors obtained in part 2), compute CbG using the q-method and QUEST method. 
# Verify that you obtain the same result as in part 4). Note that it should be exactly the same, 
# since no measurement noise has been added
 
######## q-METHOD ########
b = np.array([ne_body, ns_body])
RF = np.array([ne_eci, ns_eci])

B, k22, K11, k12, K, max_Eigenvalue, max_Eigenvector, C = adcs.q_method(b, RF)

print('\n'); print(C)
print('The RM (q-METHOD) is equal to that of the TRIAD method:', np.allclose(C, CbG))

######## QUEST-METHOD ########
S, K, p, q, C = adcs.QUEST(b, RF)

print('\n'); print(C)
print('The RM (QUEST-METHOD) is equal to that of the TRIAD method:', np.allclose(C, CbG))
