# A program to find the Euler angles that relate the final attitude
# to the original attitude of a spacecraft that performs a 45 deg 
# single principle Euler eigenaxis rotation.

#%% Import libraries
import numpy as np
from math import sin, cos, pi
from sys import path
path.append(
    "c:\\Users\\diego\\Dropbox\\Academic\\MEng Space Systems\\3. DOCA\\ADCS functions")
import ADCS_Functions as adcs

#%% Data and functions
e = 1/np.sqrt(3) * np.array([1, 1, 1])
phi = pi/4
#%% Calculations
# Skew symmetric matrix
def skew(v):
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])

# The rotation matrix for a given principle Euler eigenaxis rotation is given by:
C = adcs.Eigenaxis_rotMAT(e, phi)
print(C)

#%% Euler angles from DCM
theta = adcs.Eul_ang(C)
print(theta)
