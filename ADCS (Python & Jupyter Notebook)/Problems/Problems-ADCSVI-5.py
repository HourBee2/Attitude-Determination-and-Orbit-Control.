# A program to find the relative orientation of 2 spacecrafts

#%% Import libraries
import numpy as np
from math import sin, cos
from sys import path
path.append(
    "c:\\Users\\diego\\Dropbox\\Academic\\MEng Space Systems\\3. DOCA\\ADCS functions")
import ADCS_Functions as adcs
import ADCS_Functions_sym as adcs_sym

#%% Data and functions
# Let the orientations of two spacecraft A and B relative to an inertial frame I 
# be given through the 3-2-1 Euler angles rotation sequences:
theta_A = np.deg2rad(np.array([60, -45, 30]).T)
theta_B = np.deg2rad(np.array([-15, 25, 10]).T)

# Direction cosine matrix for 3-2-1 Euler angles
#%% Orientation matrices CAI and CBI
#C_AI = adcs.DCM_321(theta_A)
#C_BI = adcs.DCM_321(theta_B)

C_AI = adcs_sym.DCM('num', 3, 2, 1, Eul_ang=theta_A, invorder=True)
C_BI = adcs_sym.DCM('num', 3, 2, 1, Eul_ang=theta_B, invorder=True)

#%% Direction cosine matrix CAB
C_AB = np.dot(C_AI, np.linalg.inv(C_BI))

#%% Euler angles from DCM
Eul_ang_AB = adcs.Eul_ang(C_AB)
print('Euler angles of the orientation of the spacecraft A relative to spacecraft B:', (Eul_ang_AB))
