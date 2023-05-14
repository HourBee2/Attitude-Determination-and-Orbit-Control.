# A program to find the direction cosine Euler rotation matrix, the principle Euler eigenaxis, the principal Euler rotation eigenaxis,
# and the Euler parameters or quaternions

#%% Import libraries
import numpy as np
from math import sin, cos, acos

import sys
sys.path.append(
   "c:\\Users\\diego\\Dropbox\\Academic\\MEng Space Systems\\3. DOCA\\ADCS functions")
import ADCS_Functions as adcs
import ADCS_Functions_sym as adcs_sym

# The orientation of an object is given in terms of the 3-2-1 Euler angles (-15, 25, 10 ).
Eul_ang = np.array([-15, 25, 10])
Eul_ang = np.radians(Eul_ang)

#%% Write the direction cosine Euler rotation matrix C21.

# Use the direction cosine matrix for the  3-2-1 Euler angles
# Note! can also use: THIS WORKS FOR ALL TYPES OF SEQUENCES
C21 = adcs_sym.DCM('num', 3, 2, 1, Eul_ang=Eul_ang, invorder=True)
print('DCM 321 is:', C21); print('\n')


#%% Find the principle Euler eigenaxis rotation angle phi
phi = adcs.Eigenaxis_rot(C21) 
print('The principle Euler eigenaxis rotation angle is: ', np.rad2deg(phi))

#%% Find the principle Euler eigenaxis e
e = adcs.Eigenaxis_e(C21)
print(e)

#%% Verify that C21*e=e
print(np.allclose(C21@e, e))

#%% Find the corresponding Euler parameters = Quaternions
q = adcs.DCM_to_Quaternion(C21)
print(q)
## Check that it is teh same using eigenaxis and eigenaxis rotation angle
qq = adcs.Eigenaxis_to_Quaternion(e, phi)
print(qq)
#%% Is the quaternion q a unit quaternion?
print(np.allclose(q, q/np.linalg.norm(q)))
