# A program to convert between reference frames using directional cosine matrices.

#%% Import libraries
import numpy as np
# from math import sin, cos, pi

#%% Data
# Assume three reference frames A, B and I. 
# Let the two reference frames A and B be defined relative to the inertial reference frame I by the orthonormal unit base vectors.

a1, a2, a3 = np.array([1/2, np.sqrt(3)/2, 0]), np.array([0, 0, 1]), np.array([np.sqrt(3)/2, -1/2, 0])
b1, b2, b3 = np.array([0, 1, 0]), np.array([1, 0, 0]), np.array([0, 0, -1])

# where the ai and bi (i = 1, 2, 3) vector components are written in the inertial frame I. 
# Note that the unit base vectors of the inertial frame are:

i1, i2, i3 = np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])

#%% (a) Check that the unit base vectors ai respectively bi (i = 1, 2, 3) build an orthonormal reference frame.

print(np.dot(a1, a2)==0)
print(np.dot(a2, a3)==0)
print(np.dot(a1, a3)==0)

print(round(np.dot(a1, a1))==1)
print(round(np.dot(a2, a2))==1)
print(round(np.dot(a3, a3))==1)

#%% (b) Find the directional cosine matrix Cab that describes the orientation of frame A relative to frame B
Cab = np.array([[np.dot(a1, b1), np.dot(a1, b2), np.dot(a1, b3)],
                [np.dot(a2, b1), np.dot(a2, b2), np.dot(a2, b3)],
                [np.dot(a3, b1), np.dot(a3, b2), np.dot(a3, b3)]])
print(Cab)

#%% (c) Find the directional cosine matrix Cai that describes the orientation of frame A relative to frame I
Cai = np.array([[np.dot(a1, i1), np.dot(a1, i2), np.dot(a1, i3)],
                [np.dot(a2, i1), np.dot(a2, i2), np.dot(a2, i3)],
                [np.dot(a3, i1), np.dot(a3, i2), np.dot(a3, i3)]])

print(Cai)

#%% (d) Find the directional cosine matrix Cab that describes the orientation of frame B relative to frame I
Cbi = np.array([[np.dot(b1, i1), np.dot(b1, i2), np.dot(b1, i3)],
                [np.dot(b2, i1), np.dot(b2, i2), np.dot(b2, i3)],
                [np.dot(b3, i1), np.dot(b3, i2), np.dot(b3, i3)]])

print(Cbi)

#%% (e) Check if Cab = Cai*(Cbi).T holds
print(np.allclose(Cab, Cai@(Cbi.T)))

#%% (f) Check if Cab*(Cab).T = I
print(np.allclose(Cab@(Cab.T), np.eye(3)))

#%% (g) For given arbitrary matrix A and matrix B check if they do not commute (AB != BA)
A = Cai
B = (Cbi).T
print(np.allclose(A@B, B@A))

#%% (h) Is the following matrix C = [1/2, 0, 0; 0, 1, 0; 0, 0, 2] a rotation matrix?
C = np.matrix([[1/2, 0, 0],
              [0, 1, 0],
              [0, 0, 2]])
print(np.allclose(C@C.T, 1))