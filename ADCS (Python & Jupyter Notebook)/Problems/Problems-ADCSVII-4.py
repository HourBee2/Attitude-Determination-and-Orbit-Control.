""" 
The initial yaw, pitch, and roll angles (of a 3-2-1 Euler angle sequence) 
of a vehicle are (theta1, theta2, theta3) = (80, 30, 40) at time t=0.
Assume that the angular velocity vector of the vehicle is given in body
frame components as: (w1, w2, w3) = (sin(0.1*t), 0.01, cos(0.1*t)) * 5 deg/s.
The time t is given in seconds.

A program to numerically integrate the quaternions over a 
simulation time of 1 minute.
"""

#%% Import libraries
import numpy as np
from math import sin, cos
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

from sys import path
path.append(
    "c:\\Users\\diego\\Dropbox\\Academic\\MEng Space Systems\\3. DOCA\\ADCS functions")
import ADCS_Functions as adcs
import ADCS_Functions_sym as adcs_sym

#%% Data
theta = np.deg2rad([80, 30, 40])
time = np.linspace(0, 60, 244)

#%% Solve ODE

## 1) Translate from Euler angles, theta, to quaternions, q
# Compute the corresponding rotation matrix from the 3-2-1 Euler angles
#C = adcs.DCM_321(theta)
C = adcs_sym.DCM('num', 3, 2, 1, Eul_ang=theta, invorder=True) # I use the flip function because in the solutions it considers that the given angles are in the opposite order
print("The corresponding rotation matrix is =", C)

# Find the corresponding Euler parameters = Quaternions
q = adcs.DCM_to_Quaternion(C)
print("The quaternion is = ", q)

## 2) Solve the ODE
time_sol, q_sol = adcs.solve_KDE(
    q, time_range=[0, 60], time_array=time, solver='q')

#%% Plot the results
adcs.plot_quaternion(time_sol, q_sol)

#%% Given the result of the numerical integration, plot the quaternion (or
# Euler parameters) constraint |q|= 1 and comment on this constraint

adcs.plot_quaternion_constraint(time_sol, q_sol)
