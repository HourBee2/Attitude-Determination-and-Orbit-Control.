""" 
The initial yaw, pitch, and roll angles (of a 3-2-1 Euler angle sequence) 
of a vehicle are (theta1, theta2, theta3) = (80, 30, 40) at time t=0.
Assume that the angular velocity vector of the vehicle is given in body
frame components as: (w1, w2, w3) = (sin(0.1*t), 0.01, cos(0.1*t)) * 5 deg/s.
The time t is given in seconds.

A program to numerically integrate the yaw, pitch and roll angles over a 
simulation time of 1 minute. Integrate using radians but show the results in 
degrees.
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

#%% Data
theta = np.deg2rad([80, 30, 40])
time = np.linspace(0, 60, 61)

#%% Solve kinematic differential equation in matrix form
#time_sol, dot_angles = adcs.solve_KDE(theta, time_range=[0, 60], time_array=time, solver='E321')

time_sol, dot_angles = adcs.solve_KDE(
                        theta, time_range=[0, 60], time_array=time, solver='E321')
dot_angles = np.rad2deg(dot_angles)

#%% Plot the results
adcs.plot_euler(time_sol, dot_angles)
#
