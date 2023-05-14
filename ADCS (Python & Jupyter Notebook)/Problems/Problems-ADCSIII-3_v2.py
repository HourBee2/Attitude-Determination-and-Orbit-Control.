# Adapted from: https://github.com/nbelakovski/CollinearLagrangePoints/blob/master/lagrange.py

#%% Import Libraries
from math import pi
import matplotlib.pyplot as plt

# import numpy as np

# from astropy import units as u
# from astropy import constants

#%% Define values
G = 6.67408e-11 # m^3 / (kg * s^2)
m1 = 1.988409870698051e+30  # float(constants.M_sun / u.kg) # kg, mass of sun
m2 = 1.8981245973360505e+27 # float(constants.M_jup / u.kg)# kg, mass of jup
l = 778547200e3 #m, distance from sun to jup
x1 = m2 * l / (m1 + m2)
x2 = m1 * l / (m1 + m2)
velocity = (G * (m1 + m2)/ l) ** 0.5
period = 2 * pi * l / velocity
theta_dot = 2*pi / period

#%% Define function
def x_eqn(xs):
    return -G*m1/((xs + x1) * abs(xs +x1)) - G * m2 / ((xs - x2)*abs(xs-x2)) + \
            theta_dot**2 * xs

from numpy import arange

xvals = arange(-2*l, 2*l, 1e7)

yvals = x_eqn(xvals)

plt.plot(xvals, yvals)
plt.grid()
plt.ylim(-0.2, 0.2)
plt.show()

#%% Find roots

from scipy.optimize import root_scalar

L3 = root_scalar(x_eqn, bracket=[-4e11, -2e11])
print(L3.root+x1)

L1 = root_scalar(x_eqn, bracket=[2e11, 3.5e11])
print(L1.root+x1)

L2 = root_scalar(x_eqn, bracket=[3.9e11, 5e11])
print(L2.root+x1)


# unstable_Lagrange = np.array([L1.root, L2.root, L3.root])

# with plt.style.context('ggplot'):
#     plt.plot( unstable_Lagrange/1e3, np.zeros(len(unstable_Lagrange)), 'x', color='r') 