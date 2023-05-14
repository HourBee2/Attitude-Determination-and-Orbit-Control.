# Import libraries
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants

# Constants (with units)
# G = constants.G # Gravitational constant
# Re = constants.R_earth # Earth radius
# M = constants.M_earth # Earth mass
# J2 = 1.083e-3 # Pert J2 term
# mu = constants.GM_earth

# Variables (with units)
# NRR = 2*np.pi/(365.25*24*60*60) *u.rad*u.s**(-1)# Nodal regression rate (SS orbit)

# Constants
G = 6.67e-11 # Gravitational constant
Re = 6.371e6 # Earth radius
M = 5.9722e24 # Earth mass
J2 = 1.083e-3 # Pert J2 term
mu = G*M

# Variables
NRR = 2*np.pi/(365.25*24*60*60) # Nodal regression rate (SS orbit)
#%% Compute the eccentricity and radii of perigee for a geocentric 
#   sun-synchronous frozen orbits with semi-major axes a = 10000 km, 
#   a = 15000 km and a = 20000 km. Conclude that there is a range of 
#   semi-major axes for which a geocentric sun-synchronous frozen orbit 
#   is possible.

# a = [10000*1e3, 15000*1e3, 20000*1e3] * u.m
# a = [10000*1e3, 15000*1e3, 20000*1e3] 

a = np.linspace(0,50000e3,100000) 

# e, rp, hp = np.zeros((len(a))), np.zeros((len(a))) * u.m, np.zeros((len(a))) * u.m
e, rp, hp = np.zeros((len(a))), np.zeros((len(a))), np.zeros((len(a))) 

for ii in range(0,len(a)):
    # e[ii] = np.sqrt( 1 - np.sqrt( ( (3*J2*(Re**2))/(2*NRR) ) * np.sqrt( mu / (5*a[ii]**7)) ) )
    e[ii] = (1-( ( (3*J2*(Re**2))/(2*NRR) ) * np.sqrt( mu / (5*a[ii]**7) ) )**(1/2))**(1/2)
    rp[ii] = a[ii] * (1-e[ii])
    hp[ii] = rp[ii] - Re

#%% Sketch a plot of semi-major axis versus radius of perigee for a 
#   geocentric sun-synchronous frozen orbit 
with plt.style.context('ggplot'):
    plt.plot(rp, a) 
    plt.axvline(x=Re, color = 'k', linestyle='--')
    plt.xlabel("Radius of perigee, $r_p$ [m]")
    plt.ylabel("Semi-major axes, $a$ [m]")
    plt.title("Semi-major axis versus Radius of perigee for a geocentric SS frozen orbit")
    plt.tight_layout()
    plt.show()

with plt.style.context('ggplot'):
    plt.plot(a, e)
    plt.axvline(x=9809e3, color = 'k', linestyle='--')
    plt.ylabel("Eccentricity, $a$ [m]")
    plt.xlabel("Semi-major axes, $a$ [m]")
    plt.title("Semi-major axis versus Radius of perigee for a geocentric SS frozen orbit")
    plt.tight_layout()
    plt.show()

#%% Find the minimum value of a for which a geocentric sun-synchronous 
#   frozen orbit is possible.

# with a being a 500 term array -> e=nan until 98th term --> a[98]*1e-3=9819.64 km

# with a being a 1000 term array -> e=nan until 196th term --> a[196]*1e-3=9809.81 km# a = np.linspace(0,50000e3,1000) 

# # with a being a 100 000 term array -> e=nan until 19619th term --> a[19619]*1e-3=9809.59 km

def ecc_f(a):
    """ 
    Eccentricity for a given geocentric SS frozen orbit.
    Input:  a - semi major axis (m)
    Output: e - eccentricity
    """
    e = (1-( ( (3*J2*(Re**2))/(2*NRR) ) * np.sqrt( mu / (5*a**7) ) )**(1/2))**(1/2)
    return e

from scipy.optimize import fminbound # , minimize

# minimize(ecc_f(a), 0)
MIN = fminbound(ecc_f, 8e3, 10e3) # np.array([1,0])   # MIN
print(MIN)

#%% Using an iterative procedure, determine the maximum value of a for 
#   which a geocentric frozen orbit is possible, given that the orbit should 
#   stay at least 200 km above the Earth

def max_ecc_f(a):
    """ 
    Eccentricity for a given geocentric SS frozen orbit.
    Input:  a - semi major axis (m)
    Output: -e - eccentricity
    """
    e_max = -(1-( ( (3*J2*(Re**2))/(2*NRR) ) * np.sqrt( mu / (5*a**7) ) )**(1/2))**(1/2)
    return e_max

MAX = fminbound(max_ecc_f, 0, 5) # np.array([1,0])  # MAX # MAX = fmin(max_ecc_f, 1)
print(MAX)