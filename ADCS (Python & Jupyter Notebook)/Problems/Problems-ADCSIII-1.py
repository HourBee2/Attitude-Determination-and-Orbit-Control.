#%% Import libraries
# import numpy as np
# import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants
import math

#%% Functions
def F_SRP(phi_rad, A, q):
    Fsrp = ((phi_rad / constants.c) * (1+q) * A).to(u.N)
    return Fsrp

#%%  (b)
phi_rad = 1000 * u.W /u.m**2
A = 1 * u.m**2
q=0

Fsrp1_abs = F_SRP(phi_rad, A, q)
print(Fsrp1_abs)
#%% (c)
q=1

Fsrp1_ref = F_SRP(phi_rad, A, q)
print(Fsrp1_ref)
#%% (d)
L = 3.9e26 * u.W
r_se = (150.02e6 * u.km).to(u.m)
phi_Re = 1367 *u.W/u.m**2 # Wo/(4*np.pi*(r_se**2))
r = 1.76e11 * u.m   
print(r)
