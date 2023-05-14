#%% Import libraries
import numpy as np
# import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants

#%% Functions
def F_SRP(phi_rad, A, q):
    Fsrp = ((phi_rad / constants.c) * (1+q) * A).to(u.N)
    return Fsrp

def a_SRP(phi_rad, A, q, m):
    a = ((phi_rad / constants.c) * (1+q) * (A/m)).to(u.m/u.s**2)
    return a
#%%  
phi_rad = 1367 * u.W /u.m**2
A = 1e4 * u.m**2
q=1
m = 1000 * u.kg

a_rad= a_SRP(phi_rad, A, q, m)
print(a_rad)
#%%
d = 4e8 * u.m # distance earth moon

time = np.sqrt(2*d / a_rad)
time_h = time.to(u.hr)
time_d = time.to(u.day)
print(time_d)