"""
Consider a star with mass M, radius R and uniform density \rho_0=\frac{3}{4\pi \:}\frac{M}{R^3}. Use the equation of hydrostatic equilibrium ( \frac{dp}{dr}=-g\left(r\right)\rho \left(r\right) ) and show that the pressure as a function of radius in the star is given by p\left(r\right)=\frac{2\pi }{3}G\rho_0^2\left(R^2-r^2\right).
"""

#%%
import numpy as np
import matplotlib.pyplot as plt

#%%
# Define constants
G = 6.67e-11
M = 1.989e30
R = 6.96e8
rho_0 = 3*M/(4*np.pi*R**3)

#%%
# Define function for pressure
def pressure(r):
    return 2*np.pi/3*G*rho_0**2*(R**2-r**2)

#%%
# Plot pressure as a function of radius
r = np.linspace(0,R,100)
plt.plot(r,pressure(r))
plt.xlabel('Radius (m)')
plt.ylabel('Pressure (Pa)')
plt.show()

#%%
from sympy import *

M, r, rho, g, p = symbols("M r rho g p")

# Define the function
f = (rho*g*r**2)/(4*pi*M)

# Find the differential of the function
df = diff(f, r)

# Solve for p
p = solve(df, p)

print(p)