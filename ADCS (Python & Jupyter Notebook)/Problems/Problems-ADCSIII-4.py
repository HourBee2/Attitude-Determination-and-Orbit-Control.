#%% Import libraries
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants

#%% (a)
T = 86164 *u.s
omega = 2*np.pi / T

r = ((constants.G*constants.M_earth/omega**2)**(1/3)).to(u.km)
print(r)
#%% (b)
from sympy import init_session, Symbol, symbols, sin, Derivative
# from sympy import init_session
init_session(use_latex=True)

r = Symbol('r')
U = Symbol('U')
phi = symbols('Phi')
lamda = symbols('Lambda')

a_pertEW = - 1 / (r * sin(phi)) * Derivative(U, lamda)
print(a_pertEW)
print()

#%% (e)
# Orbit
r = 42164e3 * u.m

# J22 Pert term (Sectorial harmonics)
J22 = 1.816e-6

# East-West acc due to J2,2 for satellite in GEO
coef = 6*((constants.GM_earth) / r**4) * constants.R_earth**2 * J22 * np.cos(0)
print(coef)
print()

#%% (f)
# Function and data
def acc_EW(lamda):
    return coef * np.sin(2*(lamda-lamda22))

lamda = np.linspace(0, 2*np.pi, 100)
lamda22 = np.deg2rad(-14.9)

acc_ew  = acc_EW(lamda)

# Solve for the roots of the function
from scipy.optimize import root

for i in lamda:
    sol = []
    c = abs(int(round(i)))
    for j in range(0, c+1):
        y = root(acc_EW, j)
        if y.success and (round(y.x[0], 4) not in sol):
            sol.append(round(y.x[0], 3))
# Sort the roots, keep only the ones we interested in            
sol = list(set(sol)) # delete duplicates
for i in range(len(sol)):
    sol = [i for i in sol if i > 0] # Using list comprehension to delete sol<0 values
sol = sorted(sol, reverse = False)  # Sort from smallest to largest            
print('Stability points are', np.rad2deg(sol), '[deg]')
    
with plt.style.context('ggplot'):
    plt.plot(np.rad2deg(lamda), acc_ew, color='g') 
    plt.xlim([0, max(np.rad2deg(lamda))])
    # plt.axhline(y=0, color = 'k', linestyle='-')
    for i in range(len(sol)):
        plt.axvline(x=np.rad2deg(sol[i]), color = 'k', linestyle='--')
    plt.xlabel("Right Ascension [rad]")
    plt.ylabel("East-West acceleration [$m/s^2$]")
    plt.title("East-West acceleration due to J2,2 term in GEO")
    plt.tight_layout()
    plt.show()