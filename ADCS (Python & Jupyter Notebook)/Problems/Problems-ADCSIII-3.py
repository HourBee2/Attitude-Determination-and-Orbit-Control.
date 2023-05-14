#%% Import libraries
import numpy as np
import matplotlib.pyplot as plt

from astropy import units as u
from astropy import constants
# from astropy.coordinates import get_body_barycentric
# from astropy.time import Time

# from poliastro.bodies import Sun, Jupiter

from scipy.optimize import fsolve, root_scalar

#%% Define the system
m1 = constants.M_sun
m2 = constants.M_jup
r1cm = 0 * u.m      # Radius of sun CM
r2cm = (7.7834e8*u.km).to(u.m)    # Radius of Jupiter CM
CM = (m1*r1cm + m2*r2cm) / (m1 + m2) # Radius of CM for Sun-Jupyter system
r1 = CM
r2 = r2cm - CM
r12 = r1 + r2
omega = np.sqrt(constants.G * (m1+m2)/r12**3)

#%% Lagrange points using poliastro
# from poliastro.threebody import restricted
## Sun-Jupyter system ##
# Lagrange = restricted.lagrange_points_vec(constants.M_sun, 0*u.m, constants.M_jup, 
                                # 778e9*u.m,  np.array([0,1,0])*u.m/u.m)
# Lagrange = restricted.lagrange_points(778e9*u.m, constants.M_sun, constants.M_jup)

# with plt.style.context('ggplot'):
#     plt.plot( Lagrange/1e3, 'o', color='r') 

#%% Collinear Lagrange Points (y=0)

# x = np.linspace(0, 1.6*r2cm,1000)

### Method 1: Function with units ###
def col_Lagrange(x):
    return omega**2 * x*u.m - ((constants.G*m1*(x*u.m-r1))/(np.abs(x*u.m-r1)**3) 
                                + (constants.G*m2*(x*u.m+r2))/(np.abs(x*u.m+r2)**3))
L1 = fsolve(col_Lagrange, 7e11) * u.m
L2 = fsolve(col_Lagrange, 8.9e11) * u.m
L3 = fsolve(col_Lagrange, -7.78e11) * u.m

unstable_Lagrange = np.array([L1, L2, L3]) * u.m

## Method 2: lambda tells the func that x is the value to solve for ###
# func = lambda x : omega**2 * x*u.m - ((constants.G*m1*(x*u.m-r1))/(np.abs(x*u.m-r1)**3)
                                      # + (constants.G*m2*(x*u.m+r2))/(np.abs(x*u.m+r2)**3))
## Method 3: Function without units (units dont work with root_scalar) ###
# def col_Lagrange(x):
#     return (omega*u.s)**2 * x*u.m/u.m - ((constants.G*u.kg*u.s**2/u.m**3*m1/u.kg*(x/u.m*u.m-r1/u.m))/(np.abs(x/u.m*u.m-r1/u.m)**3) 
#                                 + (constants.G*u.kg*u.s**2/u.m**3*m2/u.kg*(x/u.m*u.m+r2/u.m))/(np.abs(x/u.m*u.m+r2/u.m)**3))
                                  
# L1 = root_scalar(col_Lagrange, bracket=[6e11, 8e11]) 
# L2 = root_scalar(col_Lagrange, bracket=[8e11, 9e11]) 
# L3 = root_scalar(col_Lagrange, bracket=[-8e11, -6e11]) 

# unstable_Lagrange = np.array([L1.root, L2.root, L3.root]) * u.m
                                  
### Print results ###
print(unstable_Lagrange)

#%% Plot Lagrange points
with plt.style.context('ggplot'):
    bodies = np.array([r1cm*1e-3/u.m, r2cm*1e-3/u.m])
    plt.plot( bodies, np.zeros(len(bodies)), 'o', color='y',  markersize=20) 
    plt.plot( unstable_Lagrange/1e3, np.zeros(len(unstable_Lagrange)), 'x', color='r') 
    plt.xlabel("x-coord [$km$]")
    plt.ylabel("y-coord [$km$]")
    plt.title("Position of the Lagrange points for the Sun-Jupyter system")
    plt.tight_layout()
    plt.show()