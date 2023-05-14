"""
Calculate the energy in joules for a green (λ = 525 nm) light photon as 
well as for gamma ray (λ = 10-12 m) and for radio wave (λ = 1 m).
"""

# Importing Libraries
import numpy as np

# Declaring constants
h = 6.62607004e-34
c = 299792458

# Calculating energy for green photon
E_green = h * c / (525e-9)

# Calculating energy for gamma ray
E_gamma = h * c / (10e-12)

# Calculating energy for radio wave
E_radio = h * c / (1)

# Printing results
print("Energy for green photon:", E_green, "J")
print("Energy for gamma ray:", E_gamma, "J")
print("Energy for radio wave:", E_radio, "J")