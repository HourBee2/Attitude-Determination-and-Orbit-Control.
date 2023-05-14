import numpy as np
from sys import path
path.append(
    "c:\\Users\\diego\\Dropbox\\Academic\\MEng Space Systems\\3. DOCA\\ADCS functions")
import ADCS_Functions as adcs

#%% Data
v1b = np.array([0.8273, 0.5541, -0.0920]).T
v2b = np.array([-0.8285, 0.5522, -0.0955]).T

v1i = np.array([-0.1517, -0.9669, 0.2050]).T
v2i = np.array([-0.8393, 0.4494, -0.3044]).T

#%% Functions

# TRIAD Method

#%% Calculate
rot_triad_b, rot_triad_i, Cbi_triad = adcs.triad(v1b, v2b, v1i, v2i)


t1b = rot_triad_b[:, 0]
t2b = rot_triad_b[:, 1]
t3b = rot_triad_b[:, 2]

t1i = rot_triad_i[:, 0]
t2i = rot_triad_i[:, 1]
t3i = rot_triad_i[:, 2]

print('t1b = ', t1b), print('t2b = ', t2b), print('t3b = ', t3b)
print('t1i = ', t1i), print('t2i = ', t2i), print('t3i = ', t3i)
print('\n')
print('Cbi =', Cbi_triad)

# Check if correct
print(np.allclose(Cbi_triad@v1i, v1b))

#%% Test with q-method and QUEST
b = np.array([v1b, v2b])
RF = np.array([v1i, v2i])

# q-method
B, k22, K11, k12, K, max_Eigenvalue, max_Eigenvector, q_method_C = adcs.q_method(b, RF)
print('\n'); print(q_method_C)

# QUEST
S, K, p, q, QUEST_C = adcs.QUEST(b, RF)
print('\n'); print(QUEST_C)

