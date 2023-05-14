#%% This file contains a series of functions to be used for attitude determination and control
import numpy as np
from math import sin, cos, acos
import matplotlib.pyplot as plt
from sympy.physics.mechanics import ReferenceFrame, dot, cross

#%% Reference frame
class IJKReferenceFrame(ReferenceFrame):
    def __init__(self, name):
        super().__init__(name, latexs=['\mathbf{%s}_{%s}' % (
            idx, name) for idx in ("i", "j", "k")])
        self.i = self.x
        self.j = self.y
        self.k = self.z

################################################################################
################################################################################
#%% Basic functions
################################################################################
def skew(vector):
    """
    Function to calculate a 3x3 skew-symmetric matrix
    """
    return np.array([[0, -vector[2], vector[1]],
                     [vector[2], 0, -vector[0]],
                     [-vector[1], vector[0], 0]])

################################################################################
def print_matrix(array, decimals=3):
    """
    A function to just print a matrix in Latex form. It just looks nicer.
    """
    from IPython.display import display, Latex, Math
    matrix = ''
    for row in array:
        try:
            for number in row:
                matrix += f'{round(number,decimals)}&'
        except TypeError:
            matrix += f'{round(row,decimals)}&'
        matrix = matrix[:-1] + r'\\'
    display(Math(r'\begin{bmatrix}'+matrix+r'\end{bmatrix}'))

################################################################################
def find_eig(J):
    """
    A function to find the eigenvalues and eigenvectors of a matrix.
    """
    Eigen = np.linalg.eigh(
        J)  # Other methods include la.eig (import scipy.linalg as la) or np.linalg.eig, but give a different order for eigenvectors
    Eigenvalues = Eigen[0]
    # The column v[:, i] is the normalized eigenvector corresponding to the eigenvalue w[i]. Will return a matrix object if a is a matrix object.
    Eigenvectors = Eigen[1]
    return Eigenvalues, Eigenvectors

################################################################################
################################################################################
#%% DCM and Euler angles
################################################################################
def DCM_A_B(A, B):
    """
    DCM between two reference frames, A with respect to B
    """
    return np.dot(A.T, B)
    
################################################################################
#def DCM_321(Eul_ang): !!!! This is wrong. It was wrong in the solutions
# INSTEAD USE DCM() in symbolic library with 'num' set as mode to get DCM for all 
# rotations
#    """
#    Direction cosine matrix for the  3-2-1 Euler angles
#    """
#    return np.array([[cos(Eul_ang[1])*cos(Eul_ang[2]), cos(Eul_ang[1])*sin(Eul_ang[2]), -sin(Eul_ang[1])],
#                     [sin(Eul_ang[0])*sin(Eul_ang[1])*cos(Eul_ang[2]) - cos(Eul_ang[0])*sin(Eul_ang[2]), sin(Eul_ang[0])
#                      * sin(Eul_ang[1])*sin(Eul_ang[2]) + cos(Eul_ang[0])*cos(Eul_ang[2]), sin(Eul_ang[0])*cos(Eul_ang[1])],
#                     [cos(Eul_ang[0])*sin(Eul_ang[1])*cos(Eul_ang[2]) + sin(Eul_ang[0])*sin(Eul_ang[2]), cos(Eul_ang[0])*sin(Eul_ang[1])*sin(Eul_ang[2]) - sin(Eul_ang[0])*cos(Eul_ang[2]), cos(Eul_ang[0])*cos(Eul_ang[1])]])
#
################################################################################
def Eul_ang(DCM):
    """
    Euler angles from DCM
    """
    return np.array([np.rad2deg(np.arctan(DCM[1, 2] / DCM[2, 2])),
                     np.rad2deg(-np.arcsin(DCM[0, 2])),
                     np.rad2deg(np.arctan(DCM[0, 1] / DCM[0, 0]))])

################################################################################
def diff_kinem_321_Euler(time, theta, w1=None, w2=None, w3=None):
    """
    Differential kinematics for a 3-2-1 Euler angle sequence.
    - Needs to be edited with the correct angular velocity vector for the user's 
      specific use case.
    """
    if w1 is None and w2 is None and w3 is None:
        w = np.array([np.sin(0.1*time), 0.01, np.cos(0.1*time)]) * np.deg2rad(5)
    else:
        w = np.array([w1, w2, w3])

    dot_angles = np.dot((1/cos(theta[1]) * np.array([
        [cos(theta[1]), sin(theta[0])*sin(theta[1]),
         cos(theta[0])*sin(theta[1])],
        [0,           cos(theta[0])*cos(theta[1]), -
         sin(theta[0])*cos(theta[1])],
        [0,           sin(theta[0]),             cos(theta[0])]
    ])), w)
    return dot_angles

################################################################################
def plot_euler(time, dot_angles, ylabel='Angle (deg)', ylabel1='Yaw Angle (deg)', ylabel2='Pitch Angle (deg)', ylabel3='Roll Angle (deg)'):

    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.plot(time, dot_angles[0, :], label=ylabel1)
    plt.plot(time, dot_angles[1, :], label=ylabel2)
    plt.plot(time, dot_angles[2, :], label=ylabel3)
    plt.xlabel('Time (s)')
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid()
    plt.show()

    with plt.style.context('seaborn-notebook'):
        fig = plt.figure(figsize=(10, 5))
        plt.subplot(311)
        plt.plot(time, dot_angles[0, :], label='yaw')
        plt.grid()
        plt.xlabel('Time (s)')
        plt.ylabel(ylabel1)
        plt.xlim(0, max(time))
        plt.ylim(min(dot_angles[0, :])*1.1, max(dot_angles[0, :])*1.1)

        plt.subplot(312)
        plt.plot(time, dot_angles[1, :], label='pitch', color='r')
        plt.grid()
        plt.xlabel('Time (s)')
        plt.ylabel(ylabel2)
        plt.xlim(0, max(time))
        plt.ylim(min(dot_angles[1, :])*1.1, max(dot_angles[1, :])*1.1)

        plt.subplot(313)
        plt.plot(time, dot_angles[2, :], label='roll', color='g')
        plt.grid()
        plt.xlabel('Time (s)')
        plt.ylabel(ylabel3)
        plt.xlim(0, max(time))
        plt.ylim(min(dot_angles[2, :])*1.1, max(dot_angles[2, :])*1.1)

        fig.tight_layout()
        plt.show()
################################################################################
################################################################################
#%% Quaternions
################################################################################
def Eigenaxis_rot(C):
    """
    Input: 
     - direction cosine Euler rotation matrix
    Output: 
     - principle Euler eigenaxis rotation angle phi
    """
    return acos(1/2*(C[0, 0] + C[1, 1] + C[2, 2] - 1))

################################################################################
def Eigenaxis_e(C):
    """
    Input: 
     - direction cosine Euler rotation matrix
    Output: 
     - Find the principle Euler eigenaxis e
    """
    phi = Eigenaxis_rot(C)
    return (1/(2*sin(phi))) * np.array([C[1, 2] - C[2, 1], C[2, 0] - C[0, 2], C[0, 1] - C[1, 0]])

################################################################################
def Euler_eigenaxis(C):
    """
    Find the principle Euler eigenaxis rotation angle phi and principle Euler eigenaxis e
    """
    return Eigenaxis_rot(C), Eigenaxis_e(C)

################################################################################
def Eigenaxis_rotMAT(e, phi):
    """
    Input:
     - principle Euler eigenaxis e
     - principle Euler eigenaxis rotation angle phi
    Output:
     - The rotation matrix for a given principle Euler eigenaxis rotation
    """
    C = cos(phi) * np.eye(3) + np.dot((1-cos(phi)), (e.T*e)) - sin(phi)*skew(e)
    return C

################################################################################
def DCM_to_Quaternion(C):
    """
    DCM in Euler angles to Quaternion
    """
    q = np.zeros(4)
    q[3] = 1/2 * (1 + C[0, 0] + C[1, 1] + C[2, 2])**0.5
    q[0:3] = (1/(4*q[3])) * np.array([C[1, 2] - C[2, 1],
                                  C[2, 0] - C[0, 2],
                                  C[0, 1] - C[1, 0]])
    return q#/np.linalg.norm(q)

################################################################################
def Eigenaxis_to_Quaternion(e, phi):
    """
    Eigenaxis to Quaternion
    """
    q = np.zeros(4)
    q[3] = cos(phi/2)
    q[0:3] = np.array([e[0], e[1], e[2]]) * sin(phi/2)
    return q#/np.linalg.norm(q)

################################################################################
def Quaternion_to_DCM(q):
    """
    Quaternion to DCM
    """
    q1, q2, q3, q4 = q[0], q[1], q[2], q[3]
    dcm = np.zeros((3, 3))

    dcm[0, 0] = 1 - 2*(q2**2 + q3**2)
    dcm[0, 1] = 2*(q1*q2 + q3*q4)
    dcm[0, 2] = 2*(q1*q3 - q2*q4)
    dcm[1, 0] = 2*(q2*q1 - q3*q4)
    dcm[1, 1] = 1 - 2*(q3**2 + q1**2)
    dcm[1, 2] = 2*(q2*q3 + q1*q4)
    dcm[2, 0] = 2*(q3*q1 + q2*q4)
    dcm[2, 1] = 2*(q3*q2 - q1*q4)
    dcm[2, 2] = 1 - 2*(q1**2 + q2**2)
    
    return dcm

################################################################################
def diff_kinem_Quaternion(time, q, w1=None, w2=None, w3=None, w4=None):
    """
    Differential kinematics for quaternion.
    - Needs to be edited with the correct angular velocity vector for the user's 
      specific use case.
    """
    if w1 is None and w2 is None and w3 is None and w4 is None:
        w = np.array(
                [np.sin(0.1*time), 0.01, np.cos(0.1*time), 0]) * np.deg2rad(50)
    else:
        w = np.array([w1, w2, w3, w4])

    dot_q = np.dot( (1/2) * np.array([
                [0,     w[2], -w[1], w[0]],
                [-w[2],   0,   w[0], w[1]],
                [w[1],  -w[0],   0,  w[2]],
                [-w[0], -w[1], -w[2],   0]
                ]), q)
    return dot_q

################################################################################
def plot_quaternion(t, q):
    """
    Plot the quaternion.
    """
    with plt.style.context('seaborn-notebook'):
        fig, ax = plt.subplots(2, 2, sharex=True, figsize=(10, 10))
        fig.suptitle('Quaternions over time')
        ax[0, 0].plot(t, q[0, :], label='$q_0$')
        ax[0, 1].plot(t, q[1, :], label='$q_1$', color='r')
        ax[1, 0].plot(t, q[2, :], label='$q_2$', color='g')
        ax[1, 1].plot(t, q[3, :], label='$q_3$', color='y')
        for ax in ax.flat:
            ax.set(xlabel='Time (s)', ylabel='Quaternion')
            ax.legend()
            ax.label_outer()
            ax.grid()
            ax.set_xlim(0, 60)
        fig.tight_layout()
        plt.show()

################################################################################
def plot_quaternion_constraint(t, q):
    """
    Plot the quaternion constraint | q | = 1.
    """
    q_mag = q[0, :]**2 + q[1, :]**2 + q[2, :]**2 + q[3, :]**2
    with plt.style.context('seaborn-notebook'):
        fig, ax = plt.subplots(1, 1, sharex=True, figsize=(10, 10))
        fig.suptitle('Quaternion constraint $|q|= 1$')
        ax.plot(t, q_mag, label='$|q|$')
        ax.set(xlabel='Time (s)', ylabel='$|q|$')
        ax.legend()
        ax.grid()
        ax.set_xlim(0, 60)
        fig.tight_layout()
        plt.show()

################################################################################
################################################################################
#%% Attitude Determination
################################################################################

# TRIAD Method
def triad(b1, b2, r1, r2):
    """
    Input: 
     - body frame vectors: b1 and b2  
     - reference frame vectors: r1 and r2
    Output: 
     - body frame triad: t1b, t2b and t3b  
     - reference frame triad: t1i, t2i and t3i
     - rotation matrix given by triad method
    """
    # Normalize the vectors
    b1 = b1 / np.linalg.norm(b1)
    b2 = b2 / np.linalg.norm(b2)
    r1 = r1 / np.linalg.norm(r1)
    r2 = r2 / np.linalg.norm(r2)

    # Calculate body coordinates
    t1b = b1
    t2b = np.cross(b1, b2)/np.linalg.norm(np.cross(b1, b2))
    t3b = np.cross(t1b, t2b)

    rot_tb = np.array([t1b, t2b, t3b]).T  # Rotational matrix

    # Calculate inertial coordinates
    t1r = r1
    t2r = np.cross(r1, r2)/np.linalg.norm(np.cross(r1, r2))
    t3r = np.cross(t1r, t2r)

    rot_tr = np.array([t1r, t2r, t3r]).T  # Rotational matrix

    # Calculate rotation matrix
    Cbr = np.dot(rot_tb, rot_tr.T)

    return rot_tb, rot_tr, Cbr

################################################################################
# q-METHOD
def q_method(b, RF, weights=None):
    """
    Input: 
     - body frame vectors: vb, ub, ..., un
     - reference frame vectors: vi, ui, ..., un
     - weights: weights for each vector. Input as an array
    Output: 
     - B matrix (3x3 matrix)
     - K matrix (4x4 matrix) with components: K11, k12, k22
     - Max Eigenvalue and corresponding Eigenvector
     - C Rotation matrix
    """
    if weights == None:
        weights = np.ones(len(b))

    B = np.zeros((3, 3))
    for i in range(len(b)):
        # Normalize the vectors
        b[i] = b[i] / np.linalg.norm(b[i])
        RF[i] = RF[i] / np.linalg.norm(RF[i])
        # Find B matrix
        Bb = weights[i]*np.outer(b[i], RF[i])
        B += Bb

    k22 = np.trace(B)
    K11 = B + B.T - k22*np.eye(3, 3)
    k12 = np.array(
        [(B[1, 2] - B[2, 1]), (B[2, 0] - B[0, 2]), (B[0, 1] - B[1, 0])]).T

    K = np.zeros((4, 4))
    K[0:3, 0:3] = K11
    K[0:3, 3] = k12.T
    K[3, 0:3] = k12.T
    K[3, 3] = k22

    Eigenvalues, Eigenvectors = find_eig(K)
    max_Eigenvalue = np.max(Eigenvalues)
    max_Eigenvector = Eigenvectors[:, np.where(
        Eigenvalues == np.max(Eigenvalues))]
    C = Quaternion_to_DCM(max_Eigenvector)

    return B, k22, K11, k12, K, max_Eigenvalue, max_Eigenvector, C

################################################################################
# QUEST-METHOD
def QUEST(b, RF, weights=None):
    """
    Inputs:
        b_vec: Body vectors
        RF_vec: Reference Frame vectors
        weights: Weights for the different components
    Outputs:
        S: S matrix (3x3 matrix)
        K: K matrix (4x4 matrix) with components: K11, k12, k22
        p: p vector (3x1 vector)
        q: quaternion
        C: Rotation matrix
    """
    # Weights
    if weights == None:
        weights = np.ones(len(b))

    B = np.zeros((3, 3))
    for i in range(len(b)):
        # Normalize the vectors
        b[i] = b[i] / np.linalg.norm(b[i])
        RF[i] = RF[i] / np.linalg.norm(RF[i])
        # Find B matrix
        Bb = weights[i]*np.outer(b[i], RF[i])
        B += Bb

    S = B + B.T

    k22 = np.trace(B)
    K11 = B + B.T - k22*np.eye(3, 3)
    k12 = np.array(
        [(B[1, 2] - B[2, 1]), (B[2, 0] - B[0, 2]), (B[0, 1] - B[1, 0])]).T

    K = np.zeros((4, 4))
    K[0:3, 0:3] = K11
    K[0:3, 3] = k12.T
    K[3, 0:3] = k12.T
    K[3, 3] = k22

    Eigenvalues, _ = find_eig(K)
    max_Eigenvalue = np.max(Eigenvalues)

    p = np.dot(np.linalg.inv(
        np.array([(max_Eigenvalue + k22)*np.eye(3, 3) - S])), k12).T
    p4 = np.array([p[0], p[1], p[2], 1], dtype=object).T

    q = 1 / (np.sqrt(1 + p.T @ p)) * p4
    q = np.array(q[0])

    C = Quaternion_to_DCM(q)
    return S, K, p, q, C

################################################################################
# Solve Kinematic Differential Equation
def solve_KDE(input, time_range=[0, 60], time_array=np.linspace(0, 60, 244), solver='E321', w=None):
    """
    Solve the Kinematic Differential Equation
    Input:
        input: Input data (quaternion or Euler angles)
        time_range: Time range i.e. [0,60]
        time_array: Time array i.e. time = np.linspace(0, 60, 244)
        solver: Solver: either "q" for quaternion or "E321" for Euler angles in 321 convention,
                ATM the solver is limited for just 321 Euler sequence, for a different sequence
                use the function solve_KDE_Euler()
    Output:
        output: Output data (solution time, solution_data)
    Note (!):
        w: angular velocity (3x1 vector). Needs to be edited in the diff_kimen function 
           with the correct angular velocity vector for the user's 
           specific use case.
    """
    from scipy.integrate import solve_ivp

    if solver=='E321':
        if w is None:
            sol = solve_ivp(diff_kinem_321_Euler, time_range,
                            input, t_eval=time_array)
        else:
            w1, w2, w3 = w[0, :], w[1, :], w[2, :]
            sol = solve_ivp(diff_kinem_321_Euler, time_range,
                            input, t_eval=time_array, args=(w1, w2, w3))
    #elif solver=='E':
    #    sol = solve_ivp(diff_kinem_Euler, time_range, 
    #                    input, t_eval=time_array)
    elif solver=='q':
        if w is None:
            sol = solve_ivp(diff_kinem_Quaternion, time_range, 
                            input, t_eval=time_array)
        else:
            w1, w2, w3, w4 = w[0,:], w[1,:], w[2,:], w[3,:]
            sol = solve_ivp(diff_kinem_Quaternion, time_range, 
                            input, t_eval=time_array, args=(w1, w2, w3, w4))
    else:
        print('Solver not found')
    
    return sol.t, sol.y

################################################################################
def solve_Euler_KDE(t, theta_angles, rot1, rot2, rot3, w=None, invorder=True, plot=True):
    """
    Differential kinematics for any Euler angle sequence.

    Input:
    - t: time. i.e. t = np.linspace(0, 60, 61)
    - theta_angles: Euler angles.  i.e. thetanum = np.deg2rad([80, 30, 40])
    - rot1, rot2, rot3: rotation sequence. Same convention as in diff_kinem_Euler_theta_matrix
    - w: Angular velocity vector as a sympy Matrix. i.e. Matrix([w1, w2, w3])
    - invorder: If True, the order of the rotation sequence is inverted.
                i.e. invorder=True:
                theta1 -> rotation along the 3rd rotation axis 
                theta2 -> rotation along the 2nd rotation axis 
                theta3 -> rotation along the 1st rotation axis 
    (!) Note: Need to review the solver. It is not working properly.
    """
    ## Import symbolic library
    from sys import path
    path.append(
        "c:\\Users\\diego\\Dropbox\\Academic\\MEng Space Systems\\3. DOCA\\ADCS functions")
    import ADCS_Functions_sym as adcs_sym
    import sympy as sym

    # Get the symbolic solution to the differential kinematic equation:
    sol= adcs_sym.diff_kinem_Euler(
        rot1, rot2, rot3, w=w, invorder=invorder)

    # Get the symbols needed back
    time, x, y, z = sym.symbols('t theta_1 theta_2 theta_3')
    dotx, doty, dotz = sym.symbols(r"\dot{\theta_1} \dot{\theta_2} \dot{\theta_3}")

    # Get the numerical solution to the differential kinematic equation:
    dotx = sol[dotx].subs(x, theta_angles[0]).subs(
        y, theta_angles[1]).subs(z, theta_angles[2])
    doty = sol[doty].subs(x, theta_angles[0]).subs(
        y, theta_angles[1]).subs(z, theta_angles[2])
    dotz = sol[dotz].subs(x, theta_angles[0]).subs(
        y, theta_angles[1]).subs(z, theta_angles[2])

    dot_angles = sym.Matrix([dotx, doty, dotz])

    dot_angles_num = np.zeros((3, len(t)))
    for i in range(len(t)):
        dot_angles_num[0, i] = dot_angles[0].subs(time, t[i])
        dot_angles_num[1, i] = dot_angles[1].subs(time, t[i])
        dot_angles_num[2, i] = dot_angles[2].subs(time, t[i])

    if plot==True:
        plot_euler(t, np.rad2deg(dot_angles_num))
    return dot_angles_num

################################################################################
################################################################################
#%% Attitude Dynamics
################################################################################
################################################################################
def axisym_torquefree(I, w):
    """
    This function calculates the axisymatic angular momentum, nutation angle and precession rate.
    Assumptions: axisymatic body, torque-free motion
    """
    IT = (I[0] + I[1])/2  # Since it is axisymmetric it is equal to  = Ix  = Iy
    h = I*w
    hT = np.linalg.norm(np.array([h[0], h[1]]))
    h_mag = np.linalg.norm(h)
    gamma = np.arcsin(hT/h_mag)
    Omega_p = h_mag / IT
    return h, hT, h_mag, gamma, Omega_p

################################################################################
def eulerEq(time, w, I1, I2, I3):
    """
    Differential kinematics for a 3-2-1 Euler angle sequence.

    Solve it in the form: sol = solve_ivp(eulerEq, [0, 100], wIC, t_eval=time, args=(I1,I2,I3))
    """
    #w = np.array([0.05, 0.02, -0.02])  # rad/s

    I = np.array([[I1, 0, 0], [0, I2, 0], [0, 0, I3]])
    dot_w = np.linalg.inv(I) @ (np.cross(-w, (I @ w)))

    return dot_w

################################################################################
def gravity_gradient_Tb(I, R_0b_vec, R0G_vec=None): 
    """
    This function calculates the gravity gradient torque in body coordinates.
    
    Input:
    - I:  spacecraft inertia matrix
    - R0G: spacecraft orbital position vector in ECI coordinates
    - R_0b: spacecraft position vector in body coordinates

    Note! Units are in km
    """
    ## Initialise
    mu = 3.986e5  # km^3/s^2
    # magnitude of spacecraft orbital position vector in body coordinates:
    if R0G_vec is None:
        R0 = np.linalg.norm(R_0b_vec)  # km
    if R0G_vec is not None:
        R0 = np.linalg.norm(R0G_vec)  # km
    Tgg_b = (3*mu) / (R0**5) * skew(R_0b_vec) @ I @ R_0b_vec

    return Tgg_b

################################################################################
def gravity_gradient_Tb_aprox(I, theta, R0=None, R_0b_vec=None):  
    """
    This function calculates the gravity gradient torque in body coordinates.
    Uses small angle aproximation (only for 321 sequence)
    
    Input:
    - I:  spacecraft inertia matrix
    - theta: spacecraft Euler angles
    - R0: magnitude of spacecraft orbital position vector
    - R_0b: spacecraft position vector in body coordinates

    Note! Units are in km
    """
    ## Initialise
    mu = 3.986e5  # km^3/s^2
    if R0 is None and R_0b_vec is not None:
        # magnitude of spacecraft orbital position vector in body coordinates:
        R0 = np.linalg.norm(R_0b_vec)  # km
    omega2 = mu/(R0**3)

    Tgg_b = (3*omega2) * np.array([ (I[2] - I[1])*theta[0], (I[2] - I[0])*theta[1], 0 ])

    return Tgg_b

################################################################################
################################################################################
#%% Complete Attitude Kinematics
################################################################################

def attitude_kinematics_321(time, x, I1, I2, I3):  # (time=None, ics=None):
    """
    Differential kinematics equation for a 3-2-1 Euler angle sequence.
    """
    ## Arguments
    I = np.array([I1, I2, I3])
    w = np.array([x[0], x[1], x[2]])
    theta = np.array([x[3], x[4], x[5]])

    ## Kinematics Equations
    dot_w1 = - ((I[2]-I[1]) / I[0]) * x[1] * x[2]
    dot_w2 = - ((I[0]-I[2]) / I[1]) * x[0] * x[2]
    dot_w3 = - ((I[1]-I[0]) / I[2]) * x[1] * x[0]
#
    dot_angles = np.dot((1/cos(theta[1]) * np.array([
        [cos(theta[1]), sin(theta[0])*sin(theta[1]),
         cos(theta[0])*sin(theta[1])],
        [0,           cos(theta[0])*cos(theta[1]), -
         sin(theta[0])*cos(theta[1])],
        [0,           sin(theta[0]),             cos(theta[0])]
    ])), w)

    dot_x = np.concatenate((dot_w1, dot_w2, dot_w3, dot_angles), axis=None)

    return dot_x

################################################################################
def solve_attitude_kimenatics_321(ics, I1, I2, I3, time_range=[0, 100], time_array=np.linspace(0, 100, 101), plot=True):
    """
    Solve the attitude kinematics for a 3-2-1 Euler angle sequence.

    Inputs:
    - ics: initial conditions. i.e. np.concatenate([wIC, thetaIC])
    - I1, I2, I3: inertia tensor elements
    - time_range: time vector. i.e. [0, 100]
    - time_array: time array. i.e. np.linspace(0, 100, 101)
    - plot: boolean. If True, plot the results.
    """
    from scipy.integrate import solve_ivp

    ## Solve
    sol = solve_ivp(attitude_kinematics_321, time_range, ics, t_eval=time_array, method='RK45', args=(I1, I2, I3))

    ## Get values
    time = sol.t
    omega = sol.y[0:3, :]
    theta = sol.y[3:6, :]

    if plot==True:
        plot_euler(time, omega, ylabel='$\omega$', ylabel1='$\omega_1$',
                        ylabel2='$\omega_2$', ylabel3='$\omega_3$')
        plot_euler(time, theta)

    return time, omega, theta

################################################################################
def stability_analysis(time, x, w, I1, I2, I3):
    """
    Stability analysis of attitude motion for gravity gradient stabilization. 
    
    (Three coupled linear differential equations).

    In matrix form x'' = A x' + B x

    For stability I2>I1>I3
    """
    ## Variables
    theta = np.array([x[0], x[1], x[2]])
    dottheta = np.array([x[3], x[4], x[5]])

    ## Initialize
    A = np.array([[0, 0, ((I1-I2+I3)*w)/I1],
                  [0, 0, 0], 
                  [((I1-I2+I3)*w)/I3, 0, 0]])
    dotx_arr = np.array([dottheta[0], dottheta[1], dottheta[2]])

    B = np.array([[(4*w**2*(I2-I3))/I1, 0, 0],
                  [0, -(3*w**2*(I1-I3))/I2, 0], 
                  [0, 0, (w**2*(I2-I1))/I3]])
    theta_arr = np.array([theta[0], theta[1], theta[2]])


    ## Diff equations
    dotdot_theta1 = A @ dotx_arr + B @ theta_arr

    sols = np.concatenate((dotdot_theta1, dotx_arr), axis=None)
    return sols

################################################################################
def solve_stability_analysis(ics, w, I1, I2, I3, time_range=[0, 100], time_array=np.linspace(0, 100, 101), plot=True):
    """
    Solve the stability analysis for gravity gradient stabilization.

    Inputs:
    - ics: initial conditions. i.e. np.concatenate([wIC, thetaIC])
    - w: angular velocity of the body
    - I1, I2, I3: inertia tensor elements
    - time_range: time vector. i.e. [0, 100]
    - time_array: time array. i.e. np.linspace(0, 100, 101)
    - plot: boolean. If True, plot the results.
    """
    from scipy.integrate import solve_ivp

    ## Check if stable
    if I2 > I1 > I3:
        print("Stable")
    else:
        print("Unstable")

    ## Solve
    sol = solve_ivp(stability_analysis, time_range, ics,
                    t_eval=time_array, method='RK45', args=(w, I1, I2, I3))

    ## Get values
    time = sol.t
    dottheta = sol.y[0:3, :]
    theta = sol.y[3:6, :]

    if plot == True:
        plot_euler(time, dottheta, ylabel='$\dot{\omega}$', ylabel1='$\\dot{\omega_1}$',
                        ylabel2='$\dot{\omega_2}$', ylabel3='$\dot{\omega_3}$')
        plot_euler(time, theta)

    return time, theta, dottheta

################################################################################
#%% Active Attitude Control Functions:
################################################################################
def xiomega_n_f(ts):
    """
    Use for closed-loop system only.

    Returns the product of the undamped natural frequency (omega_n) and the damping ratio (xi).
    Calculated from settling time formula (page 52 of ADCSX.pdf)
    """
    xi_omega_n = 4.4/ts
    
    return xi_omega_n

################################################################################
def close_loop_TF_solve(ts, Mp, wn):
    """"
    Solve for the TF constants for a closed loop system with PD control.
    Inputs:
        - ts: sampling time
        - Mp: maximum overshoot
    Outputs:
        - returns the TF and the constants:
        - the product of the undamped natural frequency (omega_n) and the damping ratio (xi)
        - the damped natural frequency (omega_d)
        - TF(tf): transfer function for a close-loop system with PD control
    """
    ## Import symbolic library
    from sys import path
    path.append(
        "c:\\Users\\diego\\Dropbox\\Academic\\MEng Space Systems\\3. DOCA\\ADCS functions")
    import ADCS_Functions_sym as adcs_sym
    import sympy as sym

    xiomega_n = xiomega_n_f(ts)
    omega_d = sym.solve(sym.Eq(adcs_sym.max_overshoot_CLTF(), Mp),
                        sym.Symbol('omega_d'))[0].subs(sym.Symbol('xi')*sym.Symbol('omega_n'), xiomega_n)
    TF = adcs_sym.close_loop_TF()
    TF = TF.subs(sym.Symbol('\omega_n')*sym.Symbol('xi'), xiomega_n).subs(sym.Symbol('\omega_d'),
                                                                  omega_d).subs(sym.Symbol('\omega_n'), wn)
    return xiomega_n, omega_d, TF

################################################################################
def closed_loop_poles(xiomega_n, omega_d, plot=True):
    """
    Returns the poles of a closed-loop system.

    Input:
     - xiomega_n: product of the undamped natural frequency (omega) and the damping ratio (xi).
     - omega_d: damped frequency of the system.
    """
    s = - xiomega_n + omega_d*1j
    s_conj = s.conjugate()

    poles = np.array([s, s_conj]).astype(complex)

    if plot is True:
        poles_re = [ele.real for ele in poles]
        # extract imaginary part
        poles_im = [ele.imag for ele in poles]
#
        # plot the complex numbers
        plt.scatter(poles_re, poles_im, marker='x', color='red')
        plt.axhline(0, color='k')  # x = 0
        plt.axvline(0, color='k')  # y = 0
        plt.title('s-plane')
        plt.ylabel('j')
        plt.xlabel('$\sigma$')
        plt.grid()
        plt.show()

    return poles

################################################################################
def PD_coefs_CL(I, poles=None, xiomega_n=None, xi=None, omega_n=None):
    """
    Returns the PD coefficients for the closed-loop system.
    Input:
        I: inertia, then:
        (Option 1):
            poles: list of poles
            xiomega_n: product of natural frequency and damping ratio
        (Option 2):
            xi: damping ratio
            omega_n: natural frequency
    Output:
        Kp: proportional gain
        Kd: derivative gain
    """
    if poles is not None and xiomega_n is not None and xi is None and omega_n is None:
        Kp = I * abs(poles[0])**2
    elif poles is None and xiomega_n is None and xi is not None and omega_n is not None:
        Kp = I * omega_n**2
    else:
        raise ValueError('Invalid input')

    Kd = 2*I*xiomega_n
    return Kp, Kd

################################################################################
def rise_time_CLTF(xi, omega_d):
    """"
    Determine the rise time for a closed loop system with a given xi and omega_d

    Input:
    - xi: damping ratio
    - omega_d: damped frequency
    """
    beta = np.arctan(np.sqrt(1-xi**2)/xi)
    return (np.pi-beta)/omega_d
