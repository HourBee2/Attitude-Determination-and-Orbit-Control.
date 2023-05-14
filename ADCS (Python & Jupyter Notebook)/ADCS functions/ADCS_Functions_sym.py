#%% This file contains a series of functions to be used for attitude determination and control

#%% Import libraries
import numpy as np
from sympy.physics.mechanics import ReferenceFrame, dot, cross
from sympy import *
init_printing(use_latex='mathjax')

#%% Reference frame
class IJKReferenceFrame(ReferenceFrame):
    def __init__(self, name):
        super().__init__(name, latexs=['\mathbf{%s}_{%s}' % (
            idx, name) for idx in ("i", "j", "k")])
        self.i = self.x
        self.j = self.y
        self.k = self.z

################################################################################
#%% Basic functions
################################################################################
def sym2num(sim_vec, simsubs, numsubs, simsubs2=None, numsubs2=None, simsubs3=None, numsubs3=None):
    """
    Function to convert a symbolic vector to a numerical vector (numpy array)
    """
    if simsubs2 is None and simsubs3 is None:
        num_var = np.asarray(flatten(np.array(sim_vec.subs(
            simsubs, numsubs)).astype(np.float64)))
    elif simsubs3 is None:
        num_var = np.asarray(flatten(np.array(sim_vec.subs(
            simsubs, numsubs).subs(simsubs2, numsubs2)).astype(np.float64)))
    else:
        num_var = np.asarray(flatten(np.array(sim_vec.subs(
            simsubs, numsubs).subs(simsubs2, numsubs2).subs(simsubs3, numsubs3)).astype(np.float64)))
    return num_var

################################################################################
def vector_sym(sym1=symbols('a1'), sym2=symbols('a2'), sym3=symbols('a3')):
    """
    Function to symbolically display a 3x1 vector
    """
    a1, a2, a3 = symbols('a1 a2 a3')
    Vector = Matrix([a1, a2, a3])
    return Vector.subs(a1, sym1).subs(a2, sym2).subs(a3, sym3)

################################################################################
def crossProduct(vect_A, vect_B):
    """
    Function to find cross product symbolically of 2 3x1 vectors. 
    Made out of desperation because cross(a,b) did not find with vector_sym function
    """
    cross_P = []
    cross_P.append(vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1])
    cross_P.append(vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2])
    cross_P.append(vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0])
    return Matrix(cross_P)

################################################################################
def dotProduct(vect_A, vect_B, n=3):
    """
    Function to find dot product symbolically. 
    Made out of desperation because dot(a,b) did not find with vector_sym function
    """
    product = 0
    # Loop for calculate dot product
    for i in range(0, n):
        product = product + vect_A[i] * vect_B[i]

    return product

################################################################################
def skew_sym(sym1=symbols('a1'), sym2=symbols('a2'), sym3=symbols('a3')):
    """
    Function to symbolically calculate a 3x3 skew-symmetric matrix
    """
    a1, a2, a3 = symbols('a1 a2 a3')
    Skew = Matrix([
        [0, -a3, a2],
        [a3, 0, -a1],
        [-a2, a1, 0]
    ])
    return Skew.subs(a1, sym1).subs(a2, sym2).subs(a3, sym3)

# Test the skew function
#b1, b2, b3 = symbols('b1 b2 b3')
#skew_sym(b1, b2, b3)
#skew_sym()

################################################################################
def get_matrix_terms(matrix, dotx, doty, dotz):
    """
    Expand the terms of a 3x1 matrix into a 3x3 matrix(with the coefficients).
    Example use: When deriving the kinematic differential equation get the coefficients of the 
                 3x1 matrix containing the angular velocity components (dotx, doty, dotz)
    """
    return Matrix([
        [matrix[0].coeff(dotx), matrix[0].coeff(doty), matrix[0].coeff(dotz)],
        [matrix[1].coeff(dotx), matrix[1].coeff(doty), matrix[1].coeff(dotz)],
        [matrix[2].coeff(dotx), matrix[2].coeff(doty), matrix[2].coeff(dotz)]])

################################################################################
#%% Elementary rotation matrices Functions:
################################################################################

def C1(angle=symbols("theta_1")):
    x = symbols('x')
    Rx = Matrix([
        [1, 0, 0],
        [0, cos(x), sin(x)],
        [0, -sin(x), cos(x)]])
    return Rx.subs(x, angle)


def C2(angle=symbols("theta_2")):
    y = symbols('y')
    Ry = Matrix([
        [cos(y), 0, -sin(y)],
        [0,  1, 0],
        [sin(y), 0, cos(y)]])
    return Ry.subs(y, angle)


def C3(angle=symbols("theta_3")):
    z = symbols('z')
    Rz = Matrix([
        [cos(z), sin(z), 0],
        [-sin(z),  cos(z), 0],
        [0,    0, 1]])
    return Rz.subs(z, angle)

################################################################################
def DCM(mode, rot3, rot2=None, rot1=None, Eul_ang=None, invorder=False):
    """
    Function to calculate the rotation matrix from the 3 angles of the Euler angles.
    Input:
        mode: 'sim' or 'num' <- sim for symbolic, num for numerical
        rot1, rot2, rot3: Rotation angles
        Eul_ang: Euler angles. 
        invorder: If True, the angles are in the inverse order.
    
        (!) Note the angles are in the order: theta3, theta2, theta1, where:
            theta1 -> rotation along the 1st rotation axis 
            theta2 -> rotation along the 2nd rotation axis 
            theta3 -> rotation along the 3rd rotation axis 

        if invorder == True: the angles are in the order: theta1, theta2, theta3, where:
            theta1 -> rotation along the 3rd rotation axis 
            theta2 -> rotation along the 2nd rotation axis 
            theta3 -> rotation along the 1st rotation axis 
    Output:
        R: Rotation matrix

    Example use: get_R_matrix('num', 3, 2, 1, np.array([-15, 25, 10]))
    """
    x, y, z = symbols('x y z')

    if invorder==False:
        if rot3 == 1:
            R3 = C1(x)
        elif rot3 == 2:
            R3 = C2(x)
        elif rot3 == 3:
            R3 = C3(x)

        if rot2 == 1:
            R2 = C1(y)
        elif rot2 == 2:
            R2 = C2(y)
        elif rot2 == 3:
            R2 = C3(y)

        if rot1 == 1:
            R1 = C1(z)
        elif rot1 == 2:
            R1 = C2(z)
        elif rot1 == 3:
            R1 = C3(z)

    if invorder == True:
        if rot3 == 1:
            R3 = C1(z)
        elif rot3 == 2:
            R3 = C2(z)
        elif rot3 == 3:
            R3 = C3(z)

        if rot2 == 1:
            R2 = C1(y)
        elif rot2 == 2:
            R2 = C2(y)
        elif rot2 == 3:
            R2 = C3(y)

        if rot1 == 1:
            R1 = C1(x)
        elif rot1 == 2:
            R1 = C2(x)
        elif rot1 == 3:
            R1 = C3(x)

    if rot2 == None:
        R2 = 1
    if rot1 == None:
        R1 = 1

    R = R1 * R2 * R3
   
    if mode == 'num':
        if rot2 != None and rot1 == None:
            R = R.subs(z, Eul_ang[0])
            R = R.subs(y, Eul_ang[1])
        if rot2 == None and rot1 == None:
            R = R.subs(z, Eul_ang[0])

        R = R.subs(x, Eul_ang[0])
        if rot2 != None:
            R = R.subs(y, Eul_ang[1])
        if rot1 != None:
            R = R.subs(z, Eul_ang[2])
        R = np.array(R).astype(np.float64)
    return R

################################################################################
def diff_kinem_Euler_theta_matrix(rot1, rot2, rot3, invorder=False):
    """
    Differential kinematics matrix for any Euler angle sequence.

    If the Euler kinematics differential equation is:
        d(theta) = R(theta) * d(omega) 
    [3x1 matrix] = [3x3 matrix] * [3x1 matrix]
    
    The function returns the 3x3 R(theta) matrix of the differential equation in
    symbollical form.

    This can be then used together with the diff_kinem_Euler function from the 
    numerical library to solve for any Euler kinematic differential equation.

    Input:
        rot1, rot2, rot3: Rotation angles.
        invorder: If True, the angles are in the inverse order.
    
        (!) Note the angles are in the order: theta3, theta2, theta1, where:
            theta1 -> rotation along the 1st rotation axis 
            theta2 -> rotation along the 2nd rotation axis 
            theta3 -> rotation along the 3rd rotation axis 

        if invorder == True: the angles are in the order: theta1, theta2, theta3, where:
            theta1 -> rotation along the 3rd rotation axis 
            theta2 -> rotation along the 2nd rotation axis 
            theta3 -> rotation along the 1st rotation axis 
    """

    x, y, z = symbols("theta_1 theta_2 theta_3")
    dotx, doty, dotz = symbols(r"\dot{\theta_1} \dot{\theta_2} \dot{\theta_3}")

    if invorder == False:
        if rot3 == 1:
            R3 = C1(z)
            column1 = Matrix([dotz, 0, 0])
        elif rot3 == 2:
            R3 = C2(z)
            column1 = Matrix([0, dotz, 0])
        elif rot3 == 3:
            R3 = C3(z)
            column1 = Matrix([0, 0, dotz])

        if rot2 == 1:
            R2 = C1(y)
            column2 = Matrix([doty, 0, 0])
        elif rot2 == 2:
            R2 = C2(y)
            column2 = Matrix([0, doty, 0])
        elif rot2 == 3:
            R2 = C3(y)
            column1 = Matrix([0, 0, doty])

        if rot1 == 1:
            R1 = C1(x)
            column3 = Matrix([dotx, 0, 0])
        elif rot1 == 2:
            R1 = C2(x)
            column3 = Matrix([0, dotx, 0])
        elif rot1 == 3:
            R1 = C3(x)
            column3 = Matrix([0, 0, dotx])

    if invorder==True:
        if rot3 == 1:
            R3 = C1(x)
            column1 = Matrix([dotx, 0, 0])
        elif rot3 == 2:
            R3 = C2(x)
            column1 = Matrix([0, dotx, 0])
        elif rot3 == 3:
            R3 = C3(x)
            column1 = Matrix([0, 0, dotx])

        if rot2 == 1:
            R2 = C1(y)
            column2 = Matrix([doty, 0, 0])
        elif rot2 == 2:
            R2 = C2(y)
            column2 = Matrix([0, doty, 0])
        elif rot2 == 3:
            R2 = C3(y)
            column1 = Matrix([0, 0, doty])

        if rot1 == 1:
            R1 = C1(z)
            column3 = Matrix([dotz, 0, 0])
        elif rot1 == 2:
            R1 = C2(z)
            column3 = Matrix([0, dotz, 0])
        elif rot1 == 3:
            R1 = C3(z)
            column3 = Matrix([0, 0, dotz])

    Rmatrix = simplify(column1 + R3*column2 + R3*R2*column3)
    Rmatrix = get_matrix_terms(Rmatrix, dotx, doty, dotz)
    Rmatrix = simplify(Rmatrix.inv())

    return Rmatrix

################################################################################
def diff_kinem_Euler(rot1, rot2, rot3, w=None, invorder=True):
    """
    Differential kinematics for any Euler angle sequence.

    Input:
        rot1, rot2, rot3: Rotation angles.
        w: Angular velocity vector as a sympy Matrix. i.e. Matrix([w1, w2, w3])

    (!) Note: Need to review the solver. It is not working properly.
    """
    ## Required library
    from sympy.solvers.ode.systems import dsolve_system

    ## Angular velocity vector
    time = Symbol('t')
    if w ==None:
        w = Matrix([sin(0.1*time), 0.01, cos(0.1*time)]) * np.deg2rad(5)

    ## Thetas matrix
    thetas = diff_kinem_Euler_theta_matrix(rot1, rot2, rot3, invorder=invorder)

    ## Euler angle rates
    dotx, doty, dotz = symbols(r"\dot{\theta_1} \dot{\theta_2} \dot{\theta_3}")
    dot_angles = Matrix([dotx, doty, dotz])

    ## Differential Equation
    from sympy.physics.vector import dynamicsymbols
    xt, yt, zt = dynamicsymbols("theta_1 theta_2 theta_3")
    thetanum = np.deg2rad([80, 30, 40])

    eq1 = Eq(diff(xt, time), (thetas @ w)[0])
    eq2 = Eq(diff(yt, time), (thetas @ w)[1])
    eq3 = Eq(diff(zt, time), (thetas @ w)[2])

    eqs = Matrix([eq1, eq2, eq3])

    ## Solve the differential equation
    sol = dsolve_system(eqs, t=time, ics={xt: thetanum[0], yt: thetanum[1], zt: thetanum[2]})
    sol = Matrix(sol[0])
    return sol
################################################################################
#%% Rigid Body Dynamics Functions:
################################################################################

def char_poly(J):
    """
    Inputs a matrix in sympy form.
    Finds the characteristic polynomial of a matrix in sympy form.
    Takes the coefficients of the polynomial as a numpy array and Outputs the roots of the polynomial.
    NOTE: TBH I could also just use the find_eigen function or just numpy... But this way I get the characteristic polynomial.
    """
    # J.charpoly() gives the characteristic polynomial of J. Can also write as J.charpoly().as_expr() to just get the poly equation
    char_eq = J.charpoly()
    coef = np.array(char_eq.all_coeffs())
    return J.charpoly().as_expr(), np.roots(coef)

################################################################################
def find_eigen(J):
    """
    Input: a matrix in sympy form.
    Output: the eigenvalues and eigenvectors of a matrix in sympy form as numpy arrays
    """
    # J.eigenvects() gives the eigenvectors of J. Can also write as J.eigenvects().as_expr() to just get the eigenvectors
    Eigen = np.linalg.eigh(np.array(J, dtype='float'))
    Eigenvalues = Eigen[0]
    Eigenvectors = Eigen[1]
    return Eigenvalues, Eigenvectors

################################################################################
def Inertia_cylinder():
    m = Symbol('m', positive=True, real=True)
    r = Symbol('r', positive=True, real=True)
    h = Symbol('h', positive=True, real=True)
    I1 = 1/12 * m * (3*r**2 + h**2)
    I2 = I1
    I3 = 1/2 * m * r**2
    return I1, I2, I3

################################################################################
#%% Attitude Dynamics Functions:
################################################################################

def eq_torquefree_motion(a):
    """
    Equation of torque free motion for dual spin systems in body coordinates.

    Input:
        a: wheel spin axis
    """
    # Define the symbolic variables
    omega1sym, omega2sym, omega3sym = symbols(
        '\omega_1 \omega_2 \omega_3')
    dotomega1sym, dotomega2sym, dotomega3sym = symbols(
        '\dot{\omega}_1 \dot{\omega}_2 \dot{\omega}_3')
    I1_sym, I2_sym, I3_sym = symbols('I1, I2, I3')
    hs_sym = symbols('h_s')

    # wheel spin axis vector in body coordinates
    if a == 1:
        a_vecsym = Matrix([1, 0, 0])
    elif a == 2:
        a_vecsym = Matrix([0, 1, 0])
    elif a == 3:
        a_vecsym = Matrix([0, 0, 1])
    else:
        print('a must be 1, 2 or 3')

    # Define the matrices
    omega_vecsym = Matrix([omega1sym, omega2sym, omega3sym])
    I_vecsym = Matrix([[I1_sym, 0, 0], [0, I2_sym, 0], [0, 0, I3_sym]])
    dotomega_vecsym = Matrix([dotomega1sym, dotomega2sym, dotomega3sym])

    # Define the equation
    T_vecsym = I_vecsym * dotomega_vecsym + \
        omega_vecsym.cross(I_vecsym * omega_vecsym + a_vecsym*hs_sym)

    return T_vecsym


################################################################################
#%% Active Attitude Control Functions:
################################################################################

def close_loop_TF(form='freq', subs=False, Kd=None, Kp=None, I=None, xi=None, omega_n=None):
    """"
    Input:
        - form: string with the form of the loop
                "gains" -> outputs the TF with Kp, Kd, I as cosntants
                "freq" -> outputs the TF with w_n, xi  as constants
        - subs: boolean to substitute the constants in the TF
    Output:
        - TF(tf): transfer function for a close-loop system with PD control
    """
    if form == "gains":
        TF_num = (Symbol('K_p') / (Symbol('I')*Symbol('s')
                                   ** 2 + Symbol('K_d')*Symbol('s')))
        TF_den = (Symbol('I') + Symbol('K_p') / (Symbol('I') *
                                                 Symbol('s')**2 + Symbol('K_d')*Symbol('s')))
    elif form == "freq":
        omega_n_sym = Symbol('\omega_n')
        xi_sym = Symbol('xi')
        s_sym = Symbol('s')

        TF_num = omega_n_sym**2
        TF_den = s_sym**2 + 2*xi_sym*omega_n_sym*s_sym + omega_n_sym**2
    else:
        raise ValueError("form must be 'gains' or 'freq'")

    TF = TF_num / TF_den

    if subs is True and form == "gains" and Kd is not None and Kp is not None and I is not None and xi is None and omega_n is None:
        TF = TF.subs(Symbol('K_p'), Kp).subs(
            Symbol('K_d'), Kd).subs(Symbol('I'), I).simplify()
    elif subs is True and form == "freq" and Kd is None and Kp is None and I is None and xi is not None and omega_n is not None:
        TF = TF.subs(xi_sym, xi).subs(omega_n_sym, omega_n).simplify()

    return TF
    
################################################################################
def max_overshoot_CLTF(known_real_part='y', subs=False, xi=None, xiomega_n=None):
    """
    Maximum overshoot for a close-loop system with PD control.

    Input:
        - known_real_part(string):
            - 'y': real part of the system is known
            - 'n': real part of the system is unknown
    """
    xi_sym = symbols('xi')
    omega_nsym, omega_dsym = symbols('omega_n omega_d')

    if known_real_part == 'y':
        xiomega_n_sym = xi_sym*omega_nsym
        Mp = exp(-pi*xiomega_n_sym/omega_dsym)
        if subs is True and xiomega_n is not None:
            Mp = Mp.subs(xiomega_n_sym, xiomega_n).subs(pi, np.pi).simplify()

    elif known_real_part == 'n':
        Mp = exp(-pi*xi_sym/(1-xi_sym**2) ** 0.5)
        if subs is True and xi is not None:
            Mp = Mp.subs(Symbol('xi'), xi).subs(pi, np.pi).simplify()
    else:
        raise ValueError('known_real_part must be "y" or "n"')

    
    return Mp

################################################################################
def syms2tf(TF_sym, plot=True):
    """
    Get transfer function from equation. Plot pzmap, rlocus, stepresponse, impulseresponse.
    Input:
     - eq(sym): symbolic equation
    Output
     - TF(tf): transfer function
    """
    import control as ctl
    import matplotlib.pyplot as plt
    #Get num and den of Symbolic TF
    symNum, symDen = numer(TF_sym), denom(TF_sym)

    # Convert Symbolic num and den to polynomial
    TFnum = Poly(symNum, Symbol('s')).all_coeffs()
    TFden = Poly(symDen, Symbol('s')).all_coeffs()
    TFnum = np.array(TFnum).astype(float)
    TFden = np.array(TFden).astype(float)

    TF = ctl.tf(TFnum, TFden)

    poles, zeros = ctl.pzmap(TF, plot=True, grid=True)

    if plot is True:
        plt.show()

        ctl.rlocus(TF)
        plt.show()

        T, yout = ctl.step_response(TF)  # , T=np.linspace(0, 5, 100))
        plt.plot(T, yout)
        plt.show()

        T, yout = ctl.impulse_response(TF)
        plt.plot(T, yout)
        plt.show()

    return TF, poles, zeros

################################################################################
#%% Non-spherical Earth Harmonics Functions:
################################################################################

def non_spherical_pot(n, m, subs=False, Jnmval=None, lamb_nm_val=None, h_orbit=None):
    """
    Returns the Potential and radial, North-South and East-Weest accelerations 
    based on the potential formulation for the gravity field of the Earth.
    
    Input:
        - n(int): degree of the spherical harmonic
        - m(int): order of the spherical harmonic
        - subs(bool): if True, substitute the following values for symbols, acc units in m/s^2. 
            - R_e: 6.371e6
            - G: 6.6743e-11
            - M: 5.972167867791379e+24
            - Jnm: Jnmval
            - lamb_nm: lamb_nm_val
            - h_orbit: h_orbit (input as orbital altitude, not radius)
    Output:
        - Pn: polynomial for degree n
        - Pnm: polynomial for degree n, order m (@ x= sin(phi) )
        - Unm(sym): potential formulation for the spherical harmonic
        - a_r, a_phi, a_lamb: radial, North-South and East-Weest accelerations
    """
    # Define symbols:
    from sympy.abc import x, G, M, r
    phi, R_e, Jnm, lamb, lamb_nm = symbols('phi R_e J_{nm} lambda \lambda_{nm}')

    # Calculate the polynomial and the potential:
    Pnderiv = (1-x**2)**n
    Pn = (1/((-2)**n*np.math.factorial(n)) * (Pnderiv.diff((x, n)))).simplify()
    Pnm = ((1-x**2)**(m/2)*(Pn.diff((x, m)))).subs(x, sin(phi)).simplify()

    # Calculate the potential
    Unm = G*M/r * (R_e/r)**n * Jnm * Pnm * cos(m*(lamb - lamb_nm))
    
    # Calculate the accelerations:
    a_r = -diff(Unm, r)
    a_phi = - 1/r * diff(Unm, phi)
    a_lamb = - 1/(r*cos(phi)) * diff(Unm, lamb)

    if subs is True:
        a_r = a_r.subs(R_e, 6.371e6).subs(
            G, 6.6743e-11).subs(r, h_orbit+6.371e6).subs(M, 5.972167867791379e+24).subs(Jnm, Jnmval).subs(lamb_nm, lamb_nm_val).simplify()
        a_phi = a_phi.subs(R_e, 6.371e6).subs(
            G, 6.6743e-11).subs(r, h_orbit+6.371e6).subs(M, 5.972167867791379e+24).subs(Jnm, Jnmval).subs(lamb_nm, lamb_nm_val).simplify()
        a_lamb = a_lamb.subs(R_e, 6.371e6).subs(
            G, 6.6743e-11).subs(r, h_orbit+6.371e6).subs(M, 5.972167867791379e+24).subs(Jnm, Jnmval).subs(lamb_nm, lamb_nm_val).simplify()

    return Pn, Pnm, Unm, a_r, a_phi, a_lamb
