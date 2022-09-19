from sympy import *
import numpy as np
from sympy.vector import CoordSys3D
from sympy import lambdify
import matplotlib.pyplot as plt
import matplotlib
from IPython.display import display
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
from scipy.integrate import quadrature, simps, quad
import itertools


def integrand(x, f1, f2, sign=1):

    return sign*f1(x)*f2(x)


def differential_area(x, d_af, d='drag', c=0.0254*6, S=0.0254*24, surface='upper'):

    tanx = 1/c
    tany = d_af(x)

    ds = np.sqrt(1 + tany**2)
    dA = ds*S

    nx = tany
    ny = -tanx

    norm = np.sqrt(nx**2 + ny**2)

    ux = nx/norm
    uy = ny/norm

    u = ux if d == 'drag' else uy
    d = -1 if d == 'lift' and surface == 'lower' else 1

    return u*d*dA


def force(x_c, p, d_a_int, limits, type='lift', method='simps'):

    N = CoordSys3D('N')

    # Pressure Interpolation Function
    int_u = interp1d(x_c, p[:, 0], fill_value='extrapolate')
    int_l = interp1d(x_c, p[:, 1], fill_value='extrapolate')

    # Define the integrand
    if type == 'lift':
        vec = N.j
        d = -1
    elif type == 'drag':
        vec = N.i
        d = 1
    else:
        raise ValueError('Unreocgnized force type.')

    # Differential Unit Area
    dA_u = lambda t: differential_area(t, d_a_int, d=type, surface='upper')
    dA_l = lambda t: differential_area(t, d_a_int, d=type, surface='lower')

    # Interpolated Lambda Functions
    I = lambda t: integrand(t, dA_u, int_u) + integrand(t, dA_l, int_l)

    if method == 'quad':
        F = quad(I, limits[0], limits[1])[0]
        #Fl = quadrature(Il, limits[0], limits[1])[0]

    elif method == 'simps':
        x_s = np.linspace(limits[0], limits[1], 20000)
        F = simps(I(x_s), x_s)
        #Fl = simps(Il(x_s), x_s)
    else:
        raise ValueError('Integration Method not supported.')

    #F = Fl + Fu

    return F


# Sympy Variables
x, a, b, c, d, e, t, c_a = symbols('x a b c d e t c_a')
N = CoordSys3D('N')

# Airfoil Definitions
airfoil = x*N.i + c_a*5*t*(a*sqrt(x) + b*x + c*Pow(x, 2) + d*Pow(x, 3) + e*Pow(x, 4))*N.j
airfoil = airfoil.subs(x, x/c_a)
coeffs = {
    'a':0.2969,
    'b':-0.1260,
    'c':-0.3516,
    'd':0.2843,
    'e':-0.1015,
    't':0.12,
    'c_a':6*0.0254
}
airfoil = airfoil.subs(coeffs)

# Log Discretize the Space
x_l = np.geomspace(0.0001, coeffs['c_a'], num=1000)
x_l[0] = 0.0

# Get Airfoil Coordinates
a_f = lambdify(x, airfoil.coeff(N.j))
a_l = a_f(x_l)

# Interpolate the Data to Remove the Tangent Discontinuity
a_int = InterpolatedUnivariateSpline(x_l, a_l)
d_a_int = a_int.derivative()

# Tap Locations
x_tap = np.array(
    [0.0, 0.210, 0.526, 0.842, 1.158, 1.474, 1.790, 2.422, 2.738, 3.054, 3.686, 4.002, 4.318, 4.634, 4.950, 5.265]
)
x_tap *= 0.0254

#plt.plot(x_tap, p_data[:, 0])
#plt.plot(x_tap, p_data[:, 1])
#plt.show()

#a_f = lambdify(x, airfoil.coeff(N.j))
#n_fx = lambdify(x, unit_normal.coeff(N.i))
#n_fy = lambdify(x, unit_normal.coeff(N.j))
#x_n = np.linspace(0.00001, 6*0.0254, 50)
#y_n = a_f(x_n)

#n_x = n_fx(x_n)
#n_y = n_fy(x_n)

#plt.plot(x_n, y_n)
#plt.quiver(x_n, y_n, n_x, n_y)
#plt.axis('equal')
#plt.show()

vs = ['v1', 'v2']
AoAs = [0, 4, 8, 12, 16]
files = list(itertools.product(vs, AoAs))

delH = np.array([0.1, 0.15])
delH *= 0.0254
q = dict(zip(vs, 1000*9.81*delH))

SA = 0.0254*6*0.0254*24

coeff_data = np.zeros((5, 3))
coeff_data[:, 0] = AoAs
coeff_data = {
    "v1": coeff_data,
    "v2": coeff_data
}
print(coeff_data)

for file in files:

    # Import Test data
    alpha = file[-1]
    p_data = 101325+np.load('pressure_data/{0}_{1}.npy'.format(file[0], alpha))

    # Calculate Aerodynamic Forces
    Ls = force(x_tap, p_data, a_int, (0.0, coeffs['c_a']), type='lift')
    Ds = force(x_tap, p_data, a_int, (0.0, coeffs['c_a']), type='drag')
    Fs = np.array([Ls, Ds])
    alpha = np.deg2rad(alpha)
    R = np.array([
        [np.cos(alpha), np.sin(alpha)],
        [-np.sin(alpha), np.cos(alpha)]
    ])
    F = np.matmul(R, Fs)
    L = F[0]/(q[file[0]]*SA)
    D = F[1]/(q[file[0]]*SA)

    # Store Results
    idx = AoAs.index(file[-1])
    coeff_data[file[0]][idx, 1] = L
    coeff_data[file[0]][idx, 2] = D

plt.plot(coeff_data['v2'][:, 0], coeff_data['v2'][:, 2])
plt.show()
