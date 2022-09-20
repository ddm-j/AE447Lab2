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

    return d*dA*u


def force(x_c, p, af, limits, type='lift', method='simps'):

    N = CoordSys3D('N')

    # Pressure Interpolation Function
    int_u = interp1d(x_c, p[:, 0], fill_value='extrapolate')
    int_l = interp1d(x_c, p[:, 1], fill_value='extrapolate')

    # Differential Unit Area
    dA = 24*0.0254*np.linalg.norm(af[1:] - af[:-1], axis=1)

    # Normal Vector
    t = af[1:] - af[:-1]
    n = np.zeros(t.shape)
    n[:, 0] = t[:, 1]
    n[:, 1] = -t[:, 0]
    n /= np.linalg.norm(n, axis=1)[:, None]

    # Vector Area
    dA_v = n*dA[:, None]

    # Pressure Interpolation Points
    x_int = af[1:, 0] + (af[1:, 0] - af[:-1, 0])/2

    # Upper Surface Forces
    F_upper = dA_v*int_u(x_int)[:, None]
    F_upper[1:, 1] *= -1

    # Lower Surface Forces
    F_lower = dA_v*int_l(x_int)[:, None]
    #F_lower[:, 1] *= -1

    # Total Forces
    F = F_upper + F_lower

    #plt.plot(af[:-1, 0], F[:, 1])
    #plt.show()

    F = np.sum(F, axis=0)

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

# Discretize the Space
x_l = np.geomspace(0.00000001, coeffs['c_a'], num=1000)
x_l[0] = 0.0


# Get Airfoil Coordinates
a_f = lambdify(x, airfoil.coeff(N.j))
y_l = a_f(x_l)

# Create Differential Arc Length
dx = np.diff(x_l)
dy = np.diff(y_l)
ds = np.sqrt(dx**2 + dy**2)
S_l = np.insert(np.cumsum(ds), 0, 0.0)
S = S_l[-1]

# Get x values for equally spaced arc lengths
S_int = interp1d(S_l, x_l)
s_l = np.linspace(0.0, S, 1000)
s_x = S_int(s_l)

af = np.array(list(zip(s_x, a_f(s_x))))


# Tap Locations
x_tap = np.array(
    [0.0, 0.210, 0.526, 0.842, 1.158, 1.474, 1.790, 2.422, 2.738, 3.054, 3.686, 4.002, 4.318, 4.634, 4.950, 5.265]
)
x_tap *= 0.0254

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

for file in files:

    # Import Test data
    alpha = file[-1]
    p_data = -1*np.load('pressure_data/{0}_{1}.npy'.format(file[0], alpha))

    # Calculate Aerodynamic Forces
    Fs = force(x_tap, p_data, af, (0.0, coeffs['c_a']), type='lift')
    alpha = np.deg2rad(alpha)
    R = np.array([
        [np.cos(alpha), -np.sin(alpha)],
        [np.sin(alpha), np.cos(alpha)]
    ])
    F = np.matmul(R, Fs)
    L = F[1]/(q[file[0]]*SA)
    D = F[0]/(q[file[0]]*SA)

    # Store Results
    idx = AoAs.index(file[-1])
    coeff_data[file[0]][idx, 1] = L
    coeff_data[file[0]][idx, 2] = D

plt.plot(coeff_data['v1'][:, 0], coeff_data['v1'][:, 1])
plt.twinx()
plt.plot(coeff_data['v1'][:, 0], coeff_data['v1'][:, 2])
plt.show()
