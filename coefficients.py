from sympy import *
import numpy as np
from sympy.vector import CoordSys3D
from sympy import lambdify
import matplotlib.pyplot as plt
import matplotlib
from IPython.display import display
from scipy.interpolate import interp1d
from scipy.integrate import quadrature, simps

def integrand(x, f1, f2):

    return f1(x)*f2(x)

def force(x_c, p, tangent, normal, limits, type='lift', method='simps'):
    N = CoordSys3D('N')

    # Pressure Interpolation Function
    int_u = interp1d(x_c, p[:, 0], fill_value='extrapolate')
    int_l = interp1d(x_c, p[:, 1], fill_value='extrapolate')

    # Define the integrand
    d = 1
    if type == 'lift':
        vec = N.j
        d = -1
    elif type == 'drag':
        vec = N.i
    else:
        raise ValueError('Unreocgnized force type.')

    # Differential Unit Area
    dA = sqrt(1 + Pow(tangent.coeff(N.j), 2)) * normal.coeff(vec)
    dA_f = lambdify(x, dA)

    # Interpolated Lambda Functions
    Iu = lambda t: integrand(t, dA_f, int_u)
    Il = lambda t: integrand(t, dA_f, int_l)

    #x_n = np.linspace(0.0000001, 6*0.0254, 200)
    #plt.plot(x_n, Iu(x_n))
    #plt.show()


    if method == 'quad':
        Fu = quadrature(Iu, limits[0], limits[1])[0]
        Fl = quadrature(Il, limits[0], limits[1])[0]

    elif method == 'simps':
        x_s = np.linspace(limits[0], limits[1], 20000)
        Fu = simps(Iu(x_s), x_s)
        Fl = simps(Il(x_s), x_s)
    else:
        raise ValueError('Integration Method not supported.')

    F = Fl - Fu if type == 'lift' else Fu - Fl

    return F


# Sympy Variables
x, a, b, c, d, e, t, c_a = symbols('x a b c d e t c_a')
N = CoordSys3D('N')

# Airfoil Definitions
airfoil = x*N.i + c_a*5*t*(a*sqrt(x) + b*x + c*Pow(x, 2) + d*Pow(x, 3) + e*Pow(x, 4))*N.j
airfoil = airfoil.subs(x, x/c_a)

# Airfoil Tangent
tangent = diff(airfoil, x)

# Airfoil Normal
normal = -1*tangent.coeff(N.j)*N.i + tangent.coeff(N.i)*N.j
unit_normal = normal/normal.magnitude()

# Import Test data
alpha = 16
p_data = 101325+np.load('pressure_data/v2_{0}.npy'.format(alpha))

# Tap Locations
x_tap = np.array(
    [0.0, 0.210, 0.526, 0.842, 1.158, 1.474, 1.790, 2.422, 2.738, 3.054, 3.686, 4.002, 4.318, 4.634, 4.950, 5.265]
)
x_tap *= 0.0254

#plt.plot(x_tap, p_data[:, 0])
#plt.plot(x_tap, p_data[:, 1])
#plt.show()

# Airfoil Function Constants
coeffs = {
    'a': 0.2969,
    'b': -0.1260,
    'c': -0.3516,
    'd': 0.2843,
    'e': -0.1015,
    't': 0.12,
    'c_a': 6*0.0254
}
airfoil = airfoil.subs(coeffs)
unit_normal = unit_normal.subs(coeffs)
tangent = tangent.subs(coeffs)

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

# Calculate Aerodynamic Forces
Ls = force(x_tap, p_data, tangent, unit_normal, (0.00001, coeffs['c_a']), type='lift')
Ds = force(x_tap, p_data, tangent, unit_normal, (0.00001, coeffs['c_a']), type='drag')
Fs = np.array([Ls, Ds])
alpha = np.deg2rad(alpha)
R = np.array([
    [np.cos(alpha), np.sin(alpha)],
    [-np.sin(alpha), np.cos(alpha)]
])
F = np.matmul(R, Fs)
L = F[0]
D = F[1]
ar = L/D
print("Lift Force: {0}N per unit span".format(round(F[0], 2)))
print("Drag Force: {0}N per unit span".format(round(F[1], 2)))
print("Aerodynamic Ratio: {0}".format(round(ar, 2)))


