from sympy import *
import numpy as np
from sympy.vector import CoordSys3D
from sympy import lambdify
import matplotlib.pyplot as plt
from IPython.display import display
from scipy.interpolate import interp1d
from scipy.integrate import quadrature

def integrand(x, f1, f2):

    return f1(x) + f2(x)


# Sympy Variables
x, a, b, c, d, e, t, c_a = symbols('x a b c d e t c_a')
N = CoordSys3D('N')

# Airfoil Definitions
airfoil = x*N.i + 5*t*(a*sqrt(c_a) + b*x + c*Pow(x, 2) + d*Pow(x, 3) + e*Pow(x, 4))*N.j
airfoil = airfoil.subs(x, x/c_a)

# Airfoil Tangent
tangent = diff(airfoil, x)

# Airfoil Normal
normal = -1*tangent.coeff(N.j)*N.i + tangent.coeff(N.i)*N.j
unit_normal = normal/normal.magnitude()

# Import Test data
p_data = 101325+np.load('pressure_data/v2_12.npy')

# Tap Locations
x_tap = np.array(
    [0.0, 0.210, 0.526, 0.842, 1.158, 1.474, 1.790, 2.422, 2.738, 3.054, 3.686, 4.002, 4.318, 4.634, 4.950, 5.265]
)
x_tap *= 0.0254

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

# Pressure Interpolation Function
p_interp_u = interp1d(x_tap, p_data[:, 0], fill_value='extrapolate')
p_interp_l = interp1d(x_tap, p_data[:, 1], fill_value='extrapolate')


x_c = np.linspace(0, 6, 100)*0.0254

# Define the integrand
ds = sqrt(1 + Pow(tangent.coeff(N.j), 2))*unit_normal.coeff(N.j)
ds_f = lambdify(x, ds)
Iu = lambda t: integrand(t, ds_f, p_interp_u)
Il = lambda t: integrand(t, ds_f, p_interp_l)

Lu = quadrature(Iu, 0, 6*0.0254)[0]
Ll = quadrature(Il, 0, 6*0.0254)[0]

L = Ll - Lu

print(L)


plt.plot(x_c, f_p(x_c))
plt.show()

