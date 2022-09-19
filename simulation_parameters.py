import numpy as np

# Pressure Calculation
h_lc = 1189.0
P = 101325*(1 - 2.25577e-5 * h_lc)**5.25588
print("Atmospheric Pressure: {0} Pa".format(round(P, 2)))

# Density Calcuation
T_0 = 273.15
T = 27.0
T += T_0
R = 287.052874
rho = P/(R*T)
print("Atmospheric Density: {0} kg/m^3".format(round(rho, 2)))

# Relative Pressure
p_rho = P/rho
print('Pressure over Density: {0}'.format(round(p_rho, 2)))

# Viscosity Calculation
mu_0 = 1.716e-5
S_mu = 111.0
mu = mu_0*((T/T_0)**(3.0/2))*((T_0+S_mu)/(T+S_mu))
nu = mu/rho
print(nu)
print('Dynamic Viscosity: {:.4e} kg/m*s'.format(mu))
print('Kinematic Viscosity: {:.4e} Pa*s'.format(nu))


# Flow Velocity
delH = np.array([0.1, 0.15])
q = delH*0.0254*1000*9.81
U = np.sqrt(2*q/rho)
print('Flow velocities: [{0} {1}] m/s'.format(round(U[0], 3), round(U[1], 3)))

# Reynolds Number
c = 6*0.0254
Re = rho*U*c/mu
print('Reynolds: [{:.4e} {:0.4e}] Re'.format(Re[0], Re[1]))

# Turbulent Viscosity Ratio
factor = 0.1
nut = factor*nu
nut_t = 3*nut
print('nut: {:.4e}'.format(nut))
print('nuTilda: {:.4e}'.format(nut_t))