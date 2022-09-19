import numpy as np
import matplotlib.pyplot as plt

# Speeds
vs = ['v1', 'v2']
AoAs = [0, 4, 8, 12, 16]

# Tap Locations
x_tap = np.array(
    [0.0, 0.210, 0.526, 0.842, 1.158, 1.474, 1.790, 2.422, 2.738, 3.054, 3.686, 4.002, 4.318, 4.634, 4.950, 5.265]
)
x_tap /= 6.0

# Manometer Data
delHi = np.array([0.1, 0.15])*0.0254
delHi = dict(zip(vs, delHi))
delHo = np.array([0.84-0.28, 1.05-0.28])*0.0254
delHo = dict(zip(vs, delHo))

# Imports
data_exp = {
    i: {
        j: {} for j in AoAs
    } for i in vs
}

for vel in vs:
    for a in AoAs:
        data_exp[vel][a] = (np.load('pressure_data/{0}_{1}.npy'.format(vel, a)) + 1000*9.81*delHo[vel])/(1000*9.81*delHi[vel])

# Sim Data
sim = np.loadtxt('sim_data/v1_0.csv', delimiter=',', skiprows=1, usecols=(0, 1))
dmin, dmax = min(sim[:, 1]), max(sim[:, 1])
dcord = dmax - dmin
sim[:, 1] = (sim[:, 1] - dmin)/(dmax - dmin)


plt.plot(x_tap, data_exp["v1"][16])
#plt.plot(sim[:, 1], sim[:, 0])
plt.show()

# Plotting
#fig, axs = plt.subplots(5, 2, figsize=(6.5, 9), sharex='col', sharey='row')

#for j, vel in enumerate(vs):
#    for i, a in enumerate(AoAs):
#        print(i, j)
#        axs[i, j].plot(x_tap, data_exp[vel][a])

#plt.show()
