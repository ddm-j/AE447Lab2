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
data_sim = {
    i: {
        j: {} for j in AoAs
    } for i in vs
}

for vel in vs:
    for a in AoAs:
        data_exp[vel][a] = (np.load('pressure_data/{0}_{1}.npy'.format(vel, a)) + 1000*9.81*delHo[vel])/(1000*9.81*delHi[vel])
        data_sim[vel][a] = np.loadtxt('sim_data/{0}_{1}_xflr.csv'.format(vel, a), delimiter=',', skiprows=1, usecols=(0, 1))


# Plotting
fig, axs = plt.subplots(5, 2, figsize=(6.5, 9), sharex='col', sharey='row', layout='constrained')

# Titles
Res = [5.8818e+4, 7.2037e+4]
cols = ['Re = {:.2e}'.format(col) for col in Res]
rows = ['{0}'.format(col)+r'$^{\circ}$' for col in AoAs]
fig.suptitle('Pressure Coefficient', fontsize=11)

# Styles
from matplotlib import font_manager

font_path = 'EBGaramond-Regular.ttf'  # Your font path goes here
font_manager.fontManager.addfont(font_path)
prop = font_manager.FontProperties(fname=font_path)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = prop.get_name()

for ax, col in zip(axs[0], cols):
    ax.set_title(col, fontsize=10)

for ax, row in zip(axs[:, 0], rows):
    ax.set_ylabel(row, fontsize=10)
    ax.yaxis.set_label_coords(-0.15, 0.5)

for j, vel in enumerate(vs):
    for i, a in enumerate(AoAs):

        # New Scales for Axis Reversal
        miny = min(np.min(data_exp[vel][a]), np.min(data_sim[vel][a][:, 1]))
        maxy = max(np.max(data_exp[vel][a]), np.max(data_sim[vel][a][:, 1]))

        axs[i, j].plot(x_tap, data_exp[vel][a],
                       linewidth=1,
                       color='black')

        axs[i, j].plot(data_sim[vel][a][:, 0], data_sim[vel][a][:, 1],
                       linestyle='dotted',
                       linewidth=1.5,
                       color='#A9A9A9')
        axs[i, j].set_ylim(1.2*maxy, 1.2*miny)

        # Set x labels
        if a == 16:
            axs[i, j].set_xlabel(r'$\frac{x}{c}$', fontsize=10)

plt.savefig('charts/pressure_coefficients.png', dpi=400)
plt.show()
