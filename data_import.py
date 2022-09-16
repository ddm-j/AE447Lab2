import glob
import numpy as np

# Get Paths of All Data Files
file_names = glob.glob('raw_data/*')
file_names = [i.split('\\')[-1] for i in file_names]

# File Mapping Dictionary
map = {
    'High': 'v1',
    'Low': 'v2'
}

# Load Our Files
# Create Empty Placeholders
AoAs = [0, 4, 8, 12, 16]
data = {
    'v1': {
        i: {} for i in AoAs
    },
    'v2': {
        i: {} for i in AoAs
    }
}

# Load Raw Data Files
for fname in file_names:
    AoA = fname.split('_')[1][:-7]
    velocity = map[fname.split('_')[0]]
    if AoA[0] == '-':
        surface = 'lower'
        AoA = int(AoA[1:])
    else:
        surface = 'upper'
        AoA = int(AoA)

    # Iterate Through the File
    p_data = np.zeros((16, 50))
    with open('raw_data/{0}'.format(fname), 'r') as file:
        for line in file:
            if '#' in line:
                c_idx = int(line.split(' ')[-1])
                continue
            else:
                r_idx = int(line.split(' ')[0])
                val = float(line.split(' ')[1])

            p_data[r_idx, c_idx] = val

    # Compute Row-wise Mean of the Pressure Samples
    p_avg = np.mean(p_data, axis=1)

    # Store This Data
    data[velocity][AoA].update({'{0}'.format(surface): p_avg})

# Consolidate Data & Write Out
for velocity in data.keys():
    for AoA in AoAs:
        p_data = np.zeros((16, 2))
        p_data[:, 0] = data[velocity][AoA]['upper']
        if AoA != 0:
            p_data[:, 1] = data[velocity][AoA]['lower']
        else:
            p_data[:, 1] = data[velocity][AoA]['upper']

        # Save Data to Binary NPY File
        np.save('pressure_data/{0}_{1}.npy'.format(velocity, AoA), p_data)