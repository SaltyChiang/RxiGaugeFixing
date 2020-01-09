import numpy as np

volume = 64 * 24**3
data_in = np.fromfile('data/rbc_conf_2464_m0.005_0.04_000495_hyp',
                      dtype='>f8').astype('float64').reshape(
                          4, 3, 3, 2, volume)
data_out = np.zeros((volume, 4, 3, 3, 2))

for c2 in range(3):
    for c1 in range(3):
        for l in range(volume):
            data_out[l, :, c1, c2, :] = data_in[:, c2, c1, :, l]
data_out.tofile('data/rbc_conf_2464_m0.005_0.04_000495_hyp_rearange')
