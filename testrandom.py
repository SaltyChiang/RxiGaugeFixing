import numpy as np
from numba import jit
import matplotlib.pyplot as plt

volume = 16**3 * 16
sigma2 = 1.0
sigma_range = int(7)
density = int(50)

left = -sigma_range * density
right = sigma_range * density + int(1)
data = np.fromfile('data/random.data').reshape(volume, 8)
data0 = np.around(data[:, 0] * density).astype('int')
x = np.arange(left, right) / density
y1 = np.exp((-0.5 / sigma2) * np.power(x, 2)) / \
     np.math.sqrt(2 * np.math.pi * sigma2)
y2 = np.zeros(right - left)


def fill(sample, y):
    for i in range(sample.size):
        if sample[i] >= left and sample[i] < right:
            y[sample[i] - left] += 1
    y /= sample.size / density


fill(data0, y2)
plt.plot(x, y1, 'b,-')
plt.bar(x, y2, 1 / density)
plt.show()
