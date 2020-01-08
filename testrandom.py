import numpy as np
from numba import jit
import matplotlib.pyplot as plt

volume = 16**3 * 16
data = np.fromfile('data/random.data').reshape(volume, 8)
data0 = np.around(data[:, 0] * 1000).astype('int')
sigma2 = 1.0
x = np.arange(-5000, 5001) / 1000
y1 = np.exp((-0.5 / sigma2) * np.power(x, 2)) / \
     np.math.sqrt(2 * np.math.pi * sigma2)
y2 = np.zeros(10001)
print(np.size(data0, 0))


@jit
def fill(sample, y):
    for i in range(sample.size):
        if sample[i] >= -5000 and sample[i] < 5001:
            y[sample[i] + 5000] += 1
    y /= sample.size / 1000


fill(data0, y2)
plt.plot(x, y1, 'b,-')
plt.bar(x, y2, 1 / 1000)
plt.show()
