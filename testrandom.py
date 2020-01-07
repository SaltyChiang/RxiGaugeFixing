import numpy as np
import matplotlib.pyplot as plt

volume = 16**3 * 16
data = np.fromfile('data/random.data').reshape(volume, 8)
data0 = np.around(data[:, 0] * 1000).astype('int')
print(data0.shape)
xi = 1.0
x = np.arange(-5000, 5001) / 1000
y1 = np.exp(-0.5 * np.power(x / xi, 2)) / (np.math.sqrt(2 * np.math.pi) * xi)
y2 = np.zeros(10001)
for i in range(volume):
    if data0[i] >= -5000 and data0[i] < 5001:
        y2[data0[i] + 5000] += 1
y2 /= volume / 1000
plt.plot(x, y1, 'b,-')
plt.bar(x, y2, 1 / 1000)
plt.show()
