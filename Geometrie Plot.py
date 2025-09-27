import pandas as pd

df = pd.read_csv('data/koordinatenpunkte Kreisel.py')

import matplotlib.pyplot as plt
import numpy as np

n = 100

fig = plt.figure(figsize=(12,6))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122,projection='3d')
y = df.values[:, 0]
x = df.values[:, 1]
t = np.linspace(0, np.pi*2, n)

xn = np.outer(x, np.cos(t))
yn = np.outer(x, np.sin(t))
zn = np.zeros_like(xn)

for i in range(len(x)):
    zn[i:i+1,:] = np.full_like(zn[0,:], y[i])

ax1.plot(x, y)
ax2.plot_surface(xn, yn, zn)
plt.show()