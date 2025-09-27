import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('data/koordinatenpunkte Kreisel.py') # Profil
#df = pd.read_csv('data/Koordinatenpunkte Würfel.py')

MTT = np.matrix([[sum(df.values[:, 1]**2 + df.values[:, 2]**2), sum(-df.values[:, 0]*df.values[:, 1]), sum(-df.values[:, 0]*df.values[:, 2])],
[sum(-df.values[:, 1]*df.values[:, 0]), sum(df.values[:, 0]**2 + df.values[:, 2]**2), sum(-df.values[:, 1]*df.values[:, 2])],
[sum(-df.values[:, 2]*df.values[:, 0]), sum(-df.values[:, 1]*df.values[:, 2]), sum(df.values[:, 0]**2+df.values[:, 1]**2)]])

print(MTT)

#Hauptinvarianten
HI1 = MTT[0, 0] + MTT[1, 1] + MTT[2, 2]
HI2 = MTT[0, 0]*MTT[1, 1] + MTT[2, 2]*MTT[1, 1] + MTT[0, 0]*MTT[2, 2] - MTT[0, 1]**2-MTT[1, 2]**2-MTT[0, 2]**2
HI3 = np.linalg.det(MTT)
print(f'Hauptträgheitsmomente: {HI1, HI2, HI3}')
n = 10
x = df.values[:, 0]
y = df.values[:, 1]
t = np.linspace(0, np.pi*2, n)

yn = np.outer(y, np.cos(t))
zn = np.outer(y, np.sin(t))


plt.plot(x, yn[:, :])
plt.show()



#MTTa = np.matrix([[sum(yn[:, :]**2 + yn[:, :]**2), sum(-df.values[:, 0]*yn[:, :]), sum(-df.values[:, 0]*yn[:, :])],
#[sum(-yn[:, :]*df.values[:, 0]), sum(df.values[:, 0]**2 + yn[:, :]**2), sum(-yn[:, :]*yn[:, :])],
#[sum(-yn[:, :]*df.values[:, 0]), sum(-yn[:, :]*yn[:, :]), sum(df.values[:, 0]**2+yn[:, :]**2)]])

xy = [0, 0, 0]
for i in range(n):
    for j in range(len(df)):
        xy = float(x[j]), float(yn[j, i]), float(yn[j, i])

print(xy)



