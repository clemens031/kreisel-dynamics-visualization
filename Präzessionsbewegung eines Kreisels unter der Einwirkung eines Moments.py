from sympy.abc import *
import sympy as sp
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

#Massenträgheitstensor

df = pd.read_csv('data/koordinatenpunkte Kreisel.py') # Profil
#df = pd.read_csv('data/Einfacher Kreisel.py')
df.loc[len(df)+1]=df.loc[0]
teilprofil = np.zeros([len(df)-1, 8])

for i in range(len(df)-1):
    teilprofil[i, 0] = df.values[i, 0] + (df.values[i+1, 0]-df.values[i, 0])/2 # X-Koordinate des Teilschwerpunkts
    teilprofil[i, 1] = df.values[i, 1] + (df.values[i+1, 1]-df.values[i ,1])/2 # Y-Koordinate des Teilschwerpunkts
    teilprofil[i, 2] = df.values[i, 0]
    teilprofil[i, 3] = abs(df.values[i, 0] - df.values[i+1, 0]) #höhe
    teilprofil[i, 4] = df.values[i, 1] #radius

rho =  0.008 #Dichte in g pro mm³ für Edelstahl 1.4539

#schwerpunkt = 1/sum(teilprofil[:, 3] * teilprofil[:, 4] * 2) * sum(teilprofil[:, 0] * teilprofil[:, 3] * teilprofil[:, 4] * 2) #Schwerpunkt
schwerpunkt = 30.e-3
#steiner = abs(schwerpunkt - teilprofil[:, 0]) * (teilprofil[:, 3] * teilprofil[:, 4] * 2) ** 2

#trägheitsmoment = sum(0.5 * np.pi * rho * teilprofil[:, 3] * teilprofil[:, 4] ** 4 + steiner)

#radiusaverage = sum(df.values[:, 1]) / len(df)
radiusaverage = 75 #...
#print(f'Trägheitsmoment um Rotationsachse wenn Kreisel als Zylinder vereinfacht wird: {0.5 * np.pi * rho * 435 * radiusaverage **4}')


#Massenträgheitstensor

MTT = np.matrix([[sum(df.values[:, 1]**2 + df.values[:, 2]**2), sum(-df.values[:, 0]*df.values[:, 1]), sum(-df.values[:, 0]*df.values[:, 2])],
[sum(-df.values[:, 1]*df.values[:, 0]), sum(df.values[:, 0]**2 + df.values[:, 2]**2), sum(-df.values[:, 1]*df.values[:, 2])],
[sum(-df.values[:, 2]*df.values[:, 0]), sum(-df.values[:, 1]*df.values[:, 2]), sum(df.values[:, 0]**2+df.values[:, 1]**2)]])

#Hauptinvarianten
HI1 = MTT[0, 0] + MTT[1, 1] + MTT[2, 2]
HI2 = MTT[0, 0]*MTT[1, 1] + MTT[2, 2]*MTT[1, 1] + MTT[0, 0]*MTT[2, 2] - MTT[0, 1]**2-MTT[1, 2]**2-MTT[0, 2]**2
HI3 = np.linalg.det(MTT)
print(f'Hauptträgheitsmomente: {HI1, HI2, HI3} \nSchwerpunkt: {schwerpunkt} \n ---')

#Plotten vom Querschnitt des Kreisels und einem 3d Modell
n = 100

fig = plt.figure(figsize=(12,6))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122,projection='3d')
y = df.values[:, 0]
g = df.values[:, 1]
t = np.linspace(0, np.pi*2, n)

xn = np.outer(g, np.cos(t))
yn = np.outer(g, np.sin(t))
zn = np.zeros_like(xn)

for i in range(len(g)):
    zn[i:i+1,:] = np.full_like(zn[0,:], y[i])

ax1.plot(g, y)
ax2.plot_surface(xn, yn, zn)
plt.show()


#omega = Winkelgeschwindigkeitsvektor
#phi = Präzession
#theta = Nutation
#psi = Eigenrotation

theta = sp.pi/6 #Nutationswinkel, wenn der eingegebene Winkel keine reelle Lösung liefert wird ein fehler Ausgegeben


volumen = np.pi * radiusaverage ** 2 * df.values[-2, 0]
m = rho * volumen#Masse in kg
l = schwerpunkt #Hebelarm (Schwerpunkt) in m
g = 9.81 #Gravitation auf der Erdoberfläche m/s

MTM1 = HI1 #Massenträgheitsmoment + Korrekturfaktor für sinnvolle Ergebnisse da die Daten auf der größe des ursprünglichen Bildes beruhen
MTM3 = HI3 #Hauptinvariante ist Hauptmassenträgheitsmoment
psi_punkt = x #Eigenrotationsgeschwindigkeit
phi_punkt = m*g*l /(MTM3*psi_punkt) #Präzessionsgeschwindigkeit, gilt nur wenn die Eigenrotationsgeschwindigkeit sehr groß ist


omega = [[-phi_punkt * sp.sin(theta)],[0],[(psi_punkt + phi_punkt * sp.cos(theta))]] #Winkelgeschwindigkeit im Ausgangskoordinatensystem #wird nicht benötigt

omegaf = sp.Matrix([-phi_punkt * sp.sin(theta),0.,phi_punkt * sp.cos(theta)]) #Winkelgeschwindigkeit des mitrotierenden Koordinatensystems

drall = sp.Matrix([-MTM1 * phi_punkt * sp.sin(theta),0., MTM3 * (psi_punkt + phi_punkt * sp.cos(theta))]) #Drall im gedrehten Koordinatensystem

moment = omegaf.cross(drall) #Ableitung vom Drall, Moment durch Änderung vom Drall
momenthändisch = phi_punkt * sp.sin(theta) *(MTM3 * (psi_punkt + phi_punkt*sp.cos(theta)-phi_punkt*MTM1*sp.cos(theta))) #Moment durch Ableitung vom Drall aus Skript um Abweichung festzustellen ->Abweichung ist minimal
omega3stern = psi_punkt + phi_punkt*sp.cos(theta) #Wingengeschwindigkeitsvekotor omega 3 im gedrehten Koordinatensystem bei konstanter Eigenrotationsgeschwindigkeit

A = MTM3 * omega3stern / (2 * MTM1 * sp.cos(theta))
B = sp.sqrt(1 - 4 * moment[1] * MTM1 * sp.cos(theta) / (MTM3 ** 2 * omega3stern **2 * sp.sin(theta)))
phi_punkt1 = A*(1-B) #langsame Präzession
phi_punkt2 = A*(1+B) #schnelle Präzession


# m*g*l muss so groß sein wie M0 bzw. L0_punkt, also omegaf * drall

momenteig = m*g*l*sp.sin(theta) #Moment durch Eigengewicht

expr=sp.Eq(momenteig, moment[1]) #Aufstellen der Gleichung
psi_punktl = sp.solve(expr,x) #Lösen der Gleichung
print(psi_punktl)
# Wie schnell müssen sich der Kreisel drehen damit Moment -> 0 und phi_punkt bzw. Präzessionsgeschwindigkeit -> 0

#print(f'Nicht Imaginär: {omega3stern.subs(x, psi_punktl[1]) > 4*MTM1 / ((MTM3**2) * m*g*l*sp.cos(theta))}')
print(f'Nutationswikel: {theta / (2*sp.pi) * 360} Grad \nLangsame Präzession: {phi_punkt1.subs(x, psi_punktl[1]).evalf(5)} muss gegen 0 laufen \nSchnelle Präzession: {phi_punkt2.subs(x, psi_punktl[1]).evalf(5)} Präzessionsgeschwindigkeit wenn Moment durch Eigengewicht gegen 0 läuft')
print(f'Eigenrotationsgeschwindigkeit: {psi_punktl[1]} <- So schnell muss sich der Kreisel drehen\nPräzessionsgeschwindigkeit: {phi_punkt.subs(x, psi_punktl[1])} ist gleich der langsamen Präzession')

# Mindestdrehzahl
# omega_req = 4*theta_I/theta_III^2*mass*gravity*lage_schwerpunkt*np.cos(nutation)
omega_req = sp.sqrt(4*MTM1/MTM3**2*m*9.81*schwerpunkt*sp.cos(theta))
n_req = omega_req*60. / (2*np.pi)
psi_langsam = m*9.81*schwerpunkt / (MTM3*omega_req)
n_psi_langsam = psi_langsam*60. / (2*np.pi)
psi_schnell = MTM3*omega_req / ((MTM1-MTM3)*sp.cos(theta))
n_psi_schnell = psi_schnell*60. / (2*sp.pi)

print(f'psi_langsam: {psi_langsam}, {n_psi_langsam} psi_schnell: {psi_schnell.evalf(5)}, {n_psi_schnell.evalf(5)}')

print(MTT)