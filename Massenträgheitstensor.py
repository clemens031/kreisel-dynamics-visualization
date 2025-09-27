import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Ansatz: Querschnitt wird aus Rechteckflächen zusammengesetzt
# teilprofil - Teilprofilflächen
# ftm - Flächenträgheitsmoment
# htm - Hauptflächenträgheitsmoment
# Profilteilflächen werden definiert über Schwerachsenverlauf, Profildicke und -länge

#Auswahl geschlossenes (g = 1) oder offenes (g = 0) Profil:
g = 1

rho = 0.0001 # Dichte in g/mm³
length = 1. # Länge des Profils

# Daten der Teilprofilflächen
df = pd.read_csv('data/Daten Viereck.py') # geschlossenes Profil
#df = pd.read_csv('data/koordinatenpunkte Kreisel.py') # offenes Profil

if g == 1:
    df.loc[len(df)+1]=df.loc[0]
teilprofil = np.zeros([len(df)-1, 8])

for i in range(len(df)-1):
    teilprofil[i, 0] = df.values[i, 0] + (df.values[i+1, 0]-df.values[i, 0])/2 # X-Koordinate des Teilschwerpunkts
    teilprofil[i, 1] = df.values[i, 1] + (df.values[i+1, 1]-df.values[i ,1])/2 # Y-Koordinate des Teilschwerpunkts
    teilprofil[i, 2] = np.sqrt((df.values[i+1, 0]-df.values[i, 0])**2 + (df.values[i+1, 1]-df.values[i, 1])**2) # Länge
    teilprofil[i, 3] = np.arctan2 ((df.values[i+1, 1]-df.values[i, 1]),(df.values[i+1, 0]-df.values[i, 0])) # Winkel
    J22 = df.values[i, 2]**3*teilprofil[i, 2]/12
    J33 = df.values[i, 2]*teilprofil[i, 2]**3/12
    teilprofil[i, 4] = (J22 + J33)/2+(J22-J33)*np.cos(2*teilprofil[i, 3])/2 # J22
    teilprofil[i, 5] = (J22 + J33)/2-(J22-J33)*np.cos(2*teilprofil[i, 3])/2 # J33
    teilprofil[i, 6] = -(J22-J33)*np.sin(-2*teilprofil[i, 3])/2 # J23
    teilprofil[i, 7] = teilprofil[i, 2]*df.values[i, 2] # Fläche des Teilprofils

schwerpunkt = [sum(teilprofil[:, 0]*teilprofil[:, 2])/sum(teilprofil[:, 2]), sum(teilprofil[:, 1]*teilprofil[:, 2])/sum(teilprofil[:, 2])]
masse = rho*length*sum(teilprofil[:, 7])

ftm = np.zeros(4)
ftm[0] = sum(teilprofil[:, 4] + (teilprofil[:, 1]-schwerpunkt[1])**2 * teilprofil[:, 7]) # J22
ftm[1] = sum(teilprofil[:, 5] + (teilprofil[:, 0]-schwerpunkt[0])**2 *teilprofil[:, 7]) # J33
ftm[2] = sum(teilprofil[:, 6] - (teilprofil[:, 0]-schwerpunkt[0])*(teilprofil[:, 1] - schwerpunkt[1])*teilprofil[:, 7]) # J23
ftm[3] = ftm[0] + ftm[1] # J11

htm = [0.0,0.0]
htm[0] = (ftm[0]+ftm[1])/2+np.sqrt((ftm[0]-ftm[1])**2/4+ftm[2]**2)
htm[1] = (ftm[0]+ftm[1])/2-np.sqrt((ftm[0]-ftm[1])**2/4+ftm[2]**2)

mtm_s = np.zeros(4)
mtm_s[:] = rho*length*ftm[:] # theta22, theta33, theta23, theta11

#Diagramm zeichnen
fig, ax = plt.subplots()
ax.set(xlabel='$x_2$', ylabel='$x_3$')
#ax.set_xlim(1.5*np.min(df['x2']), 1.5*np.max(df['x2']))
ax.plot(df.values[:,0],df.values[:,1])
ax.plot(schwerpunkt[0],schwerpunkt[1], 'bo')
ax.grid()
plt.show()

print('Schwerpunkt [X,Y]:')
print(schwerpunkt)
print()
print('Flächenträgheitmoment im Schwerpunkt bezogen auf Ausgangskoordinatensystem [J22, J33, J23]: ')
print(ftm)
print()
print('Hauptflächenträgheitsmoment [J_I,J_II]:')
print(htm)
print('Massenträgheitsmomente bezogen auf den Schwerpunkt [theta22, theta33, theta23, theta11]')
print(mtm_s)
print(schwerpunkt)