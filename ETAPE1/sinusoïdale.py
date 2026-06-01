import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
# Paramètres
c   = 3*10^8
mu0 = 4 * np.pi * 10-7
epsilon0 = 1 / (mu0 * c**2)
f   = 1e9 # fréquence de la source
L = c / f # longueur d'onde 
dx  = L / 20 # pas spatial (règle de bonne pratique dans la description du projet)
dt  = dx / (2 * c) # pas temporel 
M = 500 # nombre de pas spatiaux
Q = 1000 # nombre de pas temporels
msource = M / 2 # position de la source (on la place au centre)
# Champs
E  = np.zeros(M) # champ électrique Ez
B  = np.zeros(M) # champ magnétique (cBy)
E1 = np.zeros(M) # E au pas q-1  pour les CL
E2 = np.zeros(M) # E au pas q-2  pour les CL
# enregistrement pour l'animation
E_enregistrement = np.zeros((Q, M))
#  Boucle temporelle
for q in range(Q):
    #Mise à jour de B
    for m in range(0, M - 1):
        B[m] = B[m] + 0.5 * (E[m + 1] - E[m])
    #  Mise à jour de E  
    for m in range(0, M - 1):
        jz = 0
        if m == msource:
            jz = np.sin(2 * np.pi * f * q * dt)
        E[m] = E[m] + 0.5 * (B[m] - B[m - 1]) - (dt / epsilon0) * jz
# 4. Conditions limites 
    E[0]   = E2[1]
    E[M-1] = E2[M-2]
# mémorisation pour utiliser les conditions limites 
    E2[:] = E1
    E1[:] = E
    E_enregistrement[q, :] = E
    
# Animation (aidé par claude ai)
x = np.arange(M) * dx # axe x en mètres
Emax = np.max(np.abs(E_enregistrement))

fig, ax = plt.subplots(figsize=(10, 4))
line, = ax.plot(x, E_enregistrement[0], color='steelblue')
ax.set_xlim(x[0], x[-1])
ax.set_ylim(-1.5 * Emax, 1.5 * Emax)
ax.set_xlabel("Position (m)")
ax.set_ylabel("Ez (V/m)")
ax.set_title("1D  espace libre avec source sinusoïdale")
ax.grid(True, alpha=0.3)

def update(q):
    line.set_ydata(E_enregistrement[q])
    return line,

animation = FuncAnimation(fig, update, frames=range(0, Q, 5), interval=20, blit=True)
plt.tight_layout()
plt.show()