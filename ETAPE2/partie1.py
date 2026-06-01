import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import cmath

# Paramètres
c   = 3e8
mu0 = 4 * np.pi * 1e-7
eps0 = 1 / (mu0 * c**2)

f   = 1e9 # fréquence source 
L = c / f # longueur d'onde [m]
dx  = L / 20 # pas spatial
dt  = dx / (2 * c)# pas temporel 

M = 800 # nombre de points spatiaux
Q = 1500 # nombre de pas temporels
m_src = M // 2 # position de la source (milieu du domaine)

# Permittivité relative εr(m)
eps_r = np.ones(M)      # vaut 1 partout par défaut (espace libre)
eps_r[600:700] = 4.0    # bloc de diélectrique (ex : verre, εr = 4)

# Champs 
E  = np.zeros(M)# champ électrique Ez
B  = np.zeros(M)# champ magnétique
E1 = np.zeros(M)# E au pas q-1 
E2 = np.zeros(M)# E au pas q-2  

# Stockage pour l'animation 
E_stock = np.zeros((Q, M))

# Boucle temporelle 
for q in range(Q):

    # Mise à jour de B
    for m in range(0, M - 1):
        B[m] = B[m] + 0.5 * (E[m + 1] - E[m])

    # Mise à jour de E  
    for m in range(1, M - 1):
        Jz = 0
        if m == m_src:
            Jz = np.sin(2 * np.pi * f * q * dt)
        E[m] = E[m] + (1 / (2 * eps_r[m])) * (B[m] - B[m - 1]) - (dt / (eps0 * eps_r[m])) * Jz 

    # Conditions limites absorbantes 
    E[0]   = E2[1]
    E[M-1] = E2[M-2]

    # Décalage mémoire
    E2[:] = E1
    E1[:] = E

    E_stock[q, :] = E

# Animation
x = np.arange(M) * dx * 100 # axe x en centimètres
Emax = np.max(np.abs(E_stock))

fig, ax = plt.subplots(figsize=(14, 6))
line, = ax.plot(x, E_stock[0], color='steelblue')
ax.set_xlim(x[0], x[-1])
ax.set_ylim(-1.5 * Emax, 1.5 * Emax)
ax.set_xlabel("Position (cm)")
ax.set_ylabel("Ez (V/m)")
ax.set_title("FDTD 1D – étape 2 : diélectrique sans pertes")
ax.grid(True, alpha=0.3)

def update(q):
    line.set_ydata(E_stock[q])
    return line,

anim = FuncAnimation(fig, update, frames=range(0, Q, 2), interval=40, blit=True)
plt.tight_layout()
plt.show()