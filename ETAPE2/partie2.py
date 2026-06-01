import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Paramètres 
c    = 3e8
mu0  = 4 * np.pi * 1e-7
epsilon0 = 1 / (mu0 * c**2)

f   = 1e9   # fréquence source 
L = c / f   # longueur d'onde
dx  = L / 20  # pas spatial
dt  = dx / (2 * c) # pas temporel 

M     = 400   # nombre de points spatiaux
Q     = 1000  # nombre de pas temporels
m_source = M // 3  # position de la source

# Matériaux
eps_r = np.ones(M)  # permittivité relative
sigma = np.zeros(M) # conductivité 

eps_r[250:330] = 2  # bloc diélectrique 
sigma[250:330] = 0.01  # conductivité 


# Champs 
E  = np.zeros(M) 
B  = np.zeros(M) 
E1 = np.zeros(M) 
E2 = np.zeros(M) 

# stockage
E_enregistré = np.zeros((Q, M))


for q in range(Q):

    # Mise à jour des champs
    for m in range(0, M - 1):
        B[m] = B[m] + 0.5 * (E[m + 1] - E[m])
    for m in range(1, M - 1):
        Jz = 0
        if m == m_source:
            Jz = np.sin(2 * np.pi * f * q * dt)
        E[m] = (1 - sigma[m] * dt / (eps_r[m] * epsilon0)) * E[m] + (1 / (2 * eps_r[m])) * (B[m] - B[m - 1]) - (dt / (eps_r[m] * epsilon0)) * Jz

    # Conditions limites 
    E[0]   = E2[1]
    E[M-1] = E2[M-2]

    # mémorisation
    E2[:] = E1
    E1[:] = E

    E_enregistré[q, :] = E

# Animation
x    = np.arange(M) * dx  # axe x
Emax = np.max(np.abs(E_enregistré))

fig, ax = plt.subplots(figsize=(14, 6))

# Visualisation du bloc diélectrique
ax.axvspan(250 * dx * 100, 330 * dx * 100, alpha=0.15, color='orange',
           label=f'Diélectrique (eps_r={eps_r[280]:.0f}, sigma={sigma[280]} S/m)')

line, = ax.plot(x, E_enregistré[0], color='steelblue')
ax.set_xlim(x[0], x[-1])
ax.set_ylim(-1.5 * Emax, 1.5 * Emax)
ax.set_xlabel("Position (m)")
ax.set_ylabel("Ez (V/m)")
ax.set_title("FDTD 1D - etape 2 : dielectrique avec pertes")
ax.legend(loc='upper right')
ax.grid(True, alpha=0.3)

def update(q):
    line.set_ydata(E_enregistré[q])
    return line,

anim = FuncAnimation(fig, update, frames=range(0, Q, 4), interval=30, blit=True)
plt.tight_layout()
plt.show()