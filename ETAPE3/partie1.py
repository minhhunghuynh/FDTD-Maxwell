import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Paramètres
c        = 3e8
mu0      = 4 * np.pi * 1e-7
epsilon0 = 1 / (mu0 * c**2)

f  = 1e9           # fréquence source
L  = c / f         # longueur d'onde
dx = L / 20        # pas spatial
dt = dx / (2 * c)  # pas temporel

Mx = 600 # nombre de points selon x
My = 600 # nombre de points selon y
Q  = 600 # nombre de pas temporels

# Position de la source
m_source = Mx // 2
n_source = My // 2

# Champs 
Ez = np.zeros((Mx, My))   
By = np.zeros((Mx, My))  
Bx = np.zeros((Mx, My))  

# Stockage pour l'animation
E_enregistré = np.zeros((Q, Mx, My))

# mise à jour des champs
for q in range(Q):

    # Source sinusoïdale
    J = np.sin(2 * np.pi * f * q * dt)

    By[0:Mx-1, :] = By[0:Mx-1, :] + 0.5 * (Ez[1:Mx, :] - Ez[0:Mx-1, :])

    Bx[:, 0:My-1] = Bx[:, 0:My-1] - 0.5 * (Ez[:, 1:My] - Ez[:, 0:My-1])

    # 3. Mise à jour de Ez  (Ampère) + injection de la source en (m_source, n_source)
    Ez[1:Mx-1, 1:My-1] = (Ez[1:Mx-1, 1:My-1]
        + 0.5 * (By[1:Mx-1, 1:My-1] - By[0:Mx-2, 1:My-1])
        - 0.5 * (Bx[1:Mx-1, 1:My-1] - Bx[1:Mx-1, 0:My-2]))
    Ez[m_source, n_source] = Ez[m_source, n_source] - (dt / epsilon0) * J


    # Stockage
    E_enregistré[q] = Ez

# Animation
Emax = np.percentile(np.abs(E_enregistré), 99) if np.max(np.abs(E_enregistré)) > 0 else 1

fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(E_enregistré[0].T, origin='lower',
               cmap='RdBu', vmin=-Emax, vmax=Emax,
               extent=[0, Mx*dx*100, 0, My*dx*100])

plt.colorbar(im, ax=ax, label='Ez (V/m)')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('FDTD 2D – espace libre, source sinusoïdale')

def update(q):
    im.set_data(E_enregistré[q].T)
    return im,

anim = FuncAnimation(fig, update, frames=range(0, Q, 5), interval=50, blit=True)
plt.tight_layout()
plt.show()
