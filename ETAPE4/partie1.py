import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Paramètres
c        = 3e8
mu0      = 4 * np.pi * 1e-7
epsilon0 = 1 / (mu0 * c**2)

f     = 1e9
lam   = c / f
dx    = lam / 20
dt    = dx / (2 * c)


Mx = 400
My = 400
Q  = 400

# Source : centre du domaine
m_src = Mx // 2
n_src = My // 2

# Champs
Ez = np.zeros((Mx, My))
By = np.zeros((Mx, My))
Bx = np.zeros((Mx, My))

Ez_stock = np.zeros((Q, Mx, My))

# ── Obstacle PEC : un seul cylindre ──────────────────────
# Centre au milieu du domaine, rayon = λ
obstacle_x = int(0.75*Mx)
obstacle_y = My // 2
rayon  = 5         

x_min = obstacle_x - rayon
x_max = obstacle_x + rayon
y_min = obstacle_y - rayon
y_max = obstacle_y + rayon
# ─────────────────────────────────────────────────────────

# Boucle temporelle
for q in range(Q):

    J = np.sin(2 * np.pi * f * q * dt)

    By[0:Mx-1, :] = By[0:Mx-1, :] + 0.5 * (Ez[1:Mx, :] - Ez[0:Mx-1, :])
    Bx[:, 0:My-1] = Bx[:, 0:My-1] - 0.5 * (Ez[:, 1:My] - Ez[:, 0:My-1])

    Ez[1:Mx-1, 1:My-1] = (Ez[1:Mx-1, 1:My-1]
        + 0.5 * (By[1:Mx-1, 1:My-1] - By[0:Mx-2, 1:My-1])
        - 0.5 * (Bx[1:Mx-1, 1:My-1] - Bx[1:Mx-1, 0:My-2]))
    Ez[m_src, n_src] -= (dt / epsilon0) * J

 
    Ez[x_min:x_max, y_min:y_max] = 0.0

    Ez_stock[q] = Ez

# Animation
Emax = np.percentile(np.abs(Ez_stock), 99) if np.max(np.abs(Ez_stock)) > 0 else 1

fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(Ez_stock[0].T, origin='lower',
               cmap='RdBu', vmin=-Emax, vmax=Emax,
               extent=[0, Mx*dx*100, 0, My*dx*100])

plt.colorbar(im, ax=ax, label='Ez (V/m)')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('FDTD 2D avec obstacle conducteur parfait')

def update(q):
    im.set_data(Ez_stock[q].T)
    return im,

anim = FuncAnimation(fig, update, frames=range(0, Q, 5), interval=50, blit=True)
plt.tight_layout()
plt.show()