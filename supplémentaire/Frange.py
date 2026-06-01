import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Paramètres physiques
c        = 3e8
mu0      = 4 * np.pi * 1e-7
epsilon0 = 1 / (mu0 * c**2)
f = 2e9
L = c / f  # longueur d'onde

dx = L / 20
dy = dx
dt = dx / (2 * c)

Nx = 400
Ny = 400

# Positions des deux sources
d_lambda  = 3
d_cellule = int(d_lambda * L / dx)

centre_x = Nx // 2
centre_y = Ny // 2

m_source1 = centre_x - d_cellule // 2
m_source2 = centre_x + d_cellule // 2
n_source = centre_y

for n in range(-d_lambda, d_lambda+1):
    sin_t = n / d_lambda

# Champs
Ez = np.zeros((Nx, Ny))
Bx = np.zeros((Nx, Ny))
By = np.zeros((Nx, Ny))

# Boucle temporelle
Q = 400
pas_T  = round(1 / (f * dt))  # nombre de pas par période

E_enregistré = np.zeros((Q, Nx, Ny))

for q in range(Q):
    J = np.sin(2 * np.pi * f * q * dt)

    By[0:Nx-1, :] = By[0:Nx-1, :] + 0.5 * (Ez[1:Nx, :] - Ez[0:Nx-1, :])
    Bx[:, 0:Ny-1] = Bx[:, 0:Ny-1] - 0.5 * (Ez[:, 1:Ny] - Ez[:, 0:Ny-1])

    Ez[1:Nx-1, 1:Ny-1] = (Ez[1:Nx-1, 1:Ny-1]
        + 0.5 * (By[1:Nx-1, 1:Ny-1] - By[0:Nx-2, 1:Ny-1])
        - 0.5 * (Bx[1:Nx-1, 1:Ny-1] - Bx[1:Nx-1, 0:Ny-2]))
    Ez[m_source1, n_source] -= (dt / epsilon0) * J
    Ez[m_source2, n_source] -= (dt / epsilon0) * J



    E_enregistré[q] = Ez


# Amplitude sur la dernière période (méthode phaseur)
Ez_permanent = E_enregistré[Q - pas_T:, :, :]
amplitude = (np.max(Ez_permanent, axis=0) - np.min(Ez_permanent, axis=0)) / 2

# Figure 1 
fig1, ax1 = plt.subplots(figsize=(7, 6))
vmax = np.percentile(np.abs(Ez), 99) * 0.6
im1 = ax1.imshow(
    Ez.T, origin='lower',
    extent=[0, Nx*dx*100, 0, Ny*dy*100],
    cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='equal'
)
plt.colorbar(im1, ax=ax1, label='Ez [u.a.]', shrink=0.8)
ax1.set_xlabel('x [cm]')
ax1.set_ylabel('y [cm]')
ax1.set_title(f'FDTD 2D — Interférences\ndeux sources en phase, d = {d_lambda}λ')
plt.tight_layout()
plt.show()

# Figure 2 
fig2, ax2 = plt.subplots(figsize=(7, 6))
im2 = ax2.imshow(
    amplitude.T, origin='lower',
    extent=[0, Nx*dx*100, 0, Ny*dy*100],
    cmap='hot', aspect='equal'
)
plt.colorbar(im2, ax=ax2, label='|Ez| [u.a.]', shrink=0.8)

centre_x = centre_x * dx * 100
centre_y = centre_y * dy * 100
r_max = min(Nx, Ny) * dx * 100 * 0.45
for n in range(-d_lambda, d_lambda+1):
    sin_t = n / d_lambda
    if abs(sin_t) < 1:
        theta = np.arcsin(sin_t)
        xe = centre_x + r_max * np.sin(theta)
        ye = centre_y + r_max * np.cos(theta)
        ax2.annotate('', xy=(xe, ye), xytext=(centre_x, centre_y),
                     arrowprops=dict(arrowstyle='->', color='cyan', lw=1.5, alpha=0.8))
        ax2.text(xe, ye, f'n={n:+d}', color='cyan', fontsize=8, ha='center', va='bottom')

ax2.set_xlabel('x [cm]')
ax2.set_ylabel('y [cm]')
ax2.set_title(f'Amplitude |Ez| — Franges d\'interférence\nd = {d_lambda}λ  →  maxima à d·sinθ = nλ')
plt.tight_layout()
plt.show()

# Animation
Emax = np.percentile(np.abs(E_enregistré), 99) if np.max(np.abs(E_enregistré)) > 0 else 1

fig3, ax3 = plt.subplots(figsize=(7, 6))
im3 = ax3.imshow(E_enregistré[0].T, origin='lower',
                 cmap='RdBu_r', vmin=-Emax, vmax=Emax,
                 extent=[0, Nx*dx*100, 0, Ny*dy*100], aspect='equal')
plt.colorbar(im3, ax=ax3, label='Ez [u.a.]', shrink=0.8)
ax3.set_xlabel('x [cm]')
ax3.set_ylabel('y [cm]')
ax3.set_title(f'FDTD 2D — Interférences, d = {d_lambda}λ')

def update(q):
    im3.set_data(E_enregistré[q].T)
    return im3,

anim = FuncAnimation(fig3, update, frames=range(0, Q, 5), interval=50, blit=True)
plt.tight_layout()
plt.show()
