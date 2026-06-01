import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Paramètres physiques
c = 3e8
f = 2e9
L = c / f

epsilon_r1 = 2
epsilon_r2 = 4
theta_i    = 30

# Angle transmis
sin_t       = np.sqrt(epsilon_r1 / epsilon_r2) * np.sin(np.radians(theta_i))
theta_t_deg = np.degrees(np.arcsin(sin_t))

# Discrétisation
dx = L / 20
dy = dx
dt = dx / (2 * c)

Nx  = 600
Ny   = 400
n_interface = 300
n_source = 200

# Source étendue pour créer un front d'onde plan approximatif
m_source_debut = Nx // 3
m_source_fin   = 2 * Nx // 3
m_arr = np.arange(m_source_debut, m_source_fin)

beta1 = 2 * np.pi * f * np.sqrt(epsilon_r1) / c
dphi  = beta1 * np.sin(np.radians(theta_i)) * dx

# Permittivité
epsilon_r = np.ones((Nx, Ny))
epsilon_r[:, n_interface:] = epsilon_r2

# Champs
Ez = np.zeros((Nx, Ny))
Bx = np.zeros((Nx, Ny))
By = np.zeros((Nx, Ny))

# Boucle temporelle
Q  = 1000
N_periodes  = 1 # on emet qu'une fois une période pour observer la  loi
q_source_stop = int(N_periodes / (f * dt))

E_enregistré = np.zeros((Q, Nx, Ny))

for q in range(Q):
    Bx[:, :-1] -= 0.5 * (Ez[:, 1:] - Ez[:, :-1])
    By[:-1, :] += 0.5 * (Ez[1:, :] - Ez[:-1, :])
    Ez[1:-1, 1:-1] += (1.0 / epsilon_r[1:-1, 1:-1]) * (
          0.5 * (By[1:-1, 1:-1] - By[:-2,  1:-1])
        - 0.5 * (Bx[1:-1, 1:-1] - Bx[1:-1, :-2])
    )
    if q < q_source_stop:
        phase = 2 * np.pi * f * (q + 0.5) * dt - m_arr * dphi
        Ez[m_source_debut:m_source_fin, n_source] += 0.5 * np.sin(phase)

    E_enregistré[q] = Ez

# Fenêtre de zoom autour de l'interface (Aidé par l'ia)
n_zoom_debut = n_interface - 50
n_zoom_fin   = n_interface + 75

y0 = n_zoom_debut * dy
y1 = n_zoom_fin   * dy
x0 = 0
x1 = Nx * dx

def zoom(frame):
    return frame[:, n_zoom_debut:n_zoom_fin].T

# Animation
vmax = np.percentile(np.abs(E_enregistré), 99.5) * 0.8

fig, ax = plt.subplots(figsize=(8, 5))
ax.axhspan(n_interface * dy, y1, color='steelblue', alpha=0.35, zorder=0)

im = ax.imshow(
    zoom(E_enregistré[0]),
    origin='lower',
    extent=[x0, x1, y0, y1],
    cmap='RdBu_r',
    vmin=-vmax, vmax=vmax,
    aspect='equal',
    zorder=1,
    alpha=0.85
)
plt.colorbar(im, ax=ax, label='Ez', shrink=0.8)

ax.axhline(n_interface * dy, color='white', lw=3, ls='--', zorder=2,
           label=f'Interface  n₁={np.sqrt(epsilon_r1):.1f} → n₂={np.sqrt(epsilon_r2):.1f}')

f_itf = (n_interface - n_zoom_debut) / (n_zoom_fin - n_zoom_debut)
ax.text(0.02, f_itf * 0.45,
        f'Milieu 1\nεr={epsilon_r1:.0f}  n={np.sqrt(epsilon_r1):.1f}',
        transform=ax.transAxes, fontsize=9, va='center',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
ax.text(0.02, f_itf + (1 - f_itf) * 0.55,
        f'Milieu 2\nεr={epsilon_r2:.0f}  n={np.sqrt(epsilon_r2):.1f}',
        transform=ax.transAxes, fontsize=9, va='center',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.legend(loc='upper right', fontsize=9)
ax.set_title('FDTD 2D — Loi de Snell')

def update(frame):
    im.set_data(zoom(E_enregistré[frame]))
    return [im]

anim = FuncAnimation(fig, update, frames=range(0, Q, 5), interval=60, blit=True)
plt.tight_layout()
plt.show()


# pour incliner l'onde incidente, il fallait introduire un déphasage linéaire le lon de source. L'ia a été utilisé pour effectuer cela
