import numpy as np
import matplotlib.pyplot as plt

#  Constantes physiques
c    = 3e8
mu0  = 4 * np.pi * 1e-7
eps0 = 1 / (mu0 * c**2)

# Paramètres de discrétisation
f   = 1e9
L   = c / f
dx  = L / 20
dy  = dx
dt  = dx / (2 * c)

Mx = 400
My = 400
Q  = 400

m_source = Mx // 2
n_source = My // 2


Ra  = 73.0   # résistance de rayonnement

# Hauteur équivalente
theta    = np.pi / 2
h_e_perp = (L / np.pi) * np.cos(np.pi / 2 * np.cos(theta)) / np.sin(theta)


Ez = np.zeros((Mx, My))
Bx = np.zeros((Mx, My))
By = np.zeros((Mx, My))

# Stockage d'une période en régime permanent
pas_T       = round(1 / (f * dt))  # nombre de pas par période
Q_permanent = Q - pas_T # début du stockage
Ez_period   = np.zeros((pas_T, Mx, My))
# Mise à jour des champs
for q in range(Q):

    Bx[:, :-1] -= 0.5 * (Ez[:, 1:] - Ez[:, :-1])
    By[:-1, :] += 0.5 * (Ez[1:, :] - Ez[:-1, :])
    Ez[1:-1, 1:-1] += (
        + 0.5 * (By[1:-1, 1:-1] - By[:-2,  1:-1])
        - 0.5 * (Bx[1:-1, 1:-1] - Bx[1:-1,  :-2])
    )
    Ez[m_source, n_source] -= (dt / eps0) * np.sin(2 * np.pi * f * q * dt)

    if q >= Q_permanent:
        Ez_period[q - Q_permanent] = Ez


# Amplitude du phaseur : (max - min) / 2 sur la dernière période
Ez_amplitude = (np.max(Ez_period, axis=0) - np.min(Ez_period, axis=0)) / 2

# Tension en circuit ouvert  
V_oc = h_e_perp * Ez_amplitude  

# Puissance collectée maximale  
P_L = V_oc**2 / (8 * Ra)            


borne = [0, Mx * dx , 0, My * dy ]   # axes

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Ez amplitude
im0 = axes[0].imshow(Ez_amplitude.T, origin='lower', cmap='hot', extent=borne)
plt.colorbar(im0, ax=axes[0], label='Amplitude $|E_z|$ (V/m)')
axes[0].set_title('Amplitude du champ électrique $|E_z|$')
axes[0].set_xlabel('x (m)');  axes[0].set_ylabel('y (m)')

# Voc
im1 = axes[1].imshow(V_oc.T, origin='lower', cmap='hot', extent=borne)
plt.colorbar(im1, ax=axes[1], label='$|V_{oc}|$ (V)')
axes[1].set_title(r'Tension en circuit ouvert $|V_{oc}|$' )
axes[1].set_xlabel('x (m)');  axes[1].set_ylabel('y (m)')

# P_L
im2 = axes[2].imshow(P_L.T, origin='lower', cmap='hot', extent=borne,
                     vmax=np.percentile(P_L, 99.9))
plt.colorbar(im2, ax=axes[2], label='$P_L$ (W)')
axes[2].set_title(r'Puissance collectée $P_L$' )
axes[2].set_xlabel('x (m)');  axes[2].set_ylabel('y (m)')

plt.suptitle(f'FDTD 2D – Dipôle $\lambda/2$ récepteur à {f/1e9:.0f} GHz', fontsize=13)
plt.tight_layout()
plt.show()


r_max  = Mx - m_source - 1
k      = np.arange(r_max)
r_m   = k * dx          # distance depuis la source 

thetas_deg = [90, 45, 0]
colors     = ['tab:blue', 'tab:orange', 'tab:red']

fig2, ax2 = plt.subplots(figsize=(8, 5))

for theta_deg, col in zip(thetas_deg, colors):
    theta_rad = np.deg2rad(theta_deg)
    if theta_deg == 0:
        h_e = 0.0
    else:
        h_e = (L / np.pi) * np.cos(np.pi / 2 * np.cos(theta_rad)) / np.sin(theta_rad)

    P_theta = (h_e * Ez_amplitude) ** 2 / (8 * Ra)
    ax2.plot(r_m, P_theta[m_source + k, n_source],
             label=rf'$\theta = {theta_deg}°$', color=col)

ax2.set_xlabel('Distance depuis la source (m)')
ax2.set_ylabel('Puissance collectée $P_L$ (W)')
ax2.set_title(r'$P_L$ pour 3 orientations du dipôle récepteur')
ax2.legend()
ax2.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()
