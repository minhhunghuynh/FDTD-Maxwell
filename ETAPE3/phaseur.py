import numpy as np
import matplotlib.pyplot as plt

c        = 3e8
mu0      = 4 * np.pi * 1e-7
epsilon0 = 1 / (mu0 * c**2)

f  = 1e9
L  = c / f
dx = L / 20
dt = dx / (2 * c)

Mx = 600
My = 600

m_source = Mx // 2
n_source = My // 2

Ez = np.zeros((Mx, My))
By = np.zeros((Mx, My))
Bx = np.zeros((Mx, My))

# Nombre de pas par période et durée totale
pas_T = round(1 / (f * dt))   # pas par période  
Q     = 600 # nombre total de pas 

# Stockage de tous les pas
E_enregistré = np.zeros((Q, Mx, My))

for q in range(Q):

    J = np.sin(2 * np.pi * f * q * dt)

    # 1. Mise à jour de B̃y  (Faraday selon x)
    By[0:Mx-1, :] = By[0:Mx-1, :] + 0.5 * (Ez[1:Mx, :] - Ez[0:Mx-1, :])

    # 2. Mise à jour de B̃x  (Faraday selon y)
    Bx[:, 0:My-1] = Bx[:, 0:My-1] - 0.5 * (Ez[:, 1:My] - Ez[:, 0:My-1])

    # 3. Mise à jour de Ez  (Ampère) + injection de la source
    Ez[1:Mx-1, 1:My-1] = (Ez[1:Mx-1, 1:My-1]
        + 0.5 * (By[1:Mx-1, 1:My-1] - By[0:Mx-2, 1:My-1])
        - 0.5 * (Bx[1:Mx-1, 1:My-1] - Bx[1:Mx-1, 0:My-2]))
    Ez[m_source, n_source] -= (dt / epsilon0) * J

    E_enregistré[q] = Ez

# Extraction de la dernière période en régime permanent
début_régime = Q - pas_T -10
Ez_permanent = E_enregistré[début_régime:, :, :]

# Amplitude du phaseur : (max - min) / 2 sur la période stockée
amplitude = (np.max(Ez_permanent, axis=0) - np.min(Ez_permanent, axis=0)) / 2

# Graphe
fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(amplitude.T, origin='lower',
               cmap='hot', vmin=0, vmax=amplitude.max(),
               extent=[0, Mx*dx*100, 0, My*dx*100])
plt.colorbar(im, ax=ax, label='|Ez| (V/m)')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('Amplitude du phaseur Ez 2D')
plt.tight_layout()
plt.show()

# Deuxième figure : décroissance en fonction de r 
r     = np.arange(1, Mx - m_source) * dx
amp_r = amplitude[m_source + 1:, n_source]
fig2, ax2 = plt.subplots(figsize=(7, 4))
ax2.plot(r, amp_r, color='steelblue')
ax2.set_xlabel('r (m)')
ax2.set_ylabel('|Ez| (V/m)')
ax2.set_title("Décroissance de l'amplitude avec r")
ax2.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
