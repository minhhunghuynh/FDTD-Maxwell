import numpy as np
import matplotlib.pyplot as plt

#paramètres
c    = 3e8
mu0  = 4 * np.pi * 1e-7
epsilon0 = 1 / (mu0 * c**2)
f    = 1e9
L = c / f
dx   = L / 20
dt   = dx / (2 * c)
M    = 500
Q    = 3000
m_src = M // 2
# initialisation des champs
E  = np.zeros(M)
B  = np.zeros(M)
E1 = np.zeros(M)
E2 = np.zeros(M)
E_stock = np.zeros((Q, M))
# mise à jour des champs 
for q in range(Q):
    for m in range(0, M - 1):
        B[m] = B[m] + 0.5 * (E[m + 1] - E[m])
    for m in range(1, M - 1):
        Jz = 0
        if m == m_src:
            Jz = np.sin(2 * np.pi * f * q * dt)
        E[m] = E[m] + 0.5 * (B[m] - B[m - 1]) - (dt / epsilon0) * Jz 
    E_stock[q, :] = E
    # conditions limites    
    E[0]   = E2[1]
    E[M-1] = E2[M-2]
    E2[:] = E1
    E1[:] = E
# Amplitude du phaseur (max - min) / 2 sur une période en régime permanent
pas_T        = round(1 / (f * dt))  # pas par période, T/dt donc 1/f*dt
début_régime = Q - pas_T - 10 # moment fixé en régime permanent (un peu avant la fin pour pas avoir de problème avec le dernier point)
E_permanent = E_stock[début_régime: début_régime + pas_T, :]
amplitude_phaseur = (np.max(E_permanent, axis=0) - np.min(E_permanent, axis=0)) / 2 # Maximum par colonne (axis=0) - Minimum par colonne, divisé par 2 pour obtenir l'amplitude du signal sinusoïdal
# Graphe
x = np.arange(M) * dx
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(x, amplitude_phaseur, color='black')
ax.set_ylim(0, amplitude_phaseur.max() * 1.2)
ax.set_xlabel("Position (m)") 
ax.set_ylabel("|E| (V/m)")
ax.set_title("Amplitude du phaseur")
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
