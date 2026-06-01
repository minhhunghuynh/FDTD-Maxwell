import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Paramètres 
c    = 3e8
mu0  = 4 * np.pi * 1e-7
epsilon0 = 1 / (mu0 * c**2)

f   = 1e9    # fréquence source 
lam = c / f   # longueur d'onde 
dx  = lam / 20   # pas spatial
dt  = dx / (2 * c) # pas temporel 

M     = 400  # nombre de points spatiaux
Q     = 1000  # nombre de pas temporels
m_src = M // 3  # position de la source 

epsilon_r = np.ones(M) 
sigma = np.zeros(M)   

epsilon_r[250:330] = 4.0   
sigma[250:330] = 0.01  

# ── Champs ───────────────────────────────────────────────────────────────────
E  = np.zeros(M)        # champ électrique Ez
B  = np.zeros(M)        # champ magnétique B̃y = c·By
E1 = np.zeros(M)        # E au pas q-1  (pour conditions limites)
E2 = np.zeros(M)        # E au pas q-2  (pour conditions limites)

# ── Stockage pour l'animation ────────────────────────────────────────────────
E_enregistré = np.zeros((Q, M))

# ── Boucle temporelle ────────────────────────────────────────────────────────
for q in range(Q):

    # 1. Mise à jour de B̃  (éq. 16 du projet, inchangée)
    for m in range(0, M - 1):
        B[m] = B[m] + 0.5 * (E[m + 1] - E[m])

    # 2. Mise à jour de E  (approximation par moyenne, dérivée de Maxwell)
    #    E^{q+1}[m] = (1-alpha[m])/(1+alpha[m]) * E^q[m]
    #               + 1/(2*epsilon_r[m]*(1+alpha[m])) * (B[m] - B[m-1])
    #               - dt/(epsilon0*epsilon_r[m]*(1+alpha[m])) * Jz
    for m in range(1, M - 1):
        Jz = 0
        if m == m_src:
            Jz = np.sin(2 * np.pi * f * q * dt)
        E[m] = ((1 - sigma[m]*dt / (2*epsilon0*epsilon_r[m])) / (1 + sigma[m]*dt / (2*epsilon0*epsilon_r[m]))) * E[m] \
             + (1 / (2 * epsilon_r[m] * (1 + sigma[m]*dt / (2*epsilon0*epsilon_r[m])))) * (B[m] - B[m - 1]) \
             - (dt / (epsilon0 * epsilon_r[m] * (1 + sigma[m]*dt / (2*epsilon0*epsilon_r[m])))) * Jz

    # 3. Conditions limites absorbantes (éq. 17 du projet)
    #    Valables uniquement si les bords sont en espace libre
    E[0]   = E2[1]
    E[M-1] = E2[M-2]

    # Décalage mémoire
    E2[:] = E1
    E1[:] = E

    E_enregistré[q, :] = E
    

# ── Amplitude du phaseur : (max - min) / 2 sur une période en régime permanent
pas_T        = round(1 / (f * dt))  # pas par période, T/dt donc 1/f*dt
début_régime = Q - pas_T    # moment fixé en régime permanent
E_permanent = E_enregistré[début_régime: début_régime + pas_T , :]
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
