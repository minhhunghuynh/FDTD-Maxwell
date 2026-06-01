import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
# Paramètres Physiques et Numériques
c = 3*10**8  # Vitesse de la lumière                    
mu0 = 4 * np.pi * 10**(-7)  # Perméabilité du vide      
epsilon0 = 1 / (mu0 * c**2)  # Permittivité du vide
f = 1e9  # fréquence de la source (1 GHz)   
L = c / f    # Longueur d'onde              
delta_x = L / 20
delta_t = delta_x / (2 * c)
M = 400 # Nombre de points spatiaux
Q = 1000 # Nombre de points temporels
E = np.zeros(M)
B = np.zeros(M)
J = np.zeros(M)
E_1  = np.zeros(M)  # E à q-1
E_2 = np.zeros(M)  # E à q-2
m_source = M // 2  # position de la source (au centre du domaine)
# Stockage de E à chaque pas pour l'animation
E_affichage = np.zeros((Q, M))
for q in range(Q):
    #  Mise à jour de B et E
    for m in range(0, M - 1):
        B[m] = B[m] + 0.5 * (E[m+1] - E[m])
    for m in range(0, M - 1):
        #paramètres de la gaussienne
        t0    = 3e-9        
        sigma = 0.5e-9  
        jz = 0
        if m == m_source:
            jz = np.exp(-((q * delta_t - t0) ** 2) / (2 * sigma ** 2))       
        E[m] = E[m] + 0.5 * (B[m] - B[m-1]) - (delta_t / epsilon0) * jz
    # 3. Conditions limites absorbantes (éq. 17)
    E[0]   = E_2[1]
    E[M-1] = E_2[M-2]
    # Sauvegarde après application des conditions limites
    E_2[:] = E_1
    E_1[:]  = E
    
    # Stockage pour l'animation
    E_affichage[q, :] = E
# Animation 
x = np.arange(M) * delta_x # axe x 

fig, ax = plt.subplots(figsize=(10, 4))
line, = ax.plot(x, E_affichage[0, :], color='black')

ax.set_xlim(0, x[-1])
ax.set_ylim(-1.5 * np.max(np.abs(E_affichage)), 1.5 * np.max(np.abs(E_affichage)))
ax.set_xlabel("Position x (m)")
ax.set_ylabel("Ez (V/m)")
ax.set_title("Propagation du champ électrique avec une source gaussienne")
ax.grid(True)

time_text = ax.text(0.02, 0.92, '', transform=ax.transAxes)

def update(q):
    line.set_ydata(E_affichage[q, :])
    time_text.set_text(f"t = {q * delta_t * 1e9:.2f} ns")
    return line, time_text

anim = FuncAnimation(fig, update, frames=range(0, Q, 5), interval=20, blit=True)

plt.tight_layout()
plt.show()

print(f"delta_x = {delta_x:.4f} m")
print(f"L = {L:.4f} m")
print(f"largeur spatiale attendue = {c * 0.5e-9:.4f} m")