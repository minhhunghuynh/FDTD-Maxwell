import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Paramètres physiques
c = 3*10**8  # Vitesse de la lumière                    
mu0 = 4 * np.pi * 10**(-7)  # Perméabilité du vide      
epsilon0 = 1 / (mu0 * c**2)  # Permittivité du vide
f = 1e9  # fréquence de la source (1 GHz)   
L_onde = c / f   # Longueur d'onde

delta_x = L_onde / 20
delta_t = delta_x / (2 * c)

M = 400 # Nombre de points spatiaux
Q = 3000 # Nombre de points temporels

E = np.zeros(M)
B = np.zeros(M)


E_1  = np.zeros(M)  # E à q-1
E_2 = np.zeros(M)  # E à q-2

# Stockage de E à chaque pas pour l'animation
E_affichage = np.zeros((Q, M))

for q in range(Q):

    # actualisation de B et E
    for m in range(0, M - 1):
        B[m] = B[m] + 0.5 * (E[m+1] - E[m])
    for m in range(1, M - 1):
        J = 0
        if m == M // 2:
            t0    = 9e-9        
            sigma = 0.5e-9      # largeur de la gaussienne 
            J = np.sin(2 * np.pi * f * q * delta_t) # sin(2πft) 
            
        E[m] = E[m] + 0.5 * (B[m] - B[m-1]) - (delta_t / epsilon0) * J

    # 3. Conditions limites absorbantes (éq. 17)
    E[0]   = 0
    E[M-1] = 0

    # Sauvegarde après application des conditions limites
    E_2[:] = E_1
    E_1[:]  = E

    # Stockage pour l'animation
    E_affichage[q, :] = E

#  Animation 
x = np.arange(M) * delta_x / L_onde  # axe x en longueurs d'onde

fig, ax = plt.subplots(figsize=(10, 4))
line, = ax.plot(x, E_affichage[0, :], color='black')

ax.set_xlim(0, x[M-1])
ax.set_ylim(-1.5 * np.max(np.abs(E_affichage)), 1.5 * np.max(np.abs(E_affichage)))
ax.set_xlabel("Position x (en λ)")
ax.set_ylabel("Ez (V/m)")
ax.set_title("Propagation du champ électrique avec une source gaussienne")
ax.grid(True)

time_text = ax.text(0.02, 0.92, '', transform=ax.transAxes)

def update(q):
    line.set_ydata(E_affichage[q, :])
    time_text.set_text(f"t = {q * delta_t * 1e9:.2f} ns  ")
    return line, time_text

# On anime 1 frame sur 5 pour la fluidité
anim = FuncAnimation(fig, update, frames=range(0, Q, 5), interval=20, blit=True)

plt.tight_layout()
plt.show()