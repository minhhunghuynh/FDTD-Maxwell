# Projet ELEC H304 — Méthode FDTD pour les équations de Maxwell

> **École Polytechnique de Bruxelles — Physique des télécommunications**
> Année académique 2025–2026
>
> **Auteurs :** HUYNH Minhhung · DE DONCKER Phillipe

---

## Table des matières

1. [Description du projet](#description-du-projet)
2. [Fondements théoriques](#fondements-théoriques)
3. [Structure du dépôt](#structure-du-dépôt)
4. [Étape 1 — FDTD 1D en espace libre](#étape-1--fdtd-1d-en-espace-libre)
5. [Étape 2 — FDTD 1D dans un diélectrique](#étape-2--fdtd-1d-dans-un-diélectrique)
6. [Étape 3 — FDTD 2D en espace libre](#étape-3--fdtd-2d-en-espace-libre)
7. [Étape 4 — FDTD 2D dans un diélectrique](#étape-4--fdtd-2d-dans-un-diélectrique)
8. [Simulations supplémentaires](#simulations-supplémentaires)
9. [Prérequis et installation](#prérequis-et-installation)
10. [Lancer les simulations](#lancer-les-simulations)
11. [Résultats et figures](#résultats-et-figures)
12. [Références](#références)

---

## Description du projet

Ce dépôt contient l'implémentation complète de la méthode **FDTD** (*Finite-Difference Time-Domain*, ou Différences Finies dans le Domaine Temporel) pour résoudre numériquement les **équations de Maxwell** en régime transitoire.

Le projet a été réalisé dans le cadre du cours **ELEC H304 — Physique des télécommunications** à l'École Polytechnique de Bruxelles. Il couvre progressivement :

- la propagation d'ondes électromagnétiques en **1D** (espace libre puis milieu diélectrique avec et sans pertes),
- la propagation en **2D** (espace libre, puis milieu diélectrique, conducteur parfait, obstacle diélectrique),
- des **simulations supplémentaires** illustrant la loi de Snell-Descartes et les franges d'interférence (expérience de Young).

---

## Fondements théoriques

### Équations de Maxwell

Les équations de Maxwell en régime temporel (avec source de courant) s'écrivent :

$$\nabla \times \vec{E} = -\frac{\partial \vec{B}}{\partial t}$$

$$\nabla \times \vec{B} = \varepsilon_0 \mu_0 \frac{\partial \vec{E}}{\partial t} + \mu_0 \vec{J}^{\text{source}}$$

En présence d'un milieu diélectrique conducteur (permittivité relative $\varepsilon_r$, conductivité $\sigma$) :

$$\nabla \times \vec{B} = \varepsilon \mu_0 \frac{\partial \vec{E}}{\partial t} + \mu_0 \sigma \vec{E} + \mu_0 \vec{J}^{\text{source}}$$

### Schéma de Yee (algorithme FDTD)

L'algorithme FDTD repose sur le **schéma de Yee** : les champs électrique $E$ et magnétique $B$ sont décalés d'un demi-pas à la fois dans l'espace et dans le temps (*leap-frog*). Les dérivées partielles sont approchées par des **différences finies centrées d'ordre 2**.

#### Condition de stabilité (critère CFL)

$$\Delta t = \frac{\Delta x}{2c} \quad \text{(1D)}, \qquad \Delta t \leq \frac{\Delta x}{2\sqrt{2}\,c} \quad \text{(2D)}$$

avec $c = 1/\sqrt{\varepsilon_0 \mu_0}$ la vitesse de la lumière dans le vide.

### Types de sources utilisées

- **Source sinusoïdale :** $J_z^q[m_s] = \sin(2\pi f\,q\Delta t)$, placée au milieu de la grille, génère une onde monochromatique de fréquence $f$.
- **Source gaussienne :** $J_z^q[m_s] = \exp\!\left(-\frac{(q\Delta t - t_0)^2}{2\tau^2}\right)$, impulsion large bande permettant d'observer le comportement dispersif des milieux.

---

## Structure du dépôt

```
FDTD-Maxwell/
│
├── ETAPE1/                        # FDTD 1D — espace libre
│   ├── sinusoïdale.py             # Source sinusoïdale
│   ├── gaussienne.py              # Source gaussienne (impulsion)
│   ├── conditionsnulles.py        # Conditions aux bords nulles
│   └── phaseur.py                 # Calcul de l'amplitude du phaseur
│
├── ETAPE2/                        # FDTD 1D — milieu diélectrique
│   ├── partie1.py                 # Diélectrique sans pertes
│   ├── partie2.py                 # Diélectrique avec pertes (conductivité σ)
│   ├── phaseur.py                 # Phaseur en milieu diélectrique
│   ├── comparaison_interface.png
│   ├── phase_reflexion.png
│   ├── sigma_comparison.png
│   └── phaseur_debug.png
│
├── ETAPE3/                        # FDTD 2D — espace libre
│   ├── partie1.py                 # Propagation 2D, source ponctuelle
│   ├── partie3.py
│   └── phaseur.py                 # Phaseur 2D
│
├── ETAPE4/                        # FDTD 2D — milieu diélectrique / obstacles
│   ├── partie1.py                 # Obstacle conducteur parfait
│   └── partie2.py                 # Obstacle diélectrique
│
├── supplémentaire/                # Simulations supplémentaires
│   ├── loidesnell.py              # Vérification de la loi de Snell-Descartes
│   ├── Frange.py                  # Franges d'interférence (Young)
│   └── outputs/                   # Figures générées
│
├── ORAL/                          # Fichiers pour la présentation orale
│   ├── chapitre4.py               # Guide d'ondes entre deux plaques métalliques
│   ├── chapitre9/                 # Simulations complémentaires
│   ├── fdtd_animation.gif         # Animation de la propagation FDTD
│   ├── fdtd_guide_plaques.png     # Vue spatiale du guide d'ondes
│   └── fdtd_profils_transverses.png
│
├── image projet/                  # Images et illustrations du projet
│
└── rapportFDTD.pdf                # Rapport complet du projet
```

---

## Étape 1 — FDTD 1D en espace libre

**Fichiers :** `ETAPE1/`

On suppose le champ électrique polarisé selon $\hat{z}$ et se propageant selon $x$. Les équations de Maxwell se réduisent à un système couplant $E_z$ et $\tilde{B}_y = cB_y$. En posant $\Delta t = \Delta x/(2c)$, les équations discrètes de mise à jour sont :

$$E_z^{q+1}[m] = E_z^q[m] + \frac{1}{2}\left(\tilde{B}_y^{q+1/2}\!\left[m+\tfrac{1}{2}\right] - \tilde{B}_y^{q+1/2}\!\left[m-\tfrac{1}{2}\right]\right) - \frac{\Delta t}{\varepsilon_0}J_z^{q+1/2}[m]$$

$$\tilde{B}_y^{q+1/2}\!\left[m+\tfrac{1}{2}\right] = \tilde{B}_y^{q-1/2}\!\left[m+\tfrac{1}{2}\right] + \frac{1}{2}\left(E_z^q[m+1] - E_z^q[m]\right)$$

### Simulations incluses

| Script | Description |
|--------|-------------|
| `sinusoïdale.py` | Onde plane harmonique dans le vide, animation temps réel |
| `gaussienne.py` | Impulsion gaussienne large bande |
| `conditionsnulles.py` | Réflexion aux bords (conditions $E=0$) |
| `phaseur.py` | Extraction de l'amplitude locale $\|E_z(x)\|$ en régime établi |

### Paramètres de discrétisation

| Paramètre | Valeur |
|-----------|--------|
| Pas spatial $\Delta x$ | $\lambda / 20$ |
| Pas temporel $\Delta t$ | $\Delta x / (2c)$ |
| Position de la source $m_s$ | Milieu de la grille |

### Résultats

- L'amplitude du phaseur est **constante** en espace libre (≈ 3 V/m), confirmant l'absence d'atténuation.
- La source gaussienne montre clairement la propagation bidirectionnelle de l'impulsion.
- Les réflexions aux bords produisent des ondes stationnaires visibles.

---

## Étape 2 — FDTD 1D dans un diélectrique

**Fichiers :** `ETAPE2/`

### Sans pertes ($\sigma = 0$)

La permittivité $\varepsilon = \varepsilon_r \varepsilon_0$ est introduite. L'onde transmise dans le diélectrique a :
- une **longueur d'onde réduite** $\lambda/\sqrt{\varepsilon_r}$,
- une **vitesse de phase** $v_\phi = c/\sqrt{\varepsilon_r}$.

Des **ondes stationnaires** sont visibles dans le milieu 1 (avant l'interface), résultant de la superposition des ondes incidente et réfléchie.

### Avec pertes ($\sigma > 0$)

On généralise les équations en ajoutant la conductivité du milieu. L'équation de mise à jour de $E_z$ devient :

$$E_z^{q+1}[m] = \left(1 - \frac{\sigma\,\Delta t}{\varepsilon_r\varepsilon_0}\right)E_z^q[m] + \frac{1}{2\varepsilon_r}\left(\tilde{B}_y^{q+1/2}\!\left[m+\tfrac{1}{2}\right] - \tilde{B}_y^{q+1/2}\!\left[m-\tfrac{1}{2}\right]\right)$$

Le facteur $(1 - \sigma\Delta t/\varepsilon_r\varepsilon_0)$ dissipe l'énergie à chaque pas de temps. L'amplitude décroît exponentiellement dans le diélectrique avec pertes.

### Analyse du phaseur

Le champ total dans la zone avant l'interface est la superposition de l'onde incidente et de l'onde réfléchie :

$$\underline{E}(x) = A_i e^{-j\beta x} + A_r e^{+j\beta x}$$

L'amplitude du phaseur $|\underline{E}(x)|$ présente des oscillations dont l'enveloppe permet d'extraire le coefficient de réflexion.

---

## Étape 3 — FDTD 2D en espace libre

**Fichiers :** `ETAPE3/`

### Algorithme de Yee en 2D (mode TM$_z$)

On cherche les solutions avec $\vec{E} = E_z(x,y,t)\,\hat{z}$. La grille 2D de Yee place :
- $E_z$ aux points entiers $(m\Delta x,\, n\Delta y)$,
- $B_x$ aux points $(m\Delta x,\, (n+\tfrac{1}{2})\Delta y)$,
- $B_y$ aux points $((m+\tfrac{1}{2})\Delta x,\, n\Delta y)$.

Les équations de mise à jour discrètes sont :

$$E_z^{q+1}[m,n] = E_z^q[m,n] - \frac{\Delta t}{\varepsilon_0}J_z^{q+1/2}[m,n] + \frac{1}{2}\left(\tilde{B}_y^{q+1/2}\!\left[m+\tfrac{1}{2},n\right] - \tilde{B}_y^{q+1/2}\!\left[m-\tfrac{1}{2},n\right]\right) - \frac{1}{2}\left(\tilde{B}_x^{q+1/2}\!\left[m,n+\tfrac{1}{2}\right] - \tilde{B}_x^{q+1/2}\!\left[m,n-\tfrac{1}{2}\right]\right)$$

$$\tilde{B}_x^{q+1/2}\!\left[m,n+\tfrac{1}{2}\right] = \tilde{B}_x^{q-1/2}\!\left[m,n+\tfrac{1}{2}\right] - \frac{1}{2}\left(E_z^q[m,n+1] - E_z^q[m,n]\right)$$

$$\tilde{B}_y^{q+1/2}\!\left[m+\tfrac{1}{2},n\right] = \tilde{B}_y^{q-1/2}\!\left[m+\tfrac{1}{2},n\right] + \frac{1}{2}\left(E_z^q[m+1,n] - E_z^q[m,n]\right)$$

### Décroissance du champ en 2D

En 2D, par conservation de l'énergie (théorème de Poynting), le vecteur de Poynting vérifie $S \cdot 2\pi r = \text{const}$, ce qui implique :

$$S \propto \frac{1}{r} \quad \Rightarrow \quad |\vec{E}| \propto \frac{1}{\sqrt{r}}$$

Le phaseur 2D confirme cette décroissance radiale observée dans les simulations.

---

## Étape 4 — FDTD 2D dans un diélectrique

**Fichiers :** `ETAPE4/`

### Conducteur parfait (`partie1.py`)

Un obstacle rectangulaire conducteur parfait est placé aux trois quarts de la grille en $x$ et au milieu en $y$. Les conditions aux limites imposent $E_z = 0$ à l'intérieur et sur les bords du conducteur. Les points sur les faces de l'obstacle rayonnent comme des sources secondaires (principe de Huygens), produisant :
- des **interférences constructives** en avant de l'obstacle,
- une **zone d'ombre** partielle derrière l'obstacle,
- une réduction de l'intensité proportionnelle à la taille de l'obstacle.

### Obstacle diélectrique (`partie2.py`)

L'équation de mise à jour 2D avec diélectrique (permittivité $\varepsilon_r$ et conductivité $\sigma$) est :

$$E_z^{q+1}[m,n] = \left(1 - \frac{\sigma\,\Delta t}{\varepsilon_r\varepsilon_0}\right)E_z^q[m,n] + \frac{1}{2\varepsilon_r}\left[\tilde{B}_y^{q+1/2}\!\left[m+\tfrac{1}{2},n\right] - \tilde{B}_y^{q+1/2}\!\left[m-\tfrac{1}{2},n\right] - \tilde{B}_x^{q+1/2}\!\left[m,n+\tfrac{1}{2}\right] + \tilde{B}_x^{q+1/2}\!\left[m,n-\tfrac{1}{2}\right]\right]$$

On retrouve l'espace libre pour $\varepsilon_r = 1$ et $\sigma = 0$.

---

## Simulations supplémentaires

**Fichiers :** `supplémentaire/`

### Loi de Snell-Descartes (`loidesnell.py`)

Plusieurs sources sinusoïdales alignées créent un front d'onde plan arrivant sur l'interface diélectrique avec un angle d'incidence $\theta_i = 30°$ (milieu 1 : $\varepsilon_{r1} = 2$).

La loi de Snell prédit l'angle de réfraction $\theta_t$ dans le milieu 2 :

$$n_1 \sin\theta_i = n_2 \sin\theta_t$$

La simulation vérifie numériquement cette relation et illustre le changement de direction de l'onde à l'interface.

### Franges d'interférence — Expérience de Young (`Frange.py`)

Deux sources ponctuelles sinusoïdales **en phase**, séparées d'une distance $d = 3\lambda$, produisent des franges d'interférence dans le plan 2D.

Les maxima d'interférence constructive sont donnés par :

$$d \sin\theta = n\lambda, \quad n \in \mathbb{Z}$$

Pour $d = 3\lambda$, la contrainte physique $|\sin\theta| \leq 1$ impose $|n| \leq 3$, donnant **7 maxima** aux angles :

$$\theta = \arcsin\!\left(\frac{n}{3}\right), \quad n \in \{-3, -2, -1, 0, 1, 2, 3\}$$

Les ordres limites $n = \pm 3$ correspondent à $\theta = \pm 90°$ (propagation le long de l'axe des sources).

---

## Prérequis et installation

### Dépendances Python

```bash
pip install numpy matplotlib
```

| Bibliothèque | Version recommandée | Usage |
|---|---|---|
| `numpy` | ≥ 1.21 | Calculs numériques, tableaux FDTD |
| `matplotlib` | ≥ 3.4 | Tracé des champs, animations |
| `matplotlib.animation` | (inclus dans matplotlib) | Animations FDTD en temps réel |

### Cloner le dépôt

```bash
git clone https://github.com/minhhunghuynh/FDTD-Maxwell.git
cd FDTD-Maxwell
```

---

## Lancer les simulations

### Étape 1 — 1D espace libre

```bash
python ETAPE1/sinusoïdale.py      # Onde sinusoïdale
python ETAPE1/gaussienne.py       # Impulsion gaussienne
python ETAPE1/conditionsnulles.py # Conditions aux bords
python ETAPE1/phaseur.py          # Analyse du phaseur
```

### Étape 2 — 1D diélectrique

```bash
python ETAPE2/partie1.py   # Diélectrique sans pertes
python ETAPE2/partie2.py   # Diélectrique avec pertes
python ETAPE2/phaseur.py   # Analyse du phaseur
```

### Étape 3 — 2D espace libre

```bash
python ETAPE3/partie1.py   # Propagation 2D
python ETAPE3/phaseur.py   # Phaseur 2D
```

### Étape 4 — 2D avec obstacles

```bash
python ETAPE4/partie1.py   # Conducteur parfait
python ETAPE4/partie2.py   # Obstacle diélectrique
```

### Simulations supplémentaires

```bash
python supplémentaire/loidesnell.py   # Loi de Snell-Descartes
python supplémentaire/Frange.py       # Franges de Young
```

---

## Résultats et figures

### Tableau récapitulatif des figures

| Étape | Figure | Phénomène illustré |
|-------|--------|-------------------|
| 1.1 | Fig. 1.1 | Propagation sinusoïdale dans le vide |
| 1.2 | Fig. 1.2 | Réflexion aux bords (conditions nulles) |
| 1.4–1.5 | Fig. 1.4–1.5 | Amplitude du phaseur constante ≈ 3 V/m |
| 2.1 | Fig. 1.6 | Onde dans un diélectrique sans pertes |
| 2.2 | Fig. 1.8 | Phaseur avec atténuation diélectrique |
| 3 | Fig. 1.10 | Décroissance en $1/\sqrt{r}$ en 2D |
| 4.1 | Fig. 1.15 | Diffraction autour d'un conducteur parfait |
| 4.2 | Fig. 1.16 | Transmission dans un obstacle diélectrique |
| Suppl. | Fig. 2.1 | Réfraction — vérification de Snell |
| Suppl. | Fig. 2.2 | Franges d'interférence (Young) |
| Suppl. | Fig. 2.3 | Phaseur des franges avec ordres annotés |

### Phénomènes physiques couverts

- **Propagation** : vitesse de phase $v_\phi = c/\sqrt{\varepsilon_r}$ dans un diélectrique.
- **Réflexion / Transmission** : coefficients de Fresnel à une interface diélectrique.
- **Atténuation** : décroissance exponentielle dans un milieu à pertes ($\sigma > 0$).
- **Diffraction** : contournement d'un conducteur parfait et recomposition du front d'onde.
- **Interférences** : franges constructives/destructives, relation $d\sin\theta = n\lambda$.
- **Rayonnement 2D** : décroissance $|E| \propto 1/\sqrt{r}$ (vs $1/r$ en 3D).
- **Réfraction** : vérification numérique de la loi de Snell-Descartes.

---

## Références

- **[DED26]** DE DONCKER Phillipe, *Notes de cours ELEC H304 — Physique des télécommunications*, École Polytechnique de Bruxelles, 2025–2026.
- **A. Taflove & S. Hagness**, *Computational Electrodynamics: The Finite-Difference Time-Domain Method*, Artech House, 3rd ed., 2005.
- **K. S. Yee**, "Numerical solution of initial boundary value problems involving Maxwell's equations in isotropic media," *IEEE Trans. Antennas Propag.*, vol. 14, no. 3, pp. 302–307, 1966.

---

*Rapport complet disponible dans* [`rapportFDTD.pdf`](rapportFDTD.pdf).
