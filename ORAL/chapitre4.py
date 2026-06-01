#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║   Dipôle λ/2 — Animation interactive de RÉCEPTION                          ║
║   Chapitre 4 : Émission et réception des ondes EM                          ║
║                                                                              ║
║   Dépendances : numpy, matplotlib                                            ║
║   Lancement   : python dipole_reception.py                                  ║
╚══════════════════════════════════════════════════════════════════════════════╝

Paramètres jouables (sliders) :
  - θ        : angle d'incidence depuis l'axe du dipôle  (0°=axe → 90°=max)
  - Vitesse  : vitesse de défilement de l'onde incidente
  - E₀       : amplitude du champ électrique incident
  - λ visuel : longueur d'onde visuelle (espacement des fronts)

Boutons :
  - Play / Pause      : démarrer ou figer l'animation
  - Diagramme sin²θ  : superposer le diagramme de réception
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as mpatches
from matplotlib.widgets import Slider, Button

# ══════════════════════════════════════════════════════════════════════════════
#  Paramètres par défaut (modifiables ici ou via sliders)
# ══════════════════════════════════════════════════════════════════════════════
DEFAULT_THETA    = 45.0   # angle initial θ en degrés (depuis axe dipôle)
DEFAULT_SPEED    = 1.0    # vitesse de l'animation
DEFAULT_E0       = 1.5    # amplitude du champ E incident
DEFAULT_LAMBDA_V = 1.8    # longueur d'onde visuelle (unités axes)

DIPOLE_LEN = 1.10         # demi-longueur du dipôle
DIPOLE_GAP = 0.10         # demi-écart au point d'alimentation
XLIM       = (-4.4, 6.8)
YLIM       = (-3.3, 3.3)

# Couleurs
C_WAVE    = '#7F77DD'   # fronts d'onde (violet)
C_EFIELD  = '#D85A30'   # champ E⃗ (corail)
C_CURRENT = '#1D9E75'   # courant induit (vert)
C_PATTERN = '#1D9E75'   # diagramme sin²θ
C_DIPOLE  = '#3C3C3A'   # corps de l'antenne
C_BG      = '#FAFAF7'   # fond

# ══════════════════════════════════════════════════════════════════════════════
#  État interne de l'animation
# ══════════════════════════════════════════════════════════════════════════════
state = {
    'phase'       : 0.0,
    'playing'     : False,
    'show_pattern': False,
}

# ══════════════════════════════════════════════════════════════════════════════
#  Mise en page
# ══════════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(13, 8.5), facecolor=C_BG)
fig.suptitle('Antenne dipôle λ/2 — Réception  ·  Chapitre 4',
             fontsize=14, color='#3C3C3A', y=0.985, fontweight='normal')

# Axes animation (gauche)
ax = fig.add_axes([0.03, 0.23, 0.62, 0.73], facecolor=C_BG)
ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
ax.set_aspect('equal'); ax.axis('off')

# Panneau info (droite)
ax_info = fig.add_axes([0.67, 0.23, 0.31, 0.73], facecolor='#F0EEE8')
ax_info.set_xlim(0, 1); ax_info.set_ylim(0, 1); ax_info.axis('off')

# ── Sliders ──
ax_s_theta  = fig.add_axes([0.08, 0.155, 0.53, 0.022])
ax_s_speed  = fig.add_axes([0.08, 0.105, 0.53, 0.022])
ax_s_e0     = fig.add_axes([0.08, 0.055, 0.53, 0.022])
ax_s_lambda = fig.add_axes([0.08, 0.005, 0.53, 0.022])

s_theta  = Slider(ax_s_theta,  'θ (° depuis axe)', 0.5,  89.5, valinit=DEFAULT_THETA,
                  color=C_WAVE,    initcolor='none')
s_speed  = Slider(ax_s_speed,  'Vitesse          ', 0.2,  5.0,  valinit=DEFAULT_SPEED,
                  color=C_CURRENT, initcolor='none')
s_e0     = Slider(ax_s_e0,     'Amplitude E₀     ', 0.3,  3.5,  valinit=DEFAULT_E0,
                  color=C_EFIELD,  initcolor='none')
s_lambda = Slider(ax_s_lambda, 'λ visuel         ', 0.8,  3.5,  valinit=DEFAULT_LAMBDA_V,
                  color='#888780', initcolor='none')

# ── Boutons ──
ax_b_play    = fig.add_axes([0.68, 0.155, 0.13, 0.05])
ax_b_pattern = fig.add_axes([0.84, 0.155, 0.14, 0.05])

b_play    = Button(ax_b_play,    '▶  Play',    color='#E8E6DC', hovercolor='#D3D1C7')
b_pattern = Button(ax_b_pattern, 'sin²θ',      color='#E8E6DC', hovercolor='#D3D1C7')

# ══════════════════════════════════════════════════════════════════════════════
#  Fonctions physiques
# ══════════════════════════════════════════════════════════════════════════════

def get_theta():
    """θ en radians (angle depuis l'axe du dipôle)."""
    return s_theta.val * np.pi / 180.0

def get_efficiency():
    """Efficacité de réception normalisée = sin²θ."""
    return np.sin(get_theta()) ** 2

def get_lambda():
    return s_lambda.val

def prop_dir():
    """Direction de propagation de l'onde (vers le dipôle)."""
    t = get_theta()
    # θ=90°: source à droite  → prop vers gauche  = (-1,  0)
    # θ= 0°: source au-dessus → prop vers le bas  = ( 0, -1)
    return np.array([-np.sin(t), -np.cos(t)])

def perp_dir():
    """Vecteur perpendiculaire aux fronts d'onde (parallèle aux fronts)."""
    p = prop_dir()
    return np.array([-p[1], p[0]])   # rotation 90° CCW

def efield_dir():
    """Direction du champ E⃗ (perpendiculaire à la propagation)."""
    # Pour une onde plane en 2D, E ⊥ β → même que perp_dir
    return perp_dir()

# ══════════════════════════════════════════════════════════════════════════════
#  Fonctions de dessin
# ══════════════════════════════════════════════════════════════════════════════

def draw_wave_fronts(phase):
    """Fronts d'onde planes animés avec flèches E⃗."""
    lv   = get_lambda()
    p    = prop_dir()          # direction de propagation
    pv   = perp_dir()          # direction parallèle aux fronts
    ed   = efield_dir()        # direction E⃗
    E0   = s_e0.val
    osc  = np.sin(phase / lv * 2 * np.pi)   # oscillation temporelle
    T    = 12.0                # longueur de trait (plus longue que la figure)

    for k in range(-3, 14):
        # Distance signée du front par rapport au centre (>0 = approchant)
        sd = k * lv - phase
        if sd > lv * 6.0 or sd < -lv * 1.3:
            continue

        # Transparence : fondu en entrant et en sortant
        alpha = 0.82
        if sd > lv * 2.0:
            alpha *= max(0.0, 1.0 - (sd - lv * 2.0) / (lv * 3.0))
        elif sd < 0.0:
            alpha = max(0.0, alpha * (1.0 + sd / lv))
        if alpha < 0.02:
            continue

        # Centre du front
        fx, fy = sd * p[0], sd * p[1]

        # ── Trait du front ──
        ax.plot([fx - T*pv[0], fx + T*pv[0]],
                [fy - T*pv[1], fy + T*pv[1]],
                '-', color=C_WAVE, alpha=alpha, linewidth=1.8, zorder=2)

        # ── Flèches E⃗ sur les fronts approchants ──
        if 0.0 < sd < lv * 4.0 and k % 2 == 0:
            e_amp = E0 * 0.32 * osc * alpha
            if abs(e_amp) > 0.03:
                t_positions = np.linspace(-2.8, 2.8, 6)
                xs = fx + t_positions * pv[0]
                ys = fy + t_positions * pv[1]
                # Filtrer hors-cadre
                mask = ((xs > XLIM[0]+0.4) & (xs < XLIM[1]-0.4) &
                        (ys > YLIM[0]+0.4) & (ys < YLIM[1]-0.4))
                xs, ys = xs[mask], ys[mask]
                if len(xs):
                    ax.quiver(xs, ys,
                              np.full(len(xs), ed[0] * e_amp),
                              np.full(len(xs), ed[1] * e_amp),
                              color=C_EFIELD, alpha=min(0.95, alpha * 0.90),
                              scale=1, scale_units='xy', angles='xy',
                              headwidth=4, headlength=5, headaxislength=3.5,
                              width=0.007, zorder=3)


def draw_dipole(phase):
    """Dipôle avec halo d'activité et courant induit animé."""
    eff = get_efficiency()
    lv  = get_lambda()
    osc = np.sin(phase / lv * 2 * np.pi)

    # ── Halo (intensité visuelle du courant) ──
    if eff > 0.04:
        for lw_h, al_h in [(18, 0.04*eff), (12, 0.09*eff), (7, 0.14*eff)]:
            ax.plot([0, 0], [-DIPOLE_LEN, DIPOLE_LEN], '-',
                    color=C_CURRENT, lw=lw_h, alpha=al_h, zorder=4,
                    solid_capstyle='round')

    # ── Corps du dipôle ──
    ax.plot([0, 0], [ DIPOLE_GAP,  DIPOLE_LEN], '-', color=C_DIPOLE,
            lw=6, solid_capstyle='round', zorder=5)
    ax.plot([0, 0], [-DIPOLE_GAP, -DIPOLE_LEN], '-', color=C_DIPOLE,
            lw=6, solid_capstyle='round', zorder=5)

    # ── Point d'alimentation ──
    ax.plot(0, 0, 'o', color=C_CURRENT, markersize=10, zorder=6)

    # ── Flèches de courant induit (amplitude ∝ sin²θ, oscillent avec l'onde) ──
    c_amp = eff * 0.55 * osc
    y_positions = np.linspace(-DIPOLE_LEN * 0.80, DIPOLE_LEN * 0.80, 8)
    for yp in y_positions:
        if abs(yp) < DIPOLE_GAP * 0.6:
            continue
        # Distribution sinusoïdale du courant le long du dipôle λ/2
        i_local = c_amp * np.cos(yp / DIPOLE_LEN * np.pi / 2)
        if abs(i_local) > 0.04:
            ax.quiver(0, yp, 0, i_local,
                      color=C_CURRENT,
                      alpha=min(0.92, eff * 1.3),
                      scale=1, scale_units='xy', angles='xy',
                      headwidth=5, headlength=6, headaxislength=4.5,
                      width=0.012, zorder=7)

    # ── Ligne de transmission symbolique ──
    for dx_tl in [0.0, 0.14]:
        ax.plot([dx_tl, dx_tl],
                [-DIPOLE_LEN - 0.05, -DIPOLE_LEN - 0.75],
                '--', color='#AAAAAA', lw=1.2, zorder=3)
    ax.text(0.07, -DIPOLE_LEN - 0.92, 'R  (récepteur)',
            ha='center', color='#AAAAAA', fontsize=8.5, style='italic')

    # ── Étiquette V_oc ──
    voc_pct = np.sqrt(eff) * 100
    col = (C_CURRENT if eff > 0.50 else ('#EF9F27' if eff > 0.20 else '#E24B4A'))
    ax.text(0.22, 0.0,
            f'$V_{{oc}}$ ∝ {voc_pct:.0f}%',
            color=col, fontsize=10.5, va='center', zorder=8)


def draw_reception_pattern():
    """Diagramme de réception sin²θ en superposition."""
    theta_arr = np.linspace(0, 2 * np.pi, 500)
    # sin²θ_physique = sin²(angle depuis axe y) = cos²(angle standard canvas)
    scale = min(abs(XLIM[0]), YLIM[1]) * 0.68
    r_arr = scale * np.cos(theta_arr) ** 2
    xp    = r_arr * np.cos(theta_arr)
    yp    = r_arr * np.sin(theta_arr)
    ax.fill(xp, yp, color=C_PATTERN, alpha=0.13, zorder=1)
    ax.plot(xp, yp, '--', color=C_PATTERN, alpha=0.55, lw=1.6, zorder=1)
    ax.text(scale * 1.03, 0.12, 'sin²θ',
            color=C_PATTERN, fontsize=9, alpha=0.80)


def draw_labels_and_source():
    """Axes du dipôle, arc θ et indicateur de source."""
    t = get_theta()

    # Axes du dipôle (vertical)
    for y_lbl, txt in [(YLIM[1]-0.20, 'axe  θ=0°'),
                       (YLIM[0]+0.20, 'axe  θ=0°')]:
        ax.text(0, y_lbl, txt, ha='center', va='center',
                color='#CCCCCC', fontsize=8.5)
    ax.text(XLIM[0]+0.12, 0, 'θ = 90°\n(max)', ha='left',
            va='center', color='#CCCCCC', fontsize=8.5)

    # Indicateur de source (flèche pointillée depuis la source)
    src_r = 3.6
    src_x =  src_r * np.sin(t)    # position source : (sin θ, cos θ) depuis le centre
    src_y =  src_r * np.cos(t)
    mid_x =  0.55 * src_x
    mid_y =  0.55 * src_y
    ax.annotate('', xy=(mid_x, mid_y), xytext=(src_x, src_y),
                arrowprops=dict(arrowstyle='->', color=C_WAVE, lw=1.3,
                                alpha=0.60, linestyle='dashed'), zorder=2)
    ax.text(src_x, src_y + (0.28 if src_y >= 0 else -0.38),
            f'Source\n(θ = {s_theta.val:.0f}°)',
            ha='center', color=C_WAVE, fontsize=8.5, alpha=0.90,
            va='bottom' if src_y >= 0 else 'top')

    # Arc indiquant θ
    if t > 0.08:
        arc_r = 0.60
        arc_t = np.linspace(np.pi/2, np.pi/2 + t, 40)
        ax.plot(arc_r * np.cos(arc_t), arc_r * np.sin(arc_t),
                '-', color='#AAAAAA', lw=1.0, alpha=0.70, zorder=2)
        mid_arc = np.pi/2 + t / 2
        ax.text(0.88 * np.cos(mid_arc), 0.88 * np.sin(mid_arc), 'θ',
                ha='center', va='center', color='#AAAAAA', fontsize=10)


def draw_info_panel():
    """Panneau récapitulatif affiché à droite."""
    ax_info.clear()
    ax_info.set_xlim(0, 1); ax_info.set_ylim(0, 1); ax_info.axis('off')

    eff       = get_efficiency()
    theta_deg = s_theta.val
    voc_pct   = np.sqrt(eff) * 100

    # ── Titre ──
    ax_info.text(0.50, 0.965, 'Physique', ha='center', fontsize=12,
                 color='#3C3C3A', fontweight='bold')
    ax_info.axhline(0.935, xmin=0.04, xmax=0.96, color='#D3D1C7', lw=0.8)

    # ── Formules ──
    ax_info.text(0.07, 0.895,
                 r'$V_{oc} = -\vec{h}_e(\theta,\phi)\cdot\vec{E}_i$',
                 fontsize=10, color='#3C3C3A')
    ax_info.text(0.07, 0.825,
                 r'$P_L \propto |\vec{h}_e|^2 \propto \sin^2\!\theta$',
                 fontsize=10, color='#3C3C3A')
    ax_info.axhline(0.785, xmin=0.04, xmax=0.96, color='#D3D1C7', lw=0.8)

    # ── Valeurs courantes ──
    ax_info.text(0.07, 0.750, f'θ = {theta_deg:.1f}°  (depuis axe)',
                 fontsize=10, color='#5F5E5A')
    ax_info.text(0.07, 0.695, f'sin²θ = {eff:.4f}',
                 fontsize=11, color='#3C3C3A', fontweight='bold')
    ax_info.text(0.07, 0.645, f'|h_e| ∝ sin({theta_deg:.0f}°) = {np.sin(get_theta()):.3f}',
                 fontsize=9,  color='#5F5E5A')
    ax_info.axhline(0.612, xmin=0.04, xmax=0.96, color='#D3D1C7', lw=0.8)

    # ── Barre de puissance ──
    ax_info.text(0.07, 0.580, 'Puissance reçue P_L :', fontsize=9.5, color='#5F5E5A')
    bx, by, bw, bh = 0.07, 0.510, 0.76, 0.052
    ax_info.add_patch(mpatches.FancyBboxPatch(
        (bx, by), bw, bh, boxstyle='round,pad=0.005',
        facecolor='#E8E6DC', edgecolor='#B4B2A9', lw=0.8))
    bar_col = (C_CURRENT if eff > 0.50 else ('#EF9F27' if eff > 0.20 else '#E24B4A'))
    if eff * bw > 0.006:
        ax_info.add_patch(mpatches.FancyBboxPatch(
            (bx, by), eff * bw, bh, boxstyle='round,pad=0.005',
            facecolor=bar_col, edgecolor='none'))
    ax_info.text(bx + bw + 0.03, by + bh/2,
                 f'{eff*100:.1f}%', fontsize=10.5,
                 color=bar_col, va='center', fontweight='bold')
    ax_info.axhline(0.46, xmin=0.04, xmax=0.96, color='#D3D1C7', lw=0.8)

    # ── Interprétation physique ──
    if theta_deg < 10:
        msg = "E⃗ ⊥ dipôle\n→ aucun couplage\nV_oc = 0"
        ic  = '#A32D2D'
    elif theta_deg > 80:
        msg = "E⃗ ∥ dipôle\n→ couplage maximal\nV_oc ∝ h_e · E₀"
        ic  = C_CURRENT
    else:
        msg = (f"Couplage partiel\n"
               f"E_z = E₀ · sin({theta_deg:.0f}°)")
        ic  = '#BA7517'

    ax_info.add_patch(mpatches.FancyBboxPatch(
        (0.07, 0.255), 0.86, 0.175, boxstyle='round,pad=0.012',
        facecolor='white', edgecolor=ic, lw=1.0, alpha=0.90))
    ax_info.text(0.50, 0.342, msg, ha='center', va='center',
                 fontsize=9, color=ic, style='italic', linespacing=1.55)

    # ── Légende couleurs ──
    ax_info.axhline(0.22, xmin=0.04, xmax=0.96, color='#D3D1C7', lw=0.8)
    ax_info.text(0.07, 0.175, '━  Onde incidente',   fontsize=8.5, color=C_WAVE)
    ax_info.text(0.07, 0.120, '━  Champ E⃗ (proj.)', fontsize=8.5, color=C_EFIELD)
    ax_info.text(0.55, 0.175, '━  Courant I_a',      fontsize=8.5, color=C_CURRENT)
    if state['show_pattern']:
        ax_info.text(0.55, 0.120, '╌╌  sin²θ',
                     fontsize=8.5, color=C_PATTERN)


# ══════════════════════════════════════════════════════════════════════════════
#  Boucle d'animation
# ══════════════════════════════════════════════════════════════════════════════

def animate(frame):
    ax.clear()
    ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_facecolor(C_BG)

    phase = state['phase']

    if state['show_pattern']:
        draw_reception_pattern()

    draw_wave_fronts(phase)
    draw_dipole(phase)
    draw_labels_and_source()
    draw_info_panel()

    if state['playing']:
        state['phase'] = (state['phase'] + 0.60 * s_speed.val) % get_lambda()

    return []


# ══════════════════════════════════════════════════════════════════════════════
#  Callbacks
# ══════════════════════════════════════════════════════════════════════════════

def toggle_play(event):
    state['playing'] = not state['playing']
    b_play.label.set_text('⏸  Pause' if state['playing'] else '▶  Play')
    fig.canvas.draw_idle()

def toggle_pattern(event):
    state['show_pattern'] = not state['show_pattern']
    b_pattern.color = '#D0ECD8' if state['show_pattern'] else '#E8E6DC'
    if not state['playing']:
        animate(0)
        fig.canvas.draw_idle()

def on_slider_change(val):
    if not state['playing']:
        animate(0)
        fig.canvas.draw_idle()

b_play.on_clicked(toggle_play)
b_pattern.on_clicked(toggle_pattern)
s_theta.on_changed(on_slider_change)
s_speed.on_changed(on_slider_change)
s_e0.on_changed(on_slider_change)
s_lambda.on_changed(on_slider_change)

# ══════════════════════════════════════════════════════════════════════════════
#  Lancement
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    ani = animation.FuncAnimation(
        fig, animate,
        interval=48,          # ~20 fps
        blit=False,
        cache_frame_data=False
    )
    plt.show()