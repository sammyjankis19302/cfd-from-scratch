"""Spectral-analysis headline figure.

The point this figure has to make, which the first version buried: BD-1, FD-1
and CD-2 have IDENTICAL real parts. Drawn as three separate curves, two are
painted over the third and the result disappears. Here they are drawn as one
curve, named for all three.
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    _HERE, "..", "02-spectral-analysis"))
from src.spectral import FIRST, group_velocity  # noqa: E402

OUT = sys.argv[2] if len(sys.argv) > 2 else _HERE
os.makedirs(OUT, exist_ok=True)

SURFACE = "#FFFFFF"
INK = "#1A1A1A"
INK_MUTED = "#5A5A5A"
GRID = "#D8D8D8"

# validated: all six checks pass
C = {"shared": "#2E7DB8", "CD-4": "#A8322A", "BD-1": "#2E7DB8",
     "FD-1": "#B8841A", "CD-2": "#4B4B8F"}

K = np.linspace(1e-6, np.pi, 2000)
Kp = K / np.pi


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ax.spines.values():
        s.set_color(GRID)
        s.set_linewidth(0.8)
    ax.tick_params(colors=INK_MUTED, labelsize=9, length=3, width=0.8)
    ax.grid(True, color=GRID, linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    ax.set_xlabel(r"$k\Delta x/\pi$", color=INK_MUTED, fontsize=10)
    ax.set_xlim(0, 1)


fig, ax = plt.subplots(1, 3, figsize=(14.2, 4.6), facecolor=SURFACE)

# --- panel 1: dispersion -------------------------------------------------
a = ax[0]
a.plot(Kp, np.ones_like(K), ls=(0, (4, 3)), color="#111111", lw=1.4, zorder=2)
a.plot(Kp, FIRST["CD-2"](K).real, color=C["shared"], lw=2.6, zorder=3)
a.plot(Kp, FIRST["CD-4"](K).real, color=C["CD-4"], lw=2.2, zorder=4)
a.text(0.04, 1.035, "exact", color="#111111", fontsize=9.5, va="bottom")
a.text(0.045, 0.40, "BD-1, FD-1, CD-2\nall identical here",
       color=C["shared"], fontsize=10, fontweight="bold", va="top")
a.text(0.63, 0.60, "CD-4", color=C["CD-4"], fontsize=10.5, fontweight="bold")
a.set_title(r"Re$(k_{num}/k)$  —  dispersion", color=INK, fontsize=11.5, pad=8)
a.set_ylim(-0.03, 1.12)
style(a)

# --- panel 2: dissipation ------------------------------------------------
a = ax[1]
a.axhline(0, color="#111111", lw=1.4, ls=(0, (4, 3)), zorder=2)
a.plot(Kp, FIRST["FD-1"](K).imag, color=C["FD-1"], lw=2.4, zorder=3)
a.plot(Kp, FIRST["BD-1"](K).imag, color=C["BD-1"], lw=2.4, zorder=3)
a.text(0.42, 0.52, "FD-1  (downwind)\namplifies → unstable",
       color=C["FD-1"], fontsize=10, fontweight="bold")
a.text(0.055, -0.56, "BD-1  (upwind)\ndamps → stable",
       color=C["BD-1"], fontsize=10, fontweight="bold", va="top")
a.text(0.04, 0.045, "CD-2, CD-4: exactly zero, at every wavenumber",
       color=INK_MUTED, fontsize=9.5)
a.set_title(r"Im$(k_{num}/k)$  —  dissipation", color=INK, fontsize=11.5, pad=8)
a.set_ylim(-0.95, 0.95)
style(a)

# --- panel 3: group velocity --------------------------------------------
a = ax[2]
a.axhspan(-2.0, 0, color="#A8322A", alpha=0.055, zorder=1)
a.axhline(1, color="#111111", lw=1.4, ls=(0, (4, 3)), zorder=2)
a.axhline(0, color="#8A8A8A", lw=1.0, zorder=2)
a.plot(Kp, group_velocity("CD-2", K), color=C["shared"], lw=2.4, zorder=3)
a.plot(Kp, group_velocity("CD-4", K), color=C["CD-4"], lw=2.2, zorder=3)
a.axvline(0.5, color="#8A8A8A", lw=1.0, ls=":", zorder=2)
a.text(0.515, -1.45, "4 points per\nwavelength", color=INK_MUTED, fontsize=9)
a.text(0.055, -0.42, "energy travels UPSTREAM", color="#A8322A",
       fontsize=10, fontweight="bold")
a.text(0.055, 1.06, "exact", color="#111111", fontsize=9.5)
a.text(0.18, 0.46, "CD-2", color=C["shared"], fontsize=10.5, fontweight="bold")
a.text(0.40, 0.72, "CD-4", color=C["CD-4"], fontsize=10.5, fontweight="bold")
a.set_title(r"$V_g/c$  —  group velocity", color=INK, fontsize=11.5, pad=8)
a.set_ylim(-1.85, 1.25)
style(a)

fig.suptitle(
    "What a convergence study cannot see:  two stencils of different order can need identical resolution,"
    "  and under-resolved waves can run backwards",
    color=INK, fontsize=12, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(os.path.join(OUT, "spectral.png"), dpi=170,
            facecolor=SURFACE, bbox_inches="tight")
plt.close(fig)
print("spectral.png")
