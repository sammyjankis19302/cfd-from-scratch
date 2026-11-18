"""Figures for module 01, generated from src/ — nothing hand-drawn.

    python3 make_figures.py

Two panels, because this module has two separate results:
  (a) dissipation and dispersion are different failures, visible on a square wave
  (b) the measured orders of accuracy
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from src.initial_conditions import exact_solution, gaussian, make_grid, square
from src.schemes import SCHEMES
from src.solver import l2_error, solve

OUT = os.path.join(_HERE, "figures")
os.makedirs(OUT, exist_ok=True)

SURFACE, INK, MUTED, GRID = "#FFFFFF", "#1A1A1A", "#5A5A5A", "#D8D8D8"
C_UP, C_LW, C_LF, C_BD = "#2E7DB8", "#A8322A", "#B8841A", "#4B4B8F"


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ax.spines.values():
        s.set_color(GRID); s.set_linewidth(0.8)
    ax.tick_params(colors=MUTED, labelsize=9, length=3, width=0.8)
    ax.xaxis.label.set_color(MUTED); ax.yaxis.label.set_color(MUTED)
    ax.grid(True, color=GRID, linewidth=0.6, alpha=0.7); ax.set_axisbelow(True)


fig, ax = plt.subplots(1, 2, figsize=(12.6, 4.7), facecolor=SURFACE)

# ---- (a) square wave after one full revolution --------------------------
N, Cfl, T, c, L = 200, 0.8, 1.0, 1.0, 1.0
x, dx = make_grid(N, L)

for step, col, lab in ((SCHEMES["FTBS"], C_UP, "FTBS (first-order upwind)"),
                       (SCHEMES["Lax-Wendroff"], C_LW, "Lax–Wendroff")):
    _, u, tf = solve(step, square, N, Cfl, T, c, L)
    ax[0].plot(x, u, color=col, lw=2.0, label=lab, zorder=3)

u_ex = exact_solution(square, x, L, c, tf)
ax[0].plot(x, u_ex, color="#111111", lw=1.4, ls=(0, (4, 3)),
           label="exact", zorder=4)

ax[0].annotate("corners rounded off\n= dissipation", xy=(0.553, 0.50),
               xytext=(0.645, 0.52), color=C_UP, fontsize=9.5, fontweight="bold",
               va="center",
               arrowprops=dict(arrowstyle="->", color=C_UP, lw=1.3))
ax[0].annotate("overshoot and undershoot\n= dispersion", xy=(0.262, -0.24),
               xytext=(0.46, -0.33), color=C_LW, fontsize=9.5, fontweight="bold",
               va="center",
               arrowprops=dict(arrowstyle="->", color=C_LW, lw=1.3))

ax[0].set_title(f"Square wave after one revolution  (N = {N}, C = {Cfl})",
                color=INK, fontsize=11.5, pad=8)
ax[0].set_xlabel("x"); ax[0].set_ylabel("u")
ax[0].set_ylim(-0.42, 1.42)
style(ax[0])
lg = ax[0].legend(frameon=False, fontsize=9, loc="upper right")
for t in lg.get_texts():
    t.set_color(MUTED)

# ---- (b) measured order of accuracy -------------------------------------
grids = (200, 400, 800, 1600)
series = (("FTBS", C_UP, 1), ("Lax-Friedrichs", C_LF, 1),
          ("Lax-Wendroff", C_LW, 2), ("BD-2", C_BD, 2))

curves = {}
for name, col, _ in series:
    errs = []
    for n in grids:
        xx, u, tf = solve(SCHEMES[name], gaussian, n, 0.5, 1.0, 1.0, 1.0, sigma=0.1)
        ue = exact_solution(gaussian, xx, 1.0, 1.0, tf, sigma=0.1)
        _, ddx = make_grid(n, 1.0)
        errs.append(l2_error(u, ue, ddx))
    curves[name] = (np.array(errs), col, np.log2(errs[-2] / errs[-1]))

# Lax-Wendroff and BD-2 land on top of each other here, and that is a RESULT,
# not a plotting accident: at C = 1/2 their leading dispersive coefficients,
# (1-C^2) and (1-C)(2-C), are equal. Drawn as one solid line with the other
# dashed over it so both remain visible.
for name in ("FTBS", "Lax-Friedrichs"):
    e, col, p = curves[name]
    ax[1].loglog(grids, e, "o-", color=col, lw=2.0, ms=6,
                 label=f"{name}   (p = {p:.2f})")

e, col, p = curves["Lax-Wendroff"]
ax[1].loglog(grids, e, "o-", color=col, lw=2.6, ms=7,
             label=f"Lax–Wendroff   (p = {p:.2f})", zorder=3)
e2, col2, p2 = curves["BD-2"]
ax[1].loglog(grids, e2, ls=(0, (3, 3)), color=col2, lw=2.0, marker="o",
             ms=4, mfc=SURFACE, label=f"BD-2   (p = {p2:.2f})", zorder=4)

assert np.allclose(e, e2), "LW and BD-2 should coincide exactly at C = 1/2"
ax[1].annotate("Lax–Wendroff and BD-2 coincide\nexactly — and only at C = 1/2",
               xy=(780, 1.15e-4), xytext=(420, 1.1e-3),
               color=INK, fontsize=9, ha="center",
               arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.1,
                               connectionstyle="arc3,rad=-0.15"))

# reference slopes, drawn faintly behind the data
g = np.array(grids, dtype=float)
ax[1].loglog(g, 6e-2 * (g / g[0]) ** -1.0, color=MUTED, lw=1.0,
             ls=(0, (5, 4)), zorder=1)
ax[1].loglog(g, 3e-3 * (g / g[0]) ** -2.0, color=MUTED, lw=1.0,
             ls=(0, (5, 4)), zorder=1)
ax[1].text(1700, 4.0e-3, "slope −1", color=MUTED, fontsize=9)
ax[1].text(1700, 1.1e-5, "slope −2", color=MUTED, fontsize=9)

ax[1].set_title("Grid convergence  (Gaussian, C = 0.5, one revolution)",
                color=INK, fontsize=11.5, pad=8)
ax[1].set_xlabel("N  (grid points)"); ax[1].set_ylabel("L₂ error")
style(ax[1])
lg = ax[1].legend(frameon=False, fontsize=9, loc="lower left")
for t in lg.get_texts():
    t.set_color(MUTED)

fig.suptitle("Module 01 — the two ways a scheme fails, and how fast each one converges",
             color=INK, fontsize=12.5, y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(OUT, "schemes_overview.png"), dpi=170,
            facecolor=SURFACE, bbox_inches="tight")
print("figures/schemes_overview.png")
