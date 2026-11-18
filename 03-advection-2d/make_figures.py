"""Figures for module 03, generated from src/.

    python3 make_figures.py

Both panels show something that has no 1D analogue and that an axis-aligned
test cannot detect.
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from src.initial_conditions_2d import exact_solution_2d, make_grid_2d, sine_2d
from src.schemes_2d import lax_wendroff_2d, lax_wendroff_2d_no_cross, upwind_2d
from src.solver_2d import l2_error_2d, solve_2d

np.seterr(all="ignore")

OUT = os.path.join(_HERE, "figures")
os.makedirs(OUT, exist_ok=True)

SURFACE, INK, MUTED, GRID = "#FFFFFF", "#1A1A1A", "#5A5A5A", "#D8D8D8"
BLUE, GOLD, RED = "#2E7DB8", "#B8841A", "#A8322A"


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ax.spines.values():
        s.set_color(GRID); s.set_linewidth(0.8)
    ax.tick_params(colors=MUTED, labelsize=9, length=3, width=0.8)
    ax.xaxis.label.set_color(MUTED); ax.yaxis.label.set_color(MUTED)
    ax.grid(True, color=GRID, linewidth=0.6, alpha=0.7); ax.set_axisbelow(True)


fig, ax = plt.subplots(1, 2, figsize=(12.8, 4.7), facecolor=SURFACE)

# ---- (a) the stability boundary is Cx + Cy <= 1 -------------------------
# Tested on a checkerboard, which is exactly the theta_x = theta_y = pi mode
# the Von Neumann analysis identifies as the worst case. A smooth initial
# condition contains almost no energy there and stays bounded well past the
# true limit, which is how this instability stays hidden.
N, NSTEP = 64, 200
idx = np.arange(N)
u0 = ((idx[None, :] + idx[:, None]) % 2).astype(float) * 2.0 - 1.0

Cs = np.linspace(0.40, 0.60, 21)
amp = []
for C in Cs:
    u = u0.copy()
    for _ in range(NSTEP):
        u = upwind_2d(u, C, C)
    amp.append(np.abs(u).max())
amp = np.array(amp)

ax[0].axvspan(1.0, 1.25, color=RED, alpha=0.055, zorder=1)
ax[0].semilogy(2 * Cs, np.clip(amp, 1e-25, None), "o-", color=BLUE,
               lw=2.0, ms=5, zorder=3)
ax[0].axvline(1.0, color="#444444", lw=1.3, ls=(0, (4, 3)), zorder=2)
ax[0].axhline(1.0, color=MUTED, lw=1.0, zorder=2)

ax[0].text(1.012, 1e-12, "unstable\n$C_x+C_y>1$", color=RED, fontsize=10,
           fontweight="bold", va="center")
ax[0].text(0.805, 1e-12, "stable", color=BLUE, fontsize=10,
           fontweight="bold", va="center")
ax[0].text(0.805, 2.2, "neutral  |u| = 1", color=MUTED, fontsize=9)

ax[0].set_xlabel(r"$C_x + C_y$")
ax[0].set_ylabel(f"max |u| after {NSTEP} steps")
ax[0].set_title(f"2D stability limit, checkerboard mode  ({N}×{N})",
                color=INK, fontsize=11.5, pad=8)
ax[0].set_xlim(0.79, 1.21)
ax[0].set_ylim(1e-25, 1e18)
style(ax[0])

# ---- (b) Lax-Wendroff needs its cross-derivative term -------------------
grids = (40, 80, 160, 320)
T, L = 0.5, 1.0


def order(step, Cx, Cy, cx, cy):
    errs = []
    for n in grids:
        X, Y, u, tf = solve_2d(step, sine_2d, n, Cx, Cy, T, cx, cy, L)
        ue = exact_solution_2d(sine_2d, X, Y, L, cx, cy, tf)
        _, _, dx = make_grid_2d(n, L)
        errs.append(l2_error_2d(u, ue, dx))
    return np.array(errs), np.log2(errs[-2] / errs[-1])


cases = [
    ("with cross term, diagonal flow", lax_wendroff_2d, 0.25, 0.25, 1.0, 1.0, BLUE, "-"),
    ("no cross term, flow along x", lax_wendroff_2d_no_cross, 0.25, 0.0, 1.0, 0.0, GOLD, "-"),
    ("no cross term, diagonal flow", lax_wendroff_2d_no_cross, 0.25, 0.25, 1.0, 1.0, RED, "-"),
]
for lab, step, Cx, Cy, cx, cy, col, ls in cases:
    e, p = order(step, Cx, Cy, cx, cy)
    ax[1].loglog(grids, e, "o" + ls, color=col, lw=2.0, ms=6,
                 label=f"{lab}   (p = {p:.2f})")

g = np.array(grids, dtype=float)
ax[1].loglog(g, 1.3e-1 * (g / g[0]) ** -1.0, color=MUTED, lw=1.0, ls=(0, (5, 4)), zorder=1)
ax[1].loglog(g, 1.3e-2 * (g / g[0]) ** -2.0, color=MUTED, lw=1.0, ls=(0, (5, 4)), zorder=1)
ax[1].text(335, 2.1e-2, "slope −1", color=MUTED, fontsize=9)
ax[1].text(335, 1.4e-4, "slope −2", color=MUTED, fontsize=9)

ax[1].set_xlabel("N  (points per direction)")
ax[1].set_ylabel("L₂ error")
ax[1].set_title("Dropping the $u_{xy}$ term costs a full order — on diagonals only",
                color=INK, fontsize=11.5, pad=8)
style(ax[1])
lg = ax[1].legend(frameon=False, fontsize=8.5, loc="lower left")
for t in lg.get_texts():
    t.set_color(MUTED)

fig.suptitle("Module 03 — two 2D results that an axis-aligned test cannot see",
             color=INK, fontsize=12.5, y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(OUT, "two_dimensional_results.png"), dpi=170,
            facecolor=SURFACE, bbox_inches="tight")
print("figures/two_dimensional_results.png")
