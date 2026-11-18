"""Figures for module 04, generated from src/burgers.py.

    python3 make_figures.py

Time is an ordered quantity, so the profiles use a single-hue ramp running
light to dark, not arbitrary categorical colours.
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from src.burgers import solve_1d, solve_2d

np.seterr(all="ignore")

OUT = os.path.join(_HERE, "figures")
os.makedirs(OUT, exist_ok=True)

SURFACE, INK, MUTED, GRID = "#FFFFFF", "#1A1A1A", "#5A5A5A", "#D8D8D8"
RED = "#A8322A"
RAMP = LinearSegmentedColormap.from_list(
    "t", ["#A8CBE6", "#6AA6D4", "#3B82BE", "#1F5C92", "#0C3454"])
VORT = LinearSegmentedColormap.from_list(
    "d", ["#1B4F79", "#2E7DB8", "#9FC4DE", "#EFEFEE", "#E0A79F", "#A8322A", "#6E1F18"])


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ax.spines.values():
        s.set_color(GRID); s.set_linewidth(0.8)
    ax.tick_params(colors=MUTED, labelsize=9, length=3, width=0.8)
    ax.xaxis.label.set_color(MUTED); ax.yaxis.label.set_color(MUTED)
    ax.grid(True, color=GRID, linewidth=0.6, alpha=0.7); ax.set_axisbelow(True)


fig = plt.figure(figsize=(16.2, 4.7), facecolor=SURFACE)
gs = fig.add_gridspec(1, 3, width_ratios=[1.12, 1.0, 0.92], wspace=0.28)
ax0, ax1, ax2 = (fig.add_subplot(gs[i]) for i in range(3))

# ---- (a) a smooth wave steepening into a shock --------------------------
# N = 1600 to match tests/test_equivalence_burgers.py and the table in
# README.md. On a coarser grid the measured gradient departs from theory
# earlier, simply because the steepening profile is less well resolved.
N, L, nu = 1600, 1.0, 0.0
dx = L / N
x = np.linspace(0, L - dx, N)
u0 = np.sin(2 * np.pi * x)
dt = 0.2 * dx
slope0 = 2 * np.pi
t_break = 1.0 / slope0                       # = 1/(2 pi) ~ 0.1592

times = [0.0, 0.05, 0.10, 0.13, 0.150]
for k, t in enumerate(times):
    u = solve_1d(u0, dx, dt, nu, int(round(t / dt)))
    ax0.plot(x, u, color=RAMP(k / (len(times) - 1)), lw=2.0,
             label=f"t = {t:.3f}", zorder=3 + k)

ax0.set_xlabel("x"); ax0.set_ylabel("u")
ax0.set_title("Inviscid Burgers: a smooth wave steepens",
              color=INK, fontsize=11.5, pad=8)
ax0.set_xlim(0, 1)
style(ax0)
lg = ax0.legend(frameon=False, fontsize=8.5, loc="upper right")
for t_ in lg.get_texts():
    t_.set_color(MUTED)
ax0.annotate("gradient → ∞\nat $t_b = 1/2\\pi$", xy=(0.5, 0.0),
             xytext=(0.56, 0.62), color=INK, fontsize=9,
             arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.1))

# ---- (b) the gradient against its closed form ---------------------------
ts = np.linspace(0.0, 0.145, 16)
meas = []
for t in ts:
    u = solve_1d(u0, dx, dt, nu, int(round(t / dt)))
    meas.append(np.abs((np.roll(u, -1) - np.roll(u, 1)) / (2 * dx)).max())
meas = np.array(meas)

tt = np.linspace(0, 0.152, 400)
ax1.plot(tt, slope0 / (1 - tt * slope0), color="#111111", lw=1.6,
         ls=(0, (4, 3)), label=r"theory  $2\pi/(1-2\pi t)$", zorder=2)
ax1.plot(ts, meas, "o", color=RED, ms=7, mfc="none", mew=1.8,
         label="measured", zorder=3)
ax1.axvline(t_break, color=MUTED, lw=1.0, ls=":", zorder=1)
ax1.text(t_break - 0.004, 20, "$t_b = 1/2\\pi$", color=MUTED, fontsize=9,
         rotation=90, va="center", ha="right")

# Quote the agreement only over the range where the shock is still RESOLVED on
# the grid. Past t ~ 0.13 the profile is steeper than 800 points can represent,
# so the measured gradient saturates below theory -- a grid limit, not a solver
# error. Including those points would understate an agreement that is genuinely
# sub-1% where the comparison is meaningful.
theory = slope0 / (1 - ts * slope0)
err = np.abs(meas - theory) / theory
ok = ts <= 0.13

ax1.set_xlabel("t"); ax1.set_ylabel(r"max $|\partial u/\partial x|$")
ax1.set_title(f"Measured gradient vs. theory  (within {100*err[ok].max():.1f} % to t = 0.13)",
              color=INK, fontsize=11.5, pad=8)
ax1.annotate("beyond here the shock is\nsteeper than the grid resolves",
             xy=(0.137, 41), xytext=(0.035, 44), color=MUTED, fontsize=8.5,
             arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.0))
ax1.set_xlim(0, 0.158); ax1.set_ylim(0, 60)
style(ax1)
lg = ax1.legend(frameon=False, fontsize=9, loc="upper left")
for t_ in lg.get_texts():
    t_.set_color(MUTED)

# ---- (c) the 2D solver, from the original hand-written setup ------------
L2, n2, nu2 = 10.0, 101, 0.001
d2 = L2 / (n2 - 1)
xs = np.linspace(0, L2, n2)
X, Y = np.meshgrid(xs, xs)
U0 = np.sin(np.pi * X) * np.sin(np.pi * Y)
U, V = solve_2d(U0, U0.copy(), d2, d2, 0.01, nu2, 100)

lim = np.abs(U).max()
cf = ax2.contourf(X, Y, U, levels=38, cmap=VORT, vmin=-lim, vmax=lim)
try:
    cf.set_rasterized(True)
except AttributeError:
    for c in cf.collections:
        c.set_rasterized(True)
ax2.set_aspect("equal")
ax2.set_xlabel("x"); ax2.set_ylabel("y")
ax2.set_title("2D Burgers at t = 1 — each lobe steepens\nalong the diagonal it is advected",
              color=INK, fontsize=11, pad=8)
style(ax2)
cb = fig.colorbar(cf, ax=ax2, fraction=0.046, pad=0.03)
cb.set_label("u", color=MUTED, fontsize=10)
cb.ax.tick_params(colors=MUTED, labelsize=9, length=3, width=0.8)
cb.outline.set_edgecolor(GRID)

fig.suptitle("Module 04 — nonlinearity: the solution carries itself, and steepens until it breaks",
             color=INK, fontsize=12.5, y=1.0)
fig.savefig(os.path.join(OUT, "burgers_overview.png"), dpi=170,
            facecolor=SURFACE, bbox_inches="tight")
print("figures/burgers_overview.png")
