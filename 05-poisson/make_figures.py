"""Figures for module 05, generated from src/poisson.py.

    python3 make_figures.py

Panel (a) takes a minute: Gauss-Seidel and SOR are written with explicit loops
because their sequential dependency is the whole point of the method.
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from matplotlib.ticker import FixedFormatter, FixedLocator, NullFormatter, NullLocator

from src.poisson import gauss_seidel, jacobi, optimal_omega, residual, sor

OUT = os.path.join(_HERE, "figures")
os.makedirs(OUT, exist_ok=True)

SURFACE, INK, MUTED, GRID = "#FFFFFF", "#1A1A1A", "#5A5A5A", "#D8D8D8"
BLUE, GOLD, RED = "#2E7DB8", "#B8841A", "#A8322A"


def log_ticks(ax, vals):
    """Label a log axis at exactly these values.

    set_xticks alone is not enough: matplotlib's log formatter keeps drawing
    its own 2x10^1 / 3x10^1 labels underneath, and they collide.
    """
    ax.xaxis.set_major_locator(FixedLocator(vals))
    ax.xaxis.set_major_formatter(FixedFormatter([str(v) for v in vals]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ax.spines.values():
        s.set_color(GRID); s.set_linewidth(0.8)
    ax.tick_params(colors=MUTED, labelsize=9, length=3, width=0.8)
    ax.xaxis.label.set_color(MUTED); ax.yaxis.label.set_color(MUTED)
    ax.grid(True, color=GRID, linewidth=0.6, alpha=0.7); ax.set_axisbelow(True)


def manufactured(N, L=1.0):
    """Pick the answer first: s = sin(pi x) sin(pi y), so lap(s) = -2 pi^2 s."""
    dx = L / (N - 1)
    x = np.linspace(0, L, N)
    X, Y = np.meshgrid(x, x)
    exact = np.sin(np.pi * X) * np.sin(np.pi * Y)
    return np.zeros((N, N)), -2.0 * np.pi**2 * exact, exact, dx


fig, ax = plt.subplots(1, 2, figsize=(12.8, 4.7), facecolor=SURFACE)

# ---- (a) iteration count is a complexity class, not a constant ----------
grids = (11, 21, 41)
counts = {"Jacobi": [], "Gauss–Seidel": [], "SOR": []}
for N in grids:
    s0, f, _, dx = manufactured(N)
    counts["Jacobi"].append(jacobi(s0, f, dx, tol=1e-6)[1])
    counts["Gauss–Seidel"].append(gauss_seidel(s0, f, dx, tol=1e-6)[1])
    counts["SOR"].append(sor(s0, f, dx, tol=1e-6)[1])

for name, col in (("Jacobi", RED), ("Gauss–Seidel", GOLD), ("SOR", BLUE)):
    c = np.array(counts[name], dtype=float)
    p = np.log2(c[-1] / c[-2])
    ax[0].loglog(grids, c, "o-", color=col, lw=2.2, ms=7,
                 label=f"{name}   (exponent {p:.2f})")

g = np.array(grids, dtype=float)
ax[0].loglog(g, 335 * (g / g[0]) ** 2.0, color=MUTED, lw=1.0, ls=(0, (5, 4)), zorder=1)
ax[0].loglog(g, 36 * (g / g[0]) ** 1.0, color=MUTED, lw=1.0, ls=(0, (5, 4)), zorder=1)
ax[0].text(43, 6200, "$O(N^2)$", color=MUTED, fontsize=9.5)
ax[0].text(43, 120, "$O(N)$", color=MUTED, fontsize=9.5)

ax[0].set_xlabel("N  (points per direction)")
ax[0].set_ylabel("sweeps to reach residual < 1e-6")
ax[0].set_title("Jacobi is $O(N^2)$; SOR with optimal ω is $O(N)$",
                color=INK, fontsize=11.5, pad=8)
style(ax[0])
log_ticks(ax[0], grids)
lg = ax[0].legend(frameon=False, fontsize=9, loc="upper left")
for t in lg.get_texts():
    t.set_color(MUTED)

# ---- (b) the stopping test that silently lies --------------------------
# Stopping on max|s_new - s_old| is not stopping on the residual. For Jacobi
# the two differ by exactly 4/dx^2 -- a factor that GROWS as the grid refines,
# so the test is worst precisely when you refine to get a better answer.
Ns = (11, 21, 41, 81)
chg, res, pred = [], [], []
for N in Ns:
    dx = 1.0 / (N - 1)
    f = np.zeros((N, N))
    s = np.zeros((N, N))
    s[-1, :] = 1.0                      # a lid, as in Laplace
    for _ in range(200000):
        so = s.copy()
        s[1:-1, 1:-1] = 0.25 * (so[1:-1, 2:] + so[1:-1, :-2]
                                + so[2:, 1:-1] + so[:-2, 1:-1])
        c = np.abs(s - so).max()
        if c < 1e-5:
            break
    chg.append(c)
    res.append(np.abs(residual(s, f, dx)).max())
    pred.append(4.0 / dx**2)

ax[1].loglog(Ns, res, "o-", color=RED, lw=2.2, ms=7,
             label="true residual  |∇²s − f|")
ax[1].loglog(Ns, chg, "o-", color=BLUE, lw=2.2, ms=7,
             label="what was actually tested:  max|sₖ₊₁ − sₖ|")
ax[1].axhline(1e-5, color=MUTED, lw=1.0, ls=":", zorder=1)
ax[1].text(11.4, 1.25e-5, "declared tolerance 1e-5", color=MUTED, fontsize=8.5)

for N, c, r in zip(Ns, chg, res):
    ax[1].annotate("", xy=(N, r), xytext=(N, c),
                   arrowprops=dict(arrowstyle="<->", color=MUTED, lw=1.0))
ax[1].text(44, 4e-3, f"gap = 4/Δx²\n{pred[-1]:,.0f}× at N = 81".replace(",", " "),
           color=INK, fontsize=9.5, fontweight="bold")

ax[1].set_xlabel("N  (points per direction)")
ax[1].set_ylabel("magnitude")
ax[1].set_title("Converged to 1e-5 — by one measure only",
                color=INK, fontsize=11.5, pad=8)
style(ax[1])
log_ticks(ax[1], Ns)
lg = ax[1].legend(frameon=False, fontsize=8.5, loc="upper left")
for t in lg.get_texts():
    t.set_color(MUTED)

fig.suptitle("Module 05 — how long an iterative solve takes, and when it is actually finished",
             color=INK, fontsize=12.5, y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(OUT, "iterative_solvers.png"), dpi=170,
            facecolor=SURFACE, bbox_inches="tight")
print("figures/iterative_solvers.png")
print("  omega_opt:", [f"{optimal_omega(n):.4f}" for n in grids])
