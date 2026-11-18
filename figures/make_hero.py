"""Build the repository's headline figures from the saved cavity runs.

Design rules followed here:
  - vorticity is a SIGNED field  -> diverging colormap, two hues, neutral
    midpoint. Never a rainbow.
  - Re = 100 / 400 / 1000 is an ORDERED quantity -> cool-to-warm progression,
    not arbitrary categorical colours.
  - palette validated for colour-vision deficiency before use.
  - explicit light surface, because GitHub renders READMEs in dark mode too.
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

_HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    _HERE, "..", "06-lid-driven-cavity", "runs")
OUT = sys.argv[2] if len(sys.argv) > 2 else _HERE
os.makedirs(OUT, exist_ok=True)

# --- design tokens -------------------------------------------------------
SURFACE = "#FFFFFF"
INK = "#1A1A1A"
INK_MUTED = "#5A5A5A"
GRID = "#D8D8D8"

# validated ordered palette: all six checks pass, no warnings
RE_COLOR = {100: "#2E7DB8", 400: "#B8841A", 1000: "#A8322A"}

# diverging map for vorticity: blue -> neutral gray -> red
VORT = LinearSegmentedColormap.from_list(
    "vort", ["#1B4F79", "#2E7DB8", "#9FC4DE", "#EFEFEE", "#E0A79F", "#A8322A", "#6E1F18"]
)


def style(ax, grid=True):
    ax.set_facecolor(SURFACE)
    for s in ax.spines.values():
        s.set_color(GRID)
        s.set_linewidth(0.8)
    ax.tick_params(colors=INK_MUTED, labelsize=9, length=3, width=0.8)
    ax.xaxis.label.set_color(INK_MUTED)
    ax.yaxis.label.set_color(INK_MUTED)
    if grid:
        ax.grid(True, color=GRID, linewidth=0.6, alpha=0.7)
        ax.set_axisbelow(True)


def load(Re):
    return dict(np.load(os.path.join(RUNS, f"Re{Re}_N101.npz")))


runs = {Re: load(Re) for Re in (100, 400, 1000)}

# =========================================================================
# 1. HERO — vorticity field with streamlines, three Reynolds numbers
# =========================================================================
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.9), facecolor=SURFACE)
CLIP = 28.0  # symmetric clip: the lid singularity otherwise swamps the interior

# streamfunction levels chosen to expose the corner vortices as well as the core
lv_primary = np.array([-0.1, -0.09, -0.07, -0.05, -0.03, -0.015, -0.005, -1e-3, -1e-4])
lv_corner = np.array([1e-8, 1e-6, 1e-5, 5e-5, 1.5e-4, 3e-4, 6e-4, 1.2e-3])

for ax, Re in zip(axes, (100, 400, 1000)):
    r = runs[Re]
    X, Y, psi, om = r["X"], r["Y"], r["psi"], r["omega"]

    cf = ax.contourf(X, Y, np.clip(om, -CLIP, CLIP), levels=42,
                     cmap=VORT, norm=TwoSlopeNorm(vmin=-CLIP, vcenter=0.0, vmax=CLIP))
    try:                       # matplotlib >= 3.8 made ContourSet an Artist
        cf.set_rasterized(True)
    except AttributeError:     # older matplotlib exposes .collections
        for c in cf.collections:
            c.set_rasterized(True)

    ax.contour(X, Y, psi, levels=lv_primary, colors="#15202B",
               linewidths=0.75, alpha=0.85)
    ax.contour(X, Y, psi, levels=lv_corner, colors="#15202B",
               linewidths=0.6, alpha=0.9)

    ax.set_title(f"Re = {Re}", color=INK, fontsize=12.5, pad=8)
    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_yticks([0, 0.5, 1.0])
    style(ax, grid=False)
    for s in ax.spines.values():
        s.set_color("#9A9A9A")

axes[0].set_ylabel("y", fontsize=10)
for ax in axes:
    ax.set_xlabel("x", fontsize=10)

cb = fig.colorbar(cf, ax=axes, fraction=0.019, pad=0.015, extend="both")
cb.set_label("vorticity  ω", color=INK_MUTED, fontsize=10)
cb.ax.tick_params(colors=INK_MUTED, labelsize=9, length=3, width=0.8)
cb.outline.set_edgecolor(GRID)

fig.suptitle(
    "Lid-driven cavity, 101×101 — vorticity field with streamfunction contours",
    color=INK, fontsize=13.5, y=0.985)
fig.text(0.5, 0.015,
         "The lid drives the flow from above. As Re rises the primary vortex moves "
         "toward the cavity centre and the corner vortices grow.",
         ha="center", color=INK_MUTED, fontsize=9.5)
fig.savefig(os.path.join(OUT, "cavity_hero.png"), dpi=170,
            facecolor=SURFACE, bbox_inches="tight")
plt.close(fig)
print("cavity_hero.png")

# =========================================================================
# 2. VALIDATION — centreline velocities against Ghia
# =========================================================================
import csv
import io

REF = os.path.join(os.path.dirname(RUNS), "reference")


def load_ghia(table):
    fn = {"u": "ghia1982_table1_u.csv", "v": "ghia1982_table2_v.csv"}[table]
    rows = [l for l in open(os.path.join(REF, fn)) if not l.startswith("#")]
    rd = list(csv.DictReader(io.StringIO("".join(rows))))
    key = "y" if table == "u" else "x"
    coords = np.array([float(r[key]) for r in rd])
    data = {int(k[2:]): np.array([float(r[k]) for r in rd])
            for k in rd[0] if k.startswith("Re")}
    return coords, data


yc, gu = load_ghia("u")
xc, gv = load_ghia("v")

fig, ax = plt.subplots(1, 2, figsize=(11.6, 4.8), facecolor=SURFACE)

for Re in (100, 400, 1000):
    r = runs[Re]
    N = int(r["N"])
    mid = (N - 1) // 2
    c = RE_COLOR[Re]
    ax[0].plot(r["u"][:, mid], r["Y"][:, mid], "-", color=c, lw=2.0, zorder=3)
    ax[0].plot(gu[Re], yc, "o", color=c, mfc="none", ms=6.5, mew=1.5, zorder=4)
    ax[1].plot(r["X"][mid, :], r["v"][mid, :], "-", color=c, lw=2.0, zorder=3,
               label=f"Re = {Re}")
    ax[1].plot(xc, gv[Re], "o", color=c, mfc="none", ms=6.5, mew=1.5, zorder=4)

ax[0].axvline(0, color=GRID, lw=0.8, zorder=1)
ax[1].axhline(0, color=GRID, lw=0.8, zorder=1)
ax[0].set_xlabel("u", fontsize=10)
ax[0].set_ylabel("y", fontsize=10)
ax[0].set_title("u on the vertical centreline  (x = 0.5)", color=INK, fontsize=11.5, pad=7)
ax[1].set_xlabel("x", fontsize=10)
ax[1].set_ylabel("v", fontsize=10)
ax[1].set_title("v on the horizontal centreline  (y = 0.5)", color=INK, fontsize=11.5, pad=7)

# direct labels on the left panel so identity is never colour-alone.
# placed in the empty upper-left quadrant: the curves all live to the right of
# u = -0.4 above y = 0.7, so nothing collides.
for Re, y in ((1000, 0.95), (400, 0.885), (100, 0.82)):
    ax[0].text(-0.385, y, f"Re = {Re}", color=RE_COLOR[Re],
               fontsize=10, fontweight="bold", va="center", ha="left")

for a in ax:
    style(a)
leg = ax[1].legend(frameon=False, fontsize=9.5, loc="upper right")
for t in leg.get_texts():
    t.set_color(INK_MUTED)

# A second legend explaining the SYMBOLS, not the colours. Without it the only
# statement that circles are the benchmark lives in the figure title, which is
# easy to miss — the first person to read this figure asked where Ghia's data
# was.
from matplotlib.lines import Line2D  # noqa: E402

def symbol_handles():
    return [
        Line2D([0], [0], color="#444444", lw=2.0,
               label="present solver"),
        Line2D([0], [0], marker="o", color="none", mec="#444444", mew=1.6,
               ms=7, ls="none", label="Ghia et al. (1982) benchmark"),
    ]


# One on EACH panel. A reader looking at only the right-hand plot must not have
# to infer what a circle is from the other half of the figure.
for a, loc in ((ax[0], "lower right"), (ax[1], "lower left")):
    lg = a.legend(handles=symbol_handles(), frameon=True, fontsize=9,
                  loc=loc, handletextpad=0.8, borderpad=0.7)
    lg.get_frame().set_edgecolor(GRID)
    lg.get_frame().set_facecolor(SURFACE)
    for t in lg.get_texts():
        t.set_color(INK)

ax[1].add_artist(leg)   # keep the Re colour legend alongside the symbol one

fig.suptitle("Centreline velocities: present solver vs. the Ghia, Ghia & Shin (1982) benchmark",
             color=INK, fontsize=12.5, y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.945])
fig.savefig(os.path.join(OUT, "validation.png"), dpi=170,
            facecolor=SURFACE, bbox_inches="tight")
plt.close(fig)
print("validation.png")
