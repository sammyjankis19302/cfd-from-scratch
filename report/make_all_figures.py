"""
make_all_figures.py  --  every figure in the internship report, from Python only.

Put this file in the top folder of the cavity work (the folder that contains
"Ghia validation", "Grig study", compare_ghia.py, ...) and run

    python3 make_all_figures.py

It needs numpy and matplotlib (pip install numpy matplotlib).
The figures are written to  report_figures/  next to this file:

    report_figures/schematics/   8 drawings (no data): roadmap, stencils, characteristics,
                                 five-point stencil, cavity, ghost point, algorithm cycle,
                                 Ghia's configuration
    report_figures/journey/      the advection/Burgers/Poisson figures, recomputed with the
                                 parameters of my August-September scripts, and the
                                 201x201 overflow run (takes about a minute)
    report_figures/06_cavity/    every cavity figure, read from the saved
                                 solution_fields.npz and grid_independence_summary.csv

Copy the three folders into the report's figures/ folder (same names) to update the report.
The 11 September picture (figures/journey/sep11_cavity.png) is my original saved image
and is not regenerated here.

To run only some parts:  python3 make_all_figures.py schematics journey overflow cavity
"""
import os, sys, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch, Circle

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT_SCH = os.path.join(ROOT, "report_figures", "schematics")
OUT_J = os.path.join(ROOT, "report_figures", "journey") + os.sep
OUT_C = os.path.join(ROOT, "report_figures", "06_cavity") + os.sep
for d in (OUT_SCH, OUT_J, OUT_C):
    os.makedirs(d, exist_ok=True)


# =====================================================================================
# 1. SCHEMATICS (drawings only, no computed data)
# =====================================================================================
BLUE, RED, GREEN, GREY, INK = "#1f4e8c", "#a3312a", "#2e7d32", "#666666", "#111111"

def canvas(w, h, xlim, ylim):
    # the axes fill the whole figure, so one data unit is exactly w/(x1-x0) inches;
    # the height follows from the data range (h is ignored and kept only for readability)
    s = w / (xlim[1] - xlim[0])
    fig = plt.figure(figsize=(w, s * (ylim[1] - ylim[0])))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_aspect("equal"); ax.axis("off")
    return fig, ax


def box(ax, x, y, w, h, text, fc="white", ec=INK, fs=8.2, weight=None):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc=fc, ec=ec, lw=0.8, zorder=2))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, zorder=3, fontweight=weight, linespacing=1.35)


def arrow(ax, p, q, color=INK, lw=1.0, ls="-", rad=0.0, ms=10):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=ms, color=color, lw=lw, ls=ls,
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0, zorder=4))


def arc_arrow(ax, cx, cy, rx, ry, t0, t1, color, lw=1.1):
    t = np.linspace(np.radians(t0), np.radians(t1), 120)
    x, y = cx + rx * np.cos(t), cy + ry * np.sin(t)
    ax.plot(x[:-3], y[:-3], color=color, lw=lw, zorder=3)
    arrow(ax, (x[-4], y[-4]), (x[-1], y[-1]), color=color, lw=lw, ms=9)


# ------------------------------------------------------------------ 1. roadmap
def roadmap():
    fig, ax = canvas(6.4, 3.3, (-2.1, 16.6), (-6.25, 1.45))
    W, H = 3.9, 1.25
    nodes = {
        "a": (0, 0, "Day one\nNavier-Stokes derived,\nadvection solved by hand", "1-9 Aug"),
        "b": (5.9, 0, "1  Advection\nFTFS, FTBS, FTCS;\nFourier analysis", "Stage 1, 12-31 Aug"),
        "c": (11.8, 0, "2  Diffusion\n$u_t=\\nu\\,u_{xx}$\nlimit $D\\leq 1/2$", "Stage 2"),
        "d": (11.8, -2.65, "3  Nonlinearity\nupwinding on the\nsign of $u$", "Stage 3, 31 Aug"),
        "e": (5.9, -2.65, "4  Advection-diffusion\ntwo viscosities:\n$\\nu$ and $\\nu_{\\mathrm{num}}$", "Stage 4"),
        "f": (0, -2.65, "5  Burgers\nnonlinear adv.-diff.;\nthe Reynolds number", "Stage 5, 1-7 Sep"),
        "g": (0, -5.3, "6  Jacobi\nLaplace and Poisson,\n$\\nabla^2\\psi=-\\omega$", "Stage 6, 2-8 Sep"),
        "h": (5.9, -5.3, "7  Lid-driven cavity\nvorticity and\nstreamfunction", "Stage 7, 9-23 Sep"),
        "i": (11.8, -5.3, "Validation\n20 runs against Ghia,\n$Re \\leq 10^4$", "22-24 Sep"),
    }
    for k, (x, y, t, d) in nodes.items():
        box(ax, x, y, W, H, t, fc="#e8eef8" if k in ("h", "i") else "white", fs=7.4)
        ax.text(x - W / 2, y + H / 2 + 0.08, d, fontsize=7.0, color=BLUE, fontweight="bold", va="bottom")
    def h_arrow(k1, k2, label):
        (x1, y1, *_), (x2, y2, *_) = nodes[k1], nodes[k2]
        s = 1 if x2 > x1 else -1
        arrow(ax, (x1 + s * W / 2, y1), (x2 - s * W / 2, y2), color="#444444")
        ax.text((x1 + x2) / 2, y1 - 0.12, label, fontsize=6.4, color=GREY, ha="center", va="top", linespacing=1.2)
    def v_arrow(k1, k2, label):
        (x1, y1, *_), (x2, y2, *_) = nodes[k1], nodes[k2]
        arrow(ax, (x1, y1 - H / 2), (x2, y2 + H / 2 + 0.33), color="#444444")
        ax.text(x1 + 0.12, (y1 + y2) / 2 + 0.1, label, fontsize=6.4, color=GREY, ha="left", va="center", linespacing=1.2)
    h_arrow("a", "b", "exact\nsolution"); h_arrow("b", "c", "upwinding\ndiffuses")
    v_arrow("c", "d", "next: the flow\ncarries itself")
    h_arrow("d", "e", "the sign\ncheck"); h_arrow("e", "f", r"$\nu_{\mathrm{num}}/\nu$")
    v_arrow("f", "g", "the velocity\nis missing")
    h_arrow("g", "h", r"$\psi$ from $\omega$"); h_arrow("h", "i", "a tool")
    xs = 5.9 + W / 2 - 0.35
    ax.plot([xs, xs, 14.1, 14.1, 13.78], [H / 2, 1.28, 1.28, -5.3, -5.3], color=RED, lw=1.0, ls="--")
    arrow(ax, (14.1, -5.3), (13.77, -5.3), color=RED, lw=1.0)
    ax.text(14.25, -2.1, "the numerical\ndiffusion of\nStage 1 limits\nthe cavity", fontsize=6.4, color=RED, ha="left", va="center")
    fig.savefig(os.path.join(OUT_SCH, "roadmap.pdf")); plt.close(fig)


# ------------------------------------------------------------------ 2. stencils
def stencils_1d():
    fig, ax = canvas(5.8, 1.75, (-2.6, 9.0), (-2.75, 0.85))
    dx = 1.2
    rows = [(0, "FD", (0, 1), "takes information from downstream: unstable for $c>0$", RED),
            (-1, "BD", (-1, 0), "takes information from upstream: upwind, stable for $0\\leq C\\leq 1$", BLUE),
            (-2, "CD-2", (-1, 1), "symmetric: no dissipation; unstable with explicit Euler", INK)]
    for y, name, (a, b), note, col in rows:
        ax.plot([-1.7, 1.7], [y, y], color="#bbbbbb", lw=0.7, zorder=1)
        for k in (-1, 0, 1):
            ax.add_patch(Circle((k * dx, y), 0.11, fc="white", ec="#999999", lw=0.8, zorder=2))
        ax.add_patch(Circle((0, y), 0.14, fc="white", ec=INK, lw=1.1, zorder=3))
        for k in (a, b):
            ax.add_patch(Circle((k * dx, y), 0.11, fc=INK, ec=INK, zorder=4))
        ax.plot([a * dx, b * dx], [y, y], color=col, lw=1.6, zorder=3)
        ax.text(-2.1, y, name, ha="right", va="center", fontsize=8.5, fontweight="bold")
        ax.text(2.0, y, note, ha="left", va="center", fontsize=7.6, color=col)
    for k, t in [(-1, "$j-1$"), (0, "$j$"), (1, "$j+1$")]:
        ax.text(k * dx, -2.45, t, ha="center", va="top", fontsize=8)
    arrow(ax, (-1.6, 0.55), (1.6, 0.55), color=BLUE, lw=1.2)
    ax.text(1.75, 0.55, "flow, $c>0$", va="center", fontsize=7.6, color=BLUE)
    fig.savefig(os.path.join(OUT_SCH, "stencils_1d.pdf")); plt.close(fig)


# ------------------------------------------------------------------ 3. characteristics
def characteristics():
    fig, ax = plt.subplots(figsize=(4.8, 2.2))
    c, T = 1.6, 2.4
    for x0 in np.arange(-4.5, 6.1, 0.75):
        xs = np.array([x0, x0 + c * 3.1]); ts = np.array([0, 3.1])
        ax.plot(xs, ts, color="#cfcfcf", lw=0.7, zorder=1)
    x = np.linspace(0.3, 2.9, 200); f = 0.75 * np.exp(-3 * (x - 1.6) ** 2)
    ax.plot(x, f, color=BLUE, lw=1.6)
    ax.plot(x + c * T, T + f, color=RED, lw=1.6)
    ax.axhline(T, color="#999999", lw=0.6, ls="--", xmax=7.3 / 8.6)
    arrow(ax, (1.6, 0.78), (1.6 + c * T, T + 0.73), color="#333333", lw=1.1)
    ax.text(3.45, 1.55, r"$x-ct=\mathrm{const}$", rotation=np.degrees(np.arctan(1 / c)), transform_rotates_text=True,
            rotation_mode="anchor", fontsize=8, ha="center", va="center", bbox=dict(fc="white", ec="none", pad=0.5))
    ax.text(1.6, -0.3, r"$u(x,0)=f(x)$", color=BLUE, ha="center", va="top", fontsize=8)
    ax.text(7.35, T + 0.35, r"$u(x,t)=f(x-ct)$", color=RED, ha="left", fontsize=8)
    ax.set_xlim(0, 8.6); ax.set_ylim(-0.65, 3.3)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.spines["left"].set_bounds(0, 3.2); ax.spines["bottom"].set_bounds(0, 7.3)
    ax.set_xticks([]); ax.set_yticks([T]); ax.set_yticklabels(["$t$"])
    ax.set_xlabel("$x$", loc="right", labelpad=-4); ax.set_ylabel("time", loc="top")
    fig.savefig(os.path.join(OUT_SCH, "characteristics.pdf")); plt.close(fig)


# ------------------------------------------------------------------ 4. five-point stencil
def five_point():
    fig, ax = canvas(5.8, 2.05, (-2.4, 8.6), (-1.9, 1.9))
    h = 1.3
    for k in (-1, 0, 1):
        ax.plot([k * h, k * h], [-1.75, 1.75], color="#cccccc", lw=0.7)
        ax.plot([-1.75, 1.75], [k * h, k * h], color="#cccccc", lw=0.7)
    pts = {"E": (h, 0, r"$\psi_{i+1,j}$", (0.12, -0.3)), "W": (-h, 0, r"$\psi_{i-1,j}$", (-0.95, -0.3)),
           "N": (0, h, r"$\psi_{i,j+1}$", (0.12, 0.12)), "S": (0, -h, r"$\psi_{i,j-1}$", (0.12, -0.35))}
    for (x, y, t, (dxl, dyl)) in pts.values():
        ax.add_patch(Circle((x, y), 0.09, fc=INK, zorder=4))
        arrow(ax, (x * 0.88, y * 0.88), (x * 0.18, y * 0.18), color=BLUE, lw=1.1, ms=8)
        ax.text(x + dxl, y + dyl, t, fontsize=8)
    ax.add_patch(Circle((0, 0), 0.13, fc="#cfe0f5", ec=INK, lw=1.0, zorder=5))
    ax.text(-1.05, 0.28, r"$\psi_{i,j},\ \omega_{i,j}$", fontsize=8, bbox=dict(fc="white", ec="none", pad=0.6))
    ax.text(2.25, 0.55, "new value at the centre = average of the\nfour old neighbours $+\\ \\frac{1}{4}h^2\\omega_{i,j}$:",
            fontsize=8, va="center")
    ax.text(2.25, -0.45, r"$\psi^{k+1}_{i,j}=\frac{1}{4}\left(\psi^k_{i+1,j}+\psi^k_{i-1,j}+\psi^k_{i,j+1}+\psi^k_{i,j-1}+h^2\omega_{i,j}\right)$",
            fontsize=8.5, va="center")
    fig.savefig(os.path.join(OUT_SCH, "five_point.pdf")); plt.close(fig)


# ------------------------------------------------------------------ 5. cavity schematic
def cavity_schematic():
    fig, ax = canvas(5.9, 4.0, (-0.42, 2.05), (-0.34, 1.26))
    L, t = 1.0, 0.05
    # walls (hatched) and the moving lid
    for r in [(-t, -t, t, L + t), (-t, -t, L + 2 * t, t), (L, -t, t, L + t)]:
        ax.add_patch(Rectangle(r[:2], r[2], r[3], fc="white", ec="#8a8a8a", hatch="////", lw=0, zorder=1))
    ax.plot([0, 0, L, L], [L, 0, 0, L], color=INK, lw=1.2, zorder=2)
    ax.plot([-0.01, L + 0.01], [L + 0.015, L + 0.015], color=BLUE, lw=3.0, zorder=3)
    for x0 in (0.12, 0.35, 0.58, 0.81):
        arrow(ax, (x0, L + 0.075), (x0 + 0.13, L + 0.075), color=BLUE, lw=1.2, ms=8)
    ax.text(L / 2, L + 0.13, r"moving lid: $u=U=1,\ v=0$", ha="center", va="bottom", fontsize=8.5)
    # boundary conditions on the fixed walls
    ax.text(-0.075, L / 2, r"$u=v=0$", rotation=90, ha="right", va="center", fontsize=8)
    ax.text(L + 0.075, L / 2, r"$u=v=0$", rotation=90, ha="left", va="center", fontsize=8)
    ax.text(L / 2, -0.07, r"$u=v=0$", ha="center", va="top", fontsize=8)
    # primary vortex: clockwise, in the upper half of the cavity
    cx, cy = 0.55, 0.57
    arc_arrow(ax, cx, cy, 0.29, 0.27, 100, -230, BLUE, lw=1.3)
    arc_arrow(ax, cx, cy, 0.14, 0.13, 100, -230, BLUE, lw=0.8)
    ax.plot(cx, cy, "o", color=INK, ms=2.5, zorder=4)
    ax.text(0.04, 0.955, "primary vortex\n$\\psi<0$ (clockwise)", ha="left", va="top", fontsize=7.8,
            color=BLUE, linespacing=1.3)
    # corner eddies: counter-rotating (anticlockwise), in the two bottom corners
    arc_arrow(ax, 0.085, 0.085, 0.058, 0.058, -60, 250, RED, lw=0.9)
    arc_arrow(ax, L - 0.115, 0.115, 0.085, 0.085, -60, 250, RED, lw=0.9)
    ax.text(0.165, 0.085, "BL", color=RED, fontsize=8, ha="left", va="center")
    ax.text(L - 0.225, 0.115, "BR", color=RED, fontsize=8, ha="right", va="center")
    # dimension line and axes
    arrow(ax, (L / 2, -0.19), (0, -0.19), color="#555555", lw=0.6, ms=7)
    arrow(ax, (L / 2, -0.19), (L, -0.19), color="#555555", lw=0.6, ms=7)
    for x0 in (0, L):
        ax.plot([x0, x0], [-0.165, -0.215], color="#555555", lw=0.6)
    ax.text(L / 2, -0.215, "$L=1$", ha="center", va="top", fontsize=8)
    arrow(ax, (-0.36, -0.29), (-0.23, -0.29), lw=0.8, ms=7); arrow(ax, (-0.36, -0.29), (-0.36, -0.16), lw=0.8, ms=7)
    ax.text(-0.215, -0.29, "$x$", va="center", fontsize=8); ax.text(-0.36, -0.135, "$y$", ha="center", va="bottom", fontsize=8)
    # what is solved where
    ax.text(1.2, 0.93, "on all four walls:\n$\\psi=0$ (one streamline)\n$\\omega_w$ from $\\psi$ by\nthe ghost-point formula",
            fontsize=8, va="top", linespacing=1.4)
    ax.text(1.2, 0.54, "inside:\n$\\nabla^2\\psi=-\\omega$\n$\\omega_t+u\\,\\omega_x+v\\,\\omega_y=\\frac{1}{Re}\\nabla^2\\omega$",
            fontsize=8, va="top", linespacing=1.5)
    ax.text(1.2, 0.17, "corner eddies BL, BR:\n$\\psi>0$ (anticlockwise),\ngrow with $Re$", fontsize=8, va="top",
            color=RED, linespacing=1.4)
    fig.savefig(os.path.join(OUT_SCH, "cavity_schematic.pdf")); plt.close(fig)


# ------------------------------------------------------------------ 6. ghost point
def ghost_point():
    fig, ax = canvas(6.3, 2.35, (-3.0, 14.2), (-3.1, 1.75))
    ax.add_patch(Rectangle((-1.6, -2.9), 3.2, 2.9, fc="#eef3fb", ec="none", zorder=0))
    for y in (0, -1.2, -2.4): ax.plot([-1.6, 1.6], [y, y], color="#cccccc", lw=0.7)
    for x in (-1.2, 0, 1.2): ax.plot([x, x], [-2.9, 1.45], color="#cccccc", lw=0.7)
    ax.plot([-1.6, 1.6], [0, 0], color=BLUE, lw=3.0)
    arrow(ax, (-1.3, 0.3), (-0.3, 0.3), color=BLUE, lw=1.2)
    ax.text(-1.7, 0.05, "lid, $u=U_w$", ha="right", va="center", fontsize=8, color=BLUE)
    ax.add_patch(Circle((0, 1.2), 0.1, fc="white", ec=INK, lw=1.0, zorder=4))
    for y in (0, -1.2, -2.4): ax.add_patch(Circle((0, y), 0.1, fc=INK, zorder=4))
    ax.text(0.22, 1.2, r"$\psi_g$ (ghost point, outside the fluid)", va="center", fontsize=8)
    ax.text(0.22, -0.3, r"$\psi_0=0$ (on the wall)", va="center", fontsize=8)
    ax.text(0.22, -1.2, r"$\psi_1$", va="center", fontsize=8); ax.text(0.22, -2.4, r"$\psi_2$", va="center", fontsize=8)
    for y0, y1 in ((0, -1.2), (-1.2, -2.4)):
        arrow(ax, (-0.55, y0), (-0.55, y1), lw=0.6, ms=7); arrow(ax, (-0.55, y1), (-0.55, y0), lw=0.6, ms=7)
        ax.text(-0.65, (y0 + y1) / 2, "$h$", ha="right", va="center", fontsize=8)
    kw = dict(fontsize=8.6, va="center", ha="left", linespacing=2.1,
              bbox=dict(boxstyle="round,pad=0.45", fc="white", ec="#999999", lw=0.7))
    ax.text(6.2, 0.55, "1. choose $\\psi_g$ so that $u=\\partial\\psi/\\partial y=U_w$ at the wall:\n"
                       r"$\dfrac{2\psi_g+3\psi_0-6\psi_1+\psi_2}{6h}=U_w$" "\n"
                       r"$\Rightarrow\ \psi_g=3hU_w-1.5\psi_0+3\psi_1-0.5\psi_2$", **kw)
    ax.text(6.2, -1.95, "2. on the wall $\\psi_{xx}=0$, so $\\omega_w=-\\psi_{yy}$:\n"
                        r"$\omega_w=-\dfrac{\psi_g-2\psi_0+\psi_1}{h^2}=-\dfrac{8\psi_1-\psi_2}{2h^2}-\dfrac{3U_w}{h}$", **kw)
    fig.savefig(os.path.join(OUT_SCH, "ghost_point.pdf")); plt.close(fig)


# ------------------------------------------------------------------ 7. algorithm cycle
def algorithm_cycle():
    fig, ax = canvas(6.0, 3.3, (-3.5, 10.4), (-5.6, 0.75))
    W, H = 3.9, 1.05
    B = {"w": (0, 0, r"$\omega^n$ (interior)", "#e8eef8"),
         "p": (8, 0, "$\\psi$\n(warm start from previous $\\psi$)", "white"),
         "uv": (8, -2.4, r"$u,\ v$ (interior)", "white"),
         "ww": (8, -4.8, r"$\omega_w$ on the four walls", "white"),
         "wn": (0, -4.8, r"$\omega^{n+1}$", "#e8eef8"),
         "ck": (0, -2.4, "$\\delta_\\omega=\\max|\\omega^{n+1}-\\omega^n|$\n$<10^{-8}$ ?", "#e6f2e6")}
    for x, y, t, fc in B.values(): box(ax, x, y, W, H, t, fc=fc, fs=7.6)
    op = dict(fontsize=7, color="#555555", ha="center", linespacing=1.25)
    arrow(ax, (W / 2, 0), (8 - W / 2, 0))
    ax.text(4, 0.12, r"Jacobi on $\nabla^2\psi=-\omega$", va="bottom", **op)
    ax.text(4, -0.12, r"until $\max|\Delta\psi|<10^{-6}$", va="top", **op)
    arrow(ax, (8, -H / 2), (8, -2.4 + H / 2))
    ax.text(8.15, -1.2, "$u=\\psi_y$, $v=-\\psi_x$\ncentral differences", va="center", fontsize=7, color="#555555", ha="left")
    arrow(ax, (8, -2.4 - H / 2), (8, -4.8 + H / 2))
    ax.text(8.15, -3.6, "wall formula\n(ghost point)", va="center", fontsize=7, color="#555555", ha="left")
    arrow(ax, (8 - W / 2, -4.8), (W / 2, -4.8))
    ax.text(4, -4.68, "upwind advection (sign of $u$, $v$)", va="bottom", **op)
    ax.text(4, -4.92, "CD-2 diffusion, explicit Euler", va="top", **op)
    arrow(ax, (0, -4.8 + H / 2), (0, -2.4 - H / 2))
    arrow(ax, (0, -2.4 + H / 2), (0, -H / 2))
    ax.text(-0.15, -1.2, "no: $n\\to n+1$", ha="right", va="center", fontsize=7, color="#555555")
    arrow(ax, (-W / 2, -2.4), (-W / 2 - 0.9, -2.4))
    ax.text(-W / 2 - 1.0, -2.4, "yes:\nsteady\nstate", ha="right", va="center", fontsize=7.5)
    ax.text(4, -2.4, r"one time step $\Delta t$", ha="center", va="center", fontsize=8, color="#777777")
    fig.savefig(os.path.join(OUT_SCH, "algorithm_cycle.pdf")); plt.close(fig)

# ------------------------------------------------------------------ 8. benchmark configuration (after Ghia et al. 1982, Fig. 1)
def ghia_configuration():
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch
    fig, ax = canvas(5.9, 4.4, (-0.55, 2.25), (-0.42, 1.34))
    L, t = 1.0, 0.045
    for r in [(-t, -t, t, L + t), (-t, -t, L + 2 * t, t), (L, -t, t, L + t)]:
        ax.add_patch(Rectangle(r[:2], r[2], r[3], fc="white", ec="#8a8a8a", hatch="////", lw=0, zorder=1))
    ax.plot([0, 0, L, L], [L, 0, 0, L], color=INK, lw=1.2, zorder=2)
    ax.plot([-0.01, L + 0.01], [L + 0.012, L + 0.012], color=BLUE, lw=3.0, zorder=3)
    for x0 in (0.2, 0.62):
        arrow(ax, (x0, L + 0.07), (x0 + 0.2, L + 0.07), color=BLUE, lw=1.2, ms=8)
    ax.text(L / 2, L + 0.12, r"$u=1,\ v=0$", ha="center", va="bottom", fontsize=8.5)
    ax.text(-0.07, 0.465, r"$u=v=0$", rotation=90, ha="right", va="center", fontsize=8)
    ax.text(L + 0.07, 0.72, r"$u=v=0$", rotation=90, ha="left", va="center", fontsize=8)
    ax.text(0.5, -0.06, r"$u=v=0$", ha="center", va="top", fontsize=8,
            bbox=dict(fc="white", ec="none", pad=0.5), zorder=5)
    for s, (x, y, ha, va) in {"A": (-0.03, L + 0.03, "right", "bottom"), "B": (-0.03, -0.03, "right", "top"),
                              "C": (L + 0.03, -0.03, "left", "top"), "D": (L + 0.03, L + 0.03, "left", "bottom")}.items():
        ax.text(x, y, s, ha=ha, va=va, fontsize=8.5, fontweight="bold")
    # centrelines, where Ghia tabulates u(y) and v(x)
    ax.plot([0.5, 0.5], [-0.02, L + 0.02], color=GREEN, lw=0.9, ls=(0, (6, 2, 1, 2)), zorder=2)
    ax.plot([-0.02, L + 0.02], [0.5, 0.5], color=GREEN, lw=0.9, ls=(0, (6, 2, 1, 2)), zorder=2)
    ax.text(0.52, 0.93, r"$u(y)$ at $x=0.5$", color=GREEN, fontsize=7.5, va="center")
    ax.text(0.04, 0.47, r"$v(x)$ at $y=0.5$", color=GREEN, fontsize=7.5, va="top")
    # primary vortex
    arc_arrow(ax, 0.53, 0.58, 0.2, 0.19, 100, -230, BLUE, lw=1.2)
    ax.plot(0.53, 0.58, "+", color=BLUE, ms=7, mew=1.2)
    ax.text(0.55, 0.585, r"  $\psi_{\min}$", color=BLUE, fontsize=8, va="bottom")

    def sep(p0, c, p1, fill, corner):
        ax.add_patch(PathPatch(Path([p0, c, p1, corner, p0],
                                    [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.LINETO, Path.CLOSEPOLY]),
                               fc=fill, ec="none", zorder=1.5))
        ax.add_patch(PathPatch(Path([p0, c, p1], [Path.MOVETO, Path.CURVE3, Path.CURVE3]),
                               fc="none", ec=RED, lw=1.0, zorder=2))
    pink = "#f6dcda"
    sep((0, 0.30), (0.26, 0.26), (0.30, 0), pink, (0, 0))                 # BL1
    sep((0, 0.07), (0.06, 0.06), (0.07, 0), "#ecb9b5", (0, 0))            # BL2
    sep((L - 0.42, 0), (L - 0.36, 0.40), (L, 0.46), pink, (L, 0))         # BR1
    sep((L - 0.13, 0), (L - 0.11, 0.11), (L, 0.13), "#ecb9b5", (L, 0))    # BR2
    sep((0, 0.93), (0.24, 0.86), (0, 0.62), pink, (0, 0.62))                  # TL1
    for s, (x, y) in {"BL1": (0.11, 0.12), "BR1": (L - 0.2, 0.2), "TL1": (0.025, 0.78)}.items():
        ax.text(x, y, s, color=RED, fontsize=7.5, ha="center" if s != "TL1" else "left", va="center")
    ax.annotate("BL2", xy=(0.03, 0.03), xytext=(0.16, -0.19), fontsize=7.5, color=RED, ha="center",
                arrowprops=dict(arrowstyle="-", color=RED, lw=0.6))
    ax.annotate("BR2", xy=(L - 0.04, 0.04), xytext=(L - 0.12, -0.19), fontsize=7.5, color=RED, ha="center",
                arrowprops=dict(arrowstyle="-", color=RED, lw=0.6))

    def dim(p, q, label, off, ha="center", va="center"):
        arrow(ax, p, q, color="#555555", lw=0.6, ms=6); arrow(ax, q, p, color="#555555", lw=0.6, ms=6)
        ax.text((p[0] + q[0]) / 2 + off[0], (p[1] + q[1]) / 2 + off[1], label, fontsize=7.5, ha=ha, va=va)
    dim((-0.2, 0), (-0.2, 0.30), r"$V_{BL1}$", (-0.03, 0), ha="right")
    dim((0, -0.3), (0.30, -0.3), r"$H_{BL1}$", (0, -0.03), va="top")
    dim((L + 0.2, 0), (L + 0.2, 0.46), r"$V_{BR1}$", (0.03, 0), ha="left")
    dim((L - 0.42, -0.3), (L, -0.3), r"$H_{BR1}$", (0, -0.03), va="top")
    dim((-0.2, 0.62), (-0.2, 0.93), r"$V_{TL1}$", (-0.03, 0), ha="right")
    for yy in (0.62, 0.93):
        ax.plot([-0.23, -0.02], [yy, yy], color="#999999", lw=0.4)
    for xx, yy in ((0.30, 0), (L - 0.42, 0), (0, 0)):
        ax.plot([xx, xx], [-0.33, -0.02], color="#999999", lw=0.4)
    ax.plot([L, L], [-0.33, -0.05], color="#999999", lw=0.4)
    ax.plot([-0.23, -0.05], [0.30, 0.30], color="#999999", lw=0.4)
    ax.plot([L + 0.05, L + 0.23], [0.46, 0.46], color="#999999", lw=0.4)

    tx = 1.42
    ax.text(tx, 1.05, "Ghia et al. (1982) report, at each $Re$:", fontsize=8, va="top", fontweight="bold")
    ax.text(tx, 0.94, "$\\bullet$ $u$ along the vertical and $v$ along\n   the horizontal centreline (green)\n"
            "$\\bullet$ strength $\\psi_{\\min}$ and position of the\n   primary vortex (blue)\n"
            "$\\bullet$ strength, position and size $H$, $V$\n   of each secondary eddy (red)",
            fontsize=7.8, va="top", linespacing=1.45)
    ax.text(tx, 0.28, "BL1, BR1: first corner eddies\nBL2, BR2: smaller eddies inside them\n"
            "TL1: top-left eddy, only at high $Re$\neach eddy turns opposite to its neighbour",
            fontsize=7.8, va="top", color=RED, linespacing=1.45)
    ax.text(tx, -0.3, "sizes drawn for illustration, not to scale", fontsize=7, color=GREY, style="italic")
    fig.savefig(os.path.join(OUT_SCH, "ghia_configuration.pdf")); plt.close(fig)


def make_schematics():
    plt.rcdefaults()
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
                         "savefig.bbox": "tight", "savefig.pad_inches": 0.04, "hatch.linewidth": 0.6})
    for f in (roadmap, stencils_1d, characteristics, five_point, cavity_schematic, ghost_point,
              algorithm_cycle, ghia_configuration):
        f()
    print("schematics written to", OUT_SCH)


# =====================================================================================
# 2. THE JOURNEY: 1D advection, wavenumbers, 2D advection, Burgers, Laplace/Poisson
#    (recomputed with the same updates and parameters as my scripts, vectorised)
# =====================================================================================
def make_journey():
    plt.rcdefaults()
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
                         "axes.linewidth": 0.6, "legend.fontsize": 7.5,
                         "savefig.bbox": "tight", "savefig.pad_inches": 0.03})
    W = 5.8
    EX=dict(color="0.72",lw=2.2); NU=dict(color="#c0392b",lw=0.9)

    # ---------- FTBS (1D_advection_FTBS.py: C=1.005, inflow never updated) ----------
    L=10; nx=1001; dx=L/(nx-1); dt=0.001; nt=500; x=np.linspace(0,L,nx)
    def ftbs(c):
        C=c*dt/dx; u=np.sin(np.pi*x)
        for n in range(nt):
            un=u.copy(); u[1:]=un[1:]-C*(un[1:]-un[:-1])
        return u
    c=10.05; u=ftbs(c); s=np.sin(np.pi*(x-c*nt*dt))
    fig,(ax,az)=plt.subplots(1,2,figsize=(W,2.05),gridspec_kw=dict(width_ratios=[3.1,1]),constrained_layout=True)
    ax.plot(x,s,label="exact",**EX); ax.plot(x,u,label="FTBS, $C=1.005$",**NU)
    xf=c*nt*dt
    ax.axvspan(0,xf,color="#f5d0cc",alpha=0.5,lw=0); ax.text(0.2,-0.95,"region reached by the frozen inflow",fontsize=7.5)
    ax.set_xlim(0,10); ax.set_ylim(-1.15,1.15); ax.set_xlabel("$x$"); ax.set_ylabel("$u$")
    ax.legend(loc="lower right",bbox_to_anchor=(1.0,1.0),frameon=False,ncol=2,borderaxespad=0.1)
    ax.axvspan(4.86,5.12,fill=False,ec="0.35",lw=0.6)
    m=(x>4.86)&(x<5.12); az.plot(x[m],s[m],**EX); az.plot(x[m],u[m],marker=".",ms=2.5,**NU)
    az.axvspan(4.86,xf,color="#f5d0cc",alpha=0.5,lw=0)
    az.set_xlim(4.86,5.12); az.set_ylim(-0.42,0.42); az.set_xticks([4.9,5.0,5.1]); az.set_xlabel("$x$")
    az.set_title("the front, enlarged",fontsize=8)
    fig.savefig(OUT_J+"ftbs_frozen_inflow.pdf"); plt.close(fig)
    wrong=np.abs(u-s)>0.1
    print(f"FTBS: max|u|={np.abs(u).max():.4f}, max error={np.abs(u-s).max():.4f}, fraction of domain with error>0.1: {wrong.mean():.2f}")

    # ---------- FTCS (1D_advection_FTCS.py: C=0.8, both ends frozen) ----------
    C=0.8; u=np.sin(np.pi*x); amp=[]
    for n in range(nt):
        un=u.copy(); u[1:-1]=un[1:-1]-C*(un[2:]-un[:-2])/2; amp.append(np.abs(u).max())
    amp=np.array(amp); n=np.arange(1,nt+1)
    k=slice(300,500); slope=np.polyfit(n[k],np.log10(amp[k]),1)[0]
    print(f"FTCS: final max|u|={amp[-1]:.3e};  measured growth per step (steps 300-500) = {10**slope:.4f}")
    fig,axs=plt.subplots(1,2,figsize=(W,2.1))
    axs[0].semilogy(n,amp,color="#c0392b",lw=1)
    axs[0].set_xlabel("time step $n$"); axs[0].set_ylabel(r"$\max|u|$"); axs[0].set_ylim(0.5,1e52)
    axs[1].plot(x,u,color="#c0392b",lw=0.5); axs[1].set_xlabel("$x$"); axs[1].set_ylabel("$u$ at $n=500$"); axs[1].set_xlim(0,10)
    fig.tight_layout(w_pad=1.5); fig.savefig(OUT_J+"ftcs_blowup.pdf"); plt.close(fig)

    # ---------- numerical wavenumbers (Forier/order1.py, order2.py) ----------
    K=np.linspace(1e-6,np.pi,400)
    fig,axs=plt.subplots(1,2,figsize=(W,2.5))
    a=axs[0]
    a.plot(K,np.ones_like(K),"k",lw=0.7,label="exact")
    a.plot(K,np.sin(K)/K,color="#2980b9",lw=4,alpha=0.35,label="FD, BD, CD-2 (real)")
    a.plot(K,(8*np.sin(K)-np.sin(2*K))/(6*K),color="#27ae60",lw=1.4,label="CD-4 (real)")
    a.plot(K,(1-np.cos(K))/K,color="#e67e22",lw=1,ls="--",label="FD (imag)")
    a.plot(K,-(1-np.cos(K))/K,color="#8e44ad",lw=1,ls="--",label="BD (imag)")
    a.axhline(0,color="0.6",lw=0.4); a.set_xlim(0,np.pi); a.set_xticks([0,np.pi/2,np.pi]); a.set_xticklabels(["0",r"$\pi/2$",r"$\pi$"])
    a.set_xlabel(r"$K=k\Delta x$"); a.set_ylabel(r"$k_{\mathrm{num}}/k$"); a.set_title("(a) first derivative",fontsize=8.5); a.legend(frameon=False,fontsize=6.5,loc="lower left")
    b=axs[1]; S=(np.sin(K/2)/(K/2))**2
    b.plot(K,np.ones_like(K),"k",lw=0.7,label="exact")
    b.plot(K,S,color="#2980b9",lw=1.4,label="CD-2")
    b.plot(K,(30-32*np.cos(K)+2*np.cos(2*K))/(12*K**2),color="#27ae60",lw=1.4,label="CD-4")
    b.plot(K,S*np.cos(K),color="#c0392b",lw=1.2,label="FD, BD (real)")
    b.plot(K,S*np.sin(K),color="#e67e22",lw=1,ls="--",label="FD (imag)")
    b.axhline(0,color="0.6",lw=0.4); b.set_xlim(0,np.pi); b.set_xticks([0,np.pi/2,np.pi]); b.set_xticklabels(["0",r"$\pi/2$",r"$\pi$"])
    b.set_xlabel(r"$K$"); b.set_ylabel(r"$k^2_{\mathrm{num}}/k^2$"); b.set_title("(b) second derivative",fontsize=8.5); b.legend(frameon=False,fontsize=6.5,loc="lower left")
    fig.tight_layout(w_pad=1.4); fig.savefig(OUT_J+"numerical_wavenumbers.pdf"); plt.close(fig)
    Kc=1e-15; print("CD-4 second-derivative expression at K=1e-15 (as in order2.py):",(30-32*np.cos(Kc)+2*np.cos(2*Kc))/(12*Kc**2))

    # ---------- 2D linear advection (2D_advection_linear.py) ----------
    nx=ny=101; dx=dy=0.1; cx=cy=10; dt=0.0004; nt2=500; Cx=cx*dt/dx; Cy=cy*dt/dy
    xx=np.linspace(0,10,nx); X,Y=np.meshgrid(xx,xx); u=np.sin(np.pi*X)*np.sin(np.pi*Y)
    for f in range(nt2):
        uo=u.copy(); u[1:,1:]=uo[1:,1:]-Cx*(uo[1:,1:]-uo[1:,:-1])-Cy*(uo[1:,1:]-uo[:-1,1:])
    T=dt*nt2; ua=np.sin(np.pi*(X-cx*T))*np.sin(np.pi*(Y-cy*T))
    fig,axs=plt.subplots(1,3,figsize=(W,2.05),constrained_layout=True)
    for a,f,tt in zip(axs,[u,ua,u-ua],["numerical","exact","numerical $-$ exact"]):
        im=a.contourf(X,Y,f,levels=np.linspace(-1,1,21),cmap="RdBu_r",extend="both"); a.set_aspect("equal"); a.set_title(tt,fontsize=8.5); a.set_xticks([0,5,10]); a.set_yticks([0,5,10])
    for a in axs[1:]: a.set_yticklabels([])
    fig.colorbar(im,ax=axs,shrink=0.85,pad=0.02,ticks=[-1,-0.5,0,0.5,1],label="$u$")
    fig.savefig(OUT_J+"advection_2d.pdf"); plt.close(fig)
    print(f"2D advection: CFLx+CFLy={Cx+Cy:.2f}; max|error|={np.abs(u-ua).max():.3f}; error>0.1 over {np.mean(np.abs(u-ua)>0.1)*100:.0f}% of grid")

    # ---------- Burgers (1D_advection_nonlinear.py; 1D_advection_diffusion_nonlinear_Burgers_eq.py) ----------
    def burgers(nx,nt,dt,D,signcheck):
        dx=10/(nx-1); x=np.linspace(0,10,nx); u=np.sin(np.pi*x)
        for n in range(nt):
            uo=u.copy(); c=uo[1:-1]
            if signcheck: ux=np.where(c>0,(c-uo[:-2])/dx,(uo[2:]-c)/dx)
            else: ux=(c-uo[:-2])/dx
            u[1:-1]=c-c*dt*ux+D*dt/dx**2*(uo[2:]-2*c+uo[:-2])
        return x,u
    x1,u1=burgers(1001,1000,0.001,0.0,True); x2,u2=burgers(101,1000,0.001,0.1,False)
    print(f"Burgers, nu=0.1, t=1: backward difference everywhere (my script) max {u2.max():.3f} min {u2.min():.3f}")
    _,u3=burgers(101,1000,0.001,0.1,True);   print(f"  same grid, upwinding on the sign of u: max {u3.max():.3f} min {u3.min():.3f}")
    _,u4=burgers(2001,20000,0.00005,0.1,True); print(f"  2001 points, dt=5e-5, sign check:      max {u4.max():.3f} min {u4.min():.3f}")
    for (xx_,uu),t_,fn in [((x1,u1),"inviscid, upwinding on the sign of $u$, $t=1$","burgers_inviscid.pdf"),
                           ((x2,u2),"viscous, $\\nu=0.1$, $t=1$","burgers_viscous.pdf")]:
        fig,a=plt.subplots(figsize=(W,1.9))
        a.plot(xx_,np.sin(np.pi*xx_),**EX,label="$u(x,0)$"); a.plot(xx_,uu,**NU,label="$u(x,1)$")
        a.set_title(t_,fontsize=8.5); a.set_xlim(0,10); a.set_xlabel("$x$"); a.set_ylabel("$u$")
        a.legend(frameon=False,loc="lower left",fontsize=7,ncol=2)
        fig.savefig(OUT_J+fn); plt.close(fig)

    # ---------- Laplace and Poisson ----------
    fig,axs=plt.subplots(1,3,figsize=(W,2.15),constrained_layout=True)
    n=11; s=np.zeros((n,n)); s[-1,:]=1
    for it in range(500):
        so=s.copy(); sn=s.copy(); sn[1:-1,1:-1]=0.25*(so[1:-1,2:]+so[1:-1,:-2]+so[:-2,1:-1]+so[2:,1:-1]); e=np.abs(sn-so).max(); s=sn
        if e<1e-10: break
    g=np.linspace(0,10,n); Xg,Yg=np.meshgrid(g,g); im0=axs[0].contourf(Xg,Yg,s,levels=np.linspace(0,1,21),cmap="viridis"); axs[0].plot(5,5,"w+",ms=6)
    axs[0].set_title(f"(a) Laplace, $11^2$\ncentre {s[5,5]:.10f}",fontsize=7.5)
    n=101; dx=0.1; o=np.zeros((n,n)); o[20:31,20:31]=1
    def jac(k):
        s=np.zeros((n,n))
        for it in range(k):
            sn=s.copy(); sn[1:-1,1:-1]=0.25*(s[:-2,1:-1]+s[2:,1:-1]+s[1:-1,2:]+s[1:-1,:-2]+o[1:-1,1:-1]*dx**2); s=sn
        return s
    g=np.linspace(0,10,n); Xg,Yg=np.meshgrid(g,g); s5=jac(5); sc=jac(20000); lev=np.linspace(0,sc.max(),21)
    axs[1].contourf(Xg,Yg,s5,levels=lev,cmap="viridis"); axs[1].set_title(f"(b) Poisson, 5 sweeps\nmax {s5.max():.4f}",fontsize=7.5)
    im2=axs[2].contourf(Xg,Yg,sc,levels=lev,cmap="viridis"); axs[2].set_title(f"(c) Poisson, 20000 sweeps\nmax {sc.max():.4f}",fontsize=7.5)
    for a in axs: a.set_aspect("equal"); a.set_xticks([0,5,10]); a.set_yticks([0,5,10])
    fig.colorbar(im0,ax=axs[0],shrink=0.8,pad=0.02,ticks=[0,0.5,1])
    fig.colorbar(im2,ax=axs[1:],shrink=0.8,pad=0.02,ticks=[0,0.2,0.4],label="$\\psi$")
    fig.savefig(OUT_J+"laplace_poisson.pdf"); plt.close(fig)

    plt.rcParams.update({"savefig.pad_inches": 0.1})
    # --- identical to r2.py: L=10, nx=1001, c=-1, dt=0.001, T=1, FTFS, inflow at x=L ---
    L=10; nx=1001; dx=L/(nx-1); c=-1; dt=0.001; T=1; nt=int(T/dt); x=np.linspace(0,L,nx)
    fig,axs=plt.subplots(3,1,figsize=(5.8,4.6),sharex=True)
    for ax,f in zip(axs,[2,10,50]):
        u=np.sin(f*np.pi*x)
        for i in range(nt):
            u_new=u.copy()
            u[:-1]=u_new[:-1]-c*(dt/dx)*(u_new[1:]-u_new[:-1])      # same update as the j-loop in r2.py
            u[-1]=np.sin(f*np.pi*(L-c*(i+1)*dt))
        s=np.sin(f*np.pi*(x-c*T))
        m=x<8.5
        amp=np.abs(u[m]).max()
        ax.plot(x,s,color="0.75",lw=0.6,label="exact")
        ax.plot(x,u,color="#c0392b",lw=0.8,label="FTFS")
        ax.set_xlim(0,10); ax.set_ylim(-1.15,1.15); ax.set_yticks([-1,0,1])
        ax.text(0.01,0.9,f"$f={f}$: {200//f} points per wavelength, surviving amplitude {amp:.3f}",transform=ax.transAxes,fontsize=8,va="top",
                bbox=dict(fc="white",ec="none",pad=1))
        print(f"f={f}: max|u| for x<8.5 = {amp:.4f}")
    axs[0].legend(loc="lower right",fontsize=7,frameon=False,ncol=2); axs[-1].set_xlabel("$x$")
    fig.savefig(OUT_J+"r2_frequencies.pdf"); plt.close(fig)

    print("journey figures written to", OUT_J)


# =====================================================================================
# 3. THE 201x201 OVERFLOW RUN (2D_Ghia4.py recomputed; about a minute)
# =====================================================================================
def run(n, nsteps, dt=0.001, nu=0.01, U=1.0):
    h = 1.0 / (n - 1)
    s = np.zeros((n, n)); o = np.zeros((n, n)); u = np.zeros((n, n)); v = np.zeros((n, n))
    err = []
    with np.errstate(all="ignore"):
        for nn in range(nsteps):
            for _ in range(10000):                          # Jacobi, as in the script
                sn = s.copy()
                sn[1:-1, 1:-1] = 0.25 * (s[2:, 1:-1] + s[:-2, 1:-1] + s[1:-1, 2:] + s[1:-1, :-2]
                                         + o[1:-1, 1:-1] * h**2)
                e = np.max(np.abs(sn - s)); s = sn
                if not e >= 1e-6:                            # also stops on NaN
                    break
            u[1:-1, 1:-1] = (s[2:, 1:-1] - s[:-2, 1:-1]) / (2 * h)
            v[1:-1, 1:-1] = -(s[1:-1, 2:] - s[1:-1, :-2]) / (2 * h)
            # wall vorticity (ghost point), same lines as the script
            g = 3 * U * h - 1.5 * s[-1, 1:-1] + 3 * s[-2, 1:-1] - 0.5 * s[-3, 1:-1]
            o[-1, 1:-1] = -(g - 2 * s[-1, 1:-1] + s[-2, 1:-1]) / h**2
            g = -1.5 * s[0, 1:-1] + 3 * s[1, 1:-1] - 0.5 * s[2, 1:-1]
            o[0, 1:-1] = -(g - 2 * s[0, 1:-1] + s[1, 1:-1]) / h**2
            g = -1.5 * s[1:-1, 0] + 3 * s[1:-1, 1] - 0.5 * s[1:-1, 2]
            o[1:-1, 0] = -(g - 2 * s[1:-1, 0] + s[1:-1, 1]) / h**2
            g = -1.5 * s[1:-1, -1] + 3 * s[1:-1, -2] - 0.5 * s[1:-1, -3]
            o[1:-1, -1] = -(g - 2 * s[1:-1, -1] + s[1:-1, -2]) / h**2
            oc = o[1:-1, 1:-1]; uc = u[1:-1, 1:-1]; vc = v[1:-1, 1:-1]
            ox = np.where(uc > 0, (oc - o[1:-1, :-2]) / h, (o[1:-1, 2:] - oc) / h)
            oy = np.where(vc > 0, (oc - o[:-2, 1:-1]) / h, (o[2:, 1:-1] - oc) / h)
            lap = (o[1:-1, 2:] + o[1:-1, :-2] + o[2:, 1:-1] + o[:-2, 1:-1] - 4 * oc) / h**2
            on = o.copy()
            on[1:-1, 1:-1] = oc - dt * uc * ox - dt * vc * oy + nu * dt * lap
            d = np.max(np.abs(on - o)); o = on
            err.append(d)
            if not np.isfinite(d):
                break
    return np.array(err)



def make_overflow():
    plt.rcdefaults()
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
                         "axes.linewidth": 0.6, "legend.fontsize": 7.5,
                         "savefig.bbox": "tight", "savefig.pad_inches": 0.03})
    e201 = run(201, 400)
    e101 = run(101, len(e201) + 20)
    print("201^2: step 0 error", e201[0], "; steps until non-finite:", len(e201),
          "; last finite", e201[np.isfinite(e201)][-1])
    print("101^2: step 0 error", e101[0], "; error at last step", e101[-1])
    for k in [0, 1, 2, 5, 10, 20, 50]:
        if k < len(e201): print(" 201 step", k, e201[k])

    fig, ax = plt.subplots(figsize=(5.2, 2.6))
    k = np.arange(len(e201)); f = np.isfinite(e201)
    ax.semilogy(k[f], e201[f], color="#c0392b", lw=1.3, label=r"$201^2$ ($D = 0.8$)")
    ax.semilogy(np.arange(len(e101)), e101, color="#1f77b4", lw=1.3, label=r"$101^2$ ($D = 0.2$)")
    ax.set_yticks([1e0, 1e50, 1e100, 1e150, 1e200])
    ax.set_ylim(1e-2, 1e240)
    if (~f).any():
        ax.axvline(k[~f][0], color="0.4", ls=":", lw=0.8)
        ax.text(k[~f][0] + 0.5, 1e150, "step %d: overflow,\nerror is inf / NaN" % k[~f][0],
                fontsize=7.5, color="0.3", va="center")
    ax.annotate("step 0: 240,\nas my run printed", xy=(0, e201[0]), xytext=(2.5, 1e60),
                fontsize=7.5, color="#c0392b", arrowprops=dict(arrowstyle="->", color="#c0392b", lw=0.6))
    ax.set_xlabel("time step $n$"); ax.set_ylabel(r"$\max|\omega^{n+1}-\omega^n|$")
    ax.set_title(r"$\mathrm{Re} = 100$, $\Delta t = 0.001$", fontsize=8.5)
    ax.legend(frameon=False, loc="center right")
    fig.savefig(OUT_J+"overflow_201.pdf"); plt.close(fig)

    print("overflow figure written to", OUT_J)


# =====================================================================================
# 4. CAVITY FIGURES, from the saved run folders
# =====================================================================================
# ---- Ghia et al. (1982) centreline data (from compare_ghia.py) -------------
# (v at Re=1000, x=0.9063 is -0.51550 as printed by Ghia; compare_ghia.py has -0.51500)
ghia_y = np.array(
    [
        1.0000,
        0.9766,
        0.9688,
        0.9609,
        0.9531,
        0.8516,
        0.7344,
        0.6172,
        0.5000,
        0.4531,
        0.2813,
        0.1719,
        0.1016,
        0.0703,
        0.0625,
        0.0547,
        0.0000,
    ]
)
ghia_u = {
    100: np.array(
        [
            1.00000,
            0.84123,
            0.78871,
            0.73722,
            0.68717,
            0.23151,
            0.00332,
            -0.13641,
            -0.20581,
            -0.21090,
            -0.15662,
            -0.10150,
            -0.06434,
            -0.04775,
            -0.04192,
            -0.03717,
            0.00000,
        ]
    ),
    400: np.array(
        [
            1.00000,
            0.75837,
            0.68439,
            0.61756,
            0.55892,
            0.29093,
            0.16256,
            0.02135,
            -0.11477,
            -0.17119,
            -0.32726,
            -0.24299,
            -0.14612,
            -0.10338,
            -0.09266,
            -0.08186,
            0.00000,
        ]
    ),
    1000: np.array(
        [
            1.00000,
            0.65928,
            0.57492,
            0.51117,
            0.46604,
            0.33304,
            0.18719,
            0.05702,
            -0.06080,
            -0.10648,
            -0.27805,
            -0.38289,
            -0.29730,
            -0.22220,
            -0.20196,
            -0.18109,
            0.00000,
        ]
    ),
    3200: np.array(
        [
            1.00000,
            0.53236,
            0.48296,
            0.46547,
            0.46101,
            0.34682,
            0.19791,
            0.07156,
            -0.04272,
            -0.08664,
            -0.24427,
            -0.34323,
            -0.41933,
            -0.37827,
            -0.35344,
            -0.32407,
            0.00000,
        ]
    ),
    5000: np.array(
        [
            1.00000,
            0.48223,
            0.46120,
            0.45992,
            0.46036,
            0.33556,
            0.20087,
            0.08183,
            -0.03039,
            -0.07404,
            -0.22855,
            -0.33050,
            -0.40435,
            -0.43643,
            -0.42901,
            -0.41165,
            0.00000,
        ]
    ),
    7500: np.array(
        [
            1.00000,
            0.47244,
            0.47048,
            0.47323,
            0.47167,
            0.34228,
            0.20591,
            0.08342,
            -0.03800,
            -0.07503,
            -0.23176,
            -0.32393,
            -0.38324,
            -0.43025,
            -0.43590,
            -0.43154,
            0.00000,
        ]
    ),
    10000: np.array(
        [
            1.00000,
            0.47221,
            0.47783,
            0.48070,
            0.47804,
            0.34635,
            0.20673,
            0.08344,
            -0.03111,
            -0.07540,
            -0.23186,
            -0.32709,
            -0.38000,
            -0.41657,
            -0.42537,
            -0.42735,
            0.00000,
        ]
    ),
}
ghia_x = np.array(
    [
        1.0000,
        0.9688,
        0.9609,
        0.9531,
        0.9453,
        0.9063,
        0.8594,
        0.8047,
        0.5000,
        0.2344,
        0.2266,
        0.1563,
        0.0938,
        0.0781,
        0.0703,
        0.0625,
        0.0000,
    ]
)
ghia_v = {
    100: np.array(
        [
            0.00000,
            -0.05906,
            -0.07391,
            -0.08864,
            -0.10313,
            -0.16914,
            -0.22445,
            -0.24533,
            0.05454,
            0.17527,
            0.17507,
            0.16077,
            0.12317,
            0.10890,
            0.10091,
            0.09233,
            0.00000,
        ]
    ),
    400: np.array(
        [
            0.00000,
            -0.12146,
            -0.15663,
            -0.19254,
            -0.22847,
            -0.38,
            -0.44993,
            -0.38598,
            0.05186,
            0.30174,
            0.30203,
            0.28124,
            0.22965,
            0.20920,
            0.19713,
            0.18360,
            0.00000,
        ]
    ),
    1000: np.array(
        [
            0.00000,
            -0.21388,
            -0.27669,
            -0.33714,
            -0.39188,
            -0.51550,
            -0.42665,
            -0.31966,
            0.02526,
            0.32235,
            0.33075,
            0.37095,
            0.32627,
            0.30353,
            0.29012,
            0.27485,
            0.00000,
        ]
    ),
    3200: np.array(
        [
            0.00000,
            -0.39017,
            -0.47425,
            -0.52357,
            -0.54053,
            -0.44307,
            -0.37401,
            -0.31184,
            0.00999,
            0.28188,
            0.29030,
            0.37119,
            0.42768,
            0.41906,
            0.40917,
            0.39560,
            0.00000,
        ]
    ),
    5000: np.array(
        [
            0.00000,
            -0.49774,
            -0.55069,
            -0.55408,
            -0.52876,
            -0.41442,
            -0.36214,
            -0.30018,
            0.00945,
            0.27280,
            0.28066,
            0.35368,
            0.42951,
            0.43648,
            0.43329,
            0.42447,
            0.00000,
        ]
    ),
    7500: np.array(
        [
            0.00000,
            -0.53858,
            -0.55216,
            -0.52347,
            -0.48590,
            -0.41050,
            -0.36213,
            -0.30448,
            0.00824,
            0.27348,
            0.28117,
            0.35060,
            0.41824,
            0.43564,
            0.44030,
            0.43979,
            0.00000,
        ]
    ),
    10000: np.array(
        [
            0.00000,
            -0.54302,
            -0.52987,
            -0.49099,
            -0.45863,
            -0.41496,
            -0.36737,
            -0.30719,
            0.00831,
            0.27224,
            0.28003,
            0.35070,
            0.41487,
            0.43124,
            0.43733,
            0.43983,
            0.00000,
        ]
    ),
}



ghia_v[400] = ghia_v[400].copy(); ghia_v[400][5] = -0.23827   # restore printed value; point is excluded below
PSI_GHIA = {100: -0.103423, 400: -0.113909, 1000: -0.117929, 3200: -0.120377, 5000: -0.118966,
            7500: -0.119976, 10000: -0.119731}   # Ghia et al. (1982), Table V



def load_runs():
    """Every run folder */*/ with a parameters.py; fields from pictures/*/solution_fields.npz."""
    out = {}
    for pf in sorted(glob.glob(os.path.join(ROOT, "*", "*", "parameters.py"))):
        folder = os.path.dirname(pf); tag = os.path.basename(folder)
        ns = {}; exec(open(pf).read(), ns)
        npz = glob.glob(os.path.join(folder, "pictures", "*", "solution_fields.npz"))
        d = dict(tag=tag, folder=folder, Re=ns["Re"], n=ns["nx"], dt=ns["dt"])
        if npz:
            z = np.load(npz[0]); d.update({k: z[k] for k in z.files})
        out[tag] = d
    return out


def make_cavity():
    plt.rcdefaults()
    plt.rcParams.update({"font.family":"serif","mathtext.fontset":"cm","font.size":9,"axes.labelsize":9,
        "legend.fontsize":7.5,"xtick.labelsize":8,"ytick.labelsize":8,"axes.linewidth":0.6,"lines.linewidth":1.1,
        "savefig.bbox":"tight","savefig.pad_inches":0.03})
    W = 5.8
    R = load_runs()
    # the numbers each figure needs, per run
    CASES = ["P1","Q1","V1","g3","P2","Q2","V2","P3","Q3","V3","P4","Q4","V4","P5","Q5","V5","P6","Q6","P7","Q7","V6"]
    # g3 = Re 100 on 257x257 ("Grid study/g3"), finished on the internship computer
    res = {}
    for t in CASES:
        d = R[t]; Re, n = d["Re"], d["n"]; h = 1/(n-1); s = d["s"]
        res[t] = dict(Re=Re, n=n, dx=h, psi=float(s.min()), psi_g=PSI_GHIA.get(Re),
                      area=float((s[1:-1,1:-1] > 1e-6).mean()*100),   # % of cavity with psi > 0
                      ratio=Re*h/2)          # nu_num/nu = (U dx/2)/(1/Re): steady state, no time-step part
    GRID_STYLE={51:dict(color="#1f77b4",ls=":"),101:dict(color="#2ca02c",ls="--"),129:dict(color="#d62728",ls="-"),257:dict(color="#9467bd",ls="-")}
    def cl(t):
        d=R[t]; n=d["n"]; u=d["u"][:,n//2].copy(); u[-1]=1.0   # lid value is the imposed boundary condition
        return d["y"],u,d["x"],d["v"][n//2,:]
    def ghia_markers(ax,Re,which):
        if which=="u": ax.plot(ghia_u[Re],ghia_y,"ko",ms=3.2,mfc="white",mew=0.8,label="Ghia et al.",zorder=5)
        else:
            m=np.ones(17,bool)
            if Re==400: m[5]=False; ax.plot(ghia_x[5],-0.23827,"kx",ms=4.5,mew=0.9,label="excluded point",zorder=5)
            ax.plot(ghia_x[m],ghia_v[Re][m],"ko",ms=3.2,mfc="white",mew=0.8,label="Ghia et al.",zorder=5)
    def centreline_panel(fname,cases):
        fig,axs=plt.subplots(2,3,figsize=(W,4.3))
        for c,(Re,tags) in enumerate(cases):
            for t in tags:
                y,u,x,v=cl(t); n=R[t]["n"]
                axs[0,c].plot(u,y,label=f"{n}$\\times${n}",**GRID_STYLE[n]); axs[1,c].plot(x,v,label=f"{n}$\\times${n}",**GRID_STYLE[n])
            ghia_markers(axs[0,c],Re,"u"); ghia_markers(axs[1,c],Re,"v")
            axs[0,c].set_title(f"$Re = {Re}$",fontsize=9)
            axs[0,c].set_xlabel("$u(0.5,y)$"); axs[1,c].set_xlabel("$x$")
            axs[0,c].set_ylabel("$y$" if c==0 else ""); axs[1,c].set_ylabel("$v(x,0.5)$" if c==0 else "")
            for a in axs[:,c]: a.grid(lw=0.3,alpha=0.5); a.axhline(0,lw=0.4,c="gray") if a is axs[1,c] else a.axvline(0,lw=0.4,c="gray")
        axs[0,0].legend(loc="lower right",frameon=False)
        if any(Re==400 for Re,_ in cases): axs[1,[Re for Re,_ in cases].index(400)].legend(loc="lower left",frameon=False)
        fig.tight_layout(pad=0.4,w_pad=0.6,h_pad=0.8); fig.savefig(OUT_C+fname); plt.close(fig)
    centreline_panel("cavity_centreline_low_re.pdf",[(100,["P1","Q1","V1","g3"]),(400,["P2","Q2","V2"]),(1000,["P3","Q3","V3"])])



    # streamfunction topology on the finest available grid
    def psi_panel(ax,t,title):
        d=R[t]; X,Y=np.meshgrid(d["x"],d["y"]); s=d["s"]
        neg=np.linspace(s.min(),0,12)[:-1]
        pos=np.array([1e-5,1e-4,5e-4,1e-3,2e-3,4e-3,8e-3]); pos=pos[pos<s.max()]
        ax.contourf(X,Y,np.where(s>0,1,np.nan),levels=[0.5,1.5],colors=["#f4c7c3"])
        ax.contour(X,Y,s,levels=neg,colors="k",linewidths=0.45,linestyles="solid")
        if len(pos): ax.contour(X,Y,s,levels=pos,colors="#b2182b",linewidths=0.45,linestyles="--")
        j,i=np.unravel_index(np.argmin(s),s.shape); ax.plot(d["x"][i],d["y"][j],"k+",ms=5,mew=0.8)
        ax.set_aspect("equal"); ax.set_xticks([0,0.5,1]); ax.set_yticks([0,0.5,1]); ax.set_title(title,fontsize=8.5)
    fig,axs=plt.subplots(2,3,figsize=(W,4.1))
    for a,(t,ti) in zip(axs.ravel(),[("g3","$Re=100$, 257$^2$"),("V2","$Re=400$, 129$^2$"),("V3","$Re=1000$, 129$^2$"),
                                       ("V4","$Re=3200$, 129$^2$"),("V5","$Re=5000$, 257$^2$"),("V6","$Re=10000$, 257$^2$")]):
        psi_panel(a,t,ti)
    fig.tight_layout(pad=0.3,w_pad=0.4,h_pad=0.6); fig.savefig(OUT_C+"cavity_streamfunction_finest.pdf"); plt.close(fig)

    # streamlines with direction for the most extreme case
    d=R["P7"]; X,Y=np.meshgrid(d["x"],d["y"]); u=d["u"].copy(); u[-1,1:-1]=1.0
    fig,ax=plt.subplots(figsize=(2.9,2.9))
    ax.streamplot(X,Y,u,d["v"],density=1.1,color="k",linewidth=0.5,arrowsize=0.7)
    ax.set_aspect("equal"); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_xticks([0,0.5,1]); ax.set_yticks([0,0.5,1])
    ax.annotate("",xy=(0.95,1.035),xytext=(0.05,1.035),arrowprops=dict(arrowstyle="->",lw=1.0),annotation_clip=False)
    ax.text(0.5,1.06,"lid, $U = 1$",ha="center",va="bottom",fontsize=8)
    fig.savefig(OUT_C+"cavity_re10000_51_streamlines.pdf"); plt.close(fig)

    # psi_min convergence
    fig,axs=plt.subplots(1,2,figsize=(W,2.5))
    for Re,tr,mk in [(400,["P2","Q2","V2"],"s"),(1000,["P3","Q3","V3"],"o")]:
        dx=np.array([res[t]["dx"] for t in tr]); e=np.array([abs(res[t]["psi"]-res[t]["psi_g"])/abs(res[t]["psi_g"]) for t in tr])
        axs[0].loglog(dx,e*100,mk+"-",ms=4,mfc="white",label=f"$Re={Re}$")
    xx=np.array([0.0065,0.023]); axs[0].loglog(xx,1300*xx,"k--",lw=0.6); axs[0].text(0.0105,10.3,"slope 1",fontsize=7.5,rotation=24)
    axs[0].set_xlabel("$\\Delta x$"); axs[0].set_ylabel("$|\\psi_{\\min}-\\psi_{\\min}^{\\mathrm{Ghia}}|\\,/\\,|\\psi_{\\min}^{\\mathrm{Ghia}}|$  (%)")
    axs[0].legend(frameon=False); axs[0].grid(lw=0.3,alpha=0.5,which="both")
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter
    axs[0].yaxis.set_major_locator(FixedLocator([5,10,20,30])); axs[0].yaxis.set_major_formatter(FixedFormatter(["5","10","20","30"])); axs[0].yaxis.set_minor_formatter(NullFormatter())
    axs[0].xaxis.set_major_locator(FixedLocator([0.008,0.01,0.02])); axs[0].xaxis.set_major_formatter(FixedFormatter(["0.008","0.01","0.02"])); axs[0].xaxis.set_minor_formatter(NullFormatter())
    Res=[100,400,1000,3200,5000,7500,10000]
    for n,tags in [(51,["P1","P2","P3","P4","P5","P6","P7"]),(101,["Q1","Q2","Q3","Q4","Q5","Q6","Q7"]),(129,["V1","V2","V3","V4"]),(257,["V5","V6"])]:
        re=[res[t]["Re"] for t in tags]; ar=[res[t]["area"] for t in tags]
        axs[1].semilogx(re,ar,"o"+GRID_STYLE[n]["ls"],color=GRID_STYLE[n]["color"],ms=3.5,mfc="white",label=f"{n}$\\times${n}")
    axs[1].semilogx([res["g3"]["Re"]],[res["g3"]["area"]],"o",color=GRID_STYLE[257]["color"],ms=3.5,mfc="white")   # Re=100 on 257^2
    axs[1].set_xlabel("$Re$"); axs[1].set_ylabel("area with $\\psi>0$  (% of cavity)"); axs[1].legend(frameon=False,loc="upper left"); axs[1].grid(lw=0.3,alpha=0.5,which="both")
    fig.tight_layout(pad=0.4,w_pad=1.2); fig.savefig(OUT_C+"cavity_convergence_and_topology.pdf"); plt.close(fig)

    # ---- RMSE-free replacements -------------------------------------------------
    gs=np.loadtxt(os.path.join(ROOT, "Grig study", "grid_independence_summary.csv"),delimiter=",",skiprows=1)
    mk={100:"o",400:"s",1000:"^",3200:"D",5000:"v",7500:"<",10000:">"}
    PSI_G=PSI_GHIA
    fig,axs=plt.subplots(1,2,figsize=(W,2.55))
    ax=axs[0]
    for t,r in res.items():
        if r["Re"] not in PSI_G: continue
        e=abs(r["psi"]-PSI_G[r["Re"]])/abs(PSI_G[r["Re"]])*100
        ax.loglog(r["ratio"],e,mk[r["Re"]],color=GRID_STYLE[r["n"]]["color"],ms=4,mfc="white",mew=0.9)
    ax.axvline(1,c="gray",lw=0.6,ls="--")
    for Re in [100,400,1000,3200,5000,7500,10000]: ax.plot([],[],mk[Re],color="k",mfc="white",ms=4,label=f"$Re={Re}$")
    ax.legend(frameon=False,fontsize=6.5,loc="lower right")
    ax.set_xlabel("$\\nu_{\\mathrm{num}}/\\nu$"); ax.set_ylabel("$\\psi_{\\min}$ deviation from Ghia (%)"); ax.grid(lw=0.3,alpha=0.5,which="both")
    ax.set_title("(a)",fontsize=8.5)
    ax=axs[1]
    mx_val=[]; mx_u=[]
    for k,Re in enumerate([100,400,1000,3200,5000,7500,10000]):
        d=R["Q"+str(k+1)]; iv=np.arange(1,16); iu=np.arange(1,16)
        if Re==400: iv=iv[iv!=5]
        ev=np.interp(ghia_x,d["x"],d["v"][50,:])[iv]-ghia_v[Re][iv]; mx_val.append(np.abs(ev).max())
        eu=np.interp(ghia_y,d["y"],d["u"][:,50])[iu]-ghia_u[Re][iu]; mx_u.append(np.abs(eu).max())
    Rl=[100,400,1000,3200,5000,7500,10000]
    ax.semilogy(gs[:,0],gs[:,4],"o-",color="#1f77b4",ms=4,mfc="white",label="$51^2$ vs $101^2$, $v$")
    ax.semilogy(gs[:,0],gs[:,2],"o:",color="#1f77b4",ms=3.2,mfc="#1f77b4",lw=1.0,label="$51^2$ vs $101^2$, $u$")
    ax.semilogy(Rl,mx_val,"s--",color="#d62728",ms=4,mfc="white",label="$101^2$ vs Ghia, $v$")
    ax.semilogy(Rl,mx_u,"s:",color="#d62728",ms=3.2,mfc="#d62728",lw=1.0,label="$101^2$ vs Ghia, $u$")
    ax.set_xscale("log"); ax.set_xlabel("$Re$"); ax.set_ylabel("largest centreline difference"); ax.set_ylim(3e-3,2.0); ax.legend(frameon=False,fontsize=6.6,loc="upper left"); ax.grid(lw=0.3,alpha=0.5,which="both")
    ax.set_title("(b)",fontsize=8.5)
    fig.tight_layout(pad=0.4,w_pad=1.2); fig.savefig(OUT_C+"cavity_viscosity_and_gridstudy.pdf"); plt.close(fig)
    print("max |v-Ghia| on 101:",[round(x,4) for x in mx_val])

    # results at a glance: vortex strength and centreline error against Re
    grids = [(51, ["P1","P2","P3","P4","P5","P6","P7"]), (101, ["Q1","Q2","Q3","Q4","Q5","Q6","Q7"]),
             (129, ["V1","V2","V3","V4"]), (257, ["V5","V6"])]
    def cl_err(t):
        d = R[t]; Re, n = d["Re"], d["n"]
        iu = np.arange(1, 16); iv = iu[iu != 5] if Re == 400 else iu
        du = np.interp(ghia_y, d["y"], d["u"][:, n // 2])[iu] - ghia_u[Re][iu]
        dv = np.interp(ghia_x, d["x"], d["v"][n // 2, :])[iv] - ghia_v[Re][iv]
        return max(np.abs(du).max(), np.abs(dv).max())
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.6))
    for a in axs:
        a.axvspan(60, 200, color="#e3f1e3", lw=0); a.axvspan(200, 2000, color="#fdf1dc", lw=0)
        a.axvspan(2000, 16000, color="#f9e0dd", lw=0); a.set_xscale("log"); a.set_xlim(60, 16000)
    rg = sorted(PSI_GHIA)
    axs[0].plot(rg, [-PSI_GHIA[r] for r in rg], "k-o", ms=4, lw=1.4, mfc="k", label="Ghia et al.", zorder=5)
    for n, tags in grids:
        st = GRID_STYLE[n]
        axs[0].plot([R[t]["Re"] for t in tags], [-R[t]["s"].min() for t in tags], "o", ls=st["ls"],
                    color=st["color"], ms=3.5, mfc="white", label=f"{n}$\\times${n}")
        axs[1].plot([R[t]["Re"] for t in tags], [cl_err(t) for t in tags], "o", ls=st["ls"],
                    color=st["color"], ms=3.5, mfc="white")
    st = GRID_STYLE[257]   # Re = 100 on 257^2 (g3): a point of its own, not joined to Re >= 5000
    axs[0].plot([100], [-R["g3"]["s"].min()], "o", color=st["color"], ms=3.5, mfc="white")
    axs[1].plot([100], [cl_err("g3")], "o", color=st["color"], ms=3.5, mfc="white")
    axs[0].set_ylabel(r"vortex strength $|\psi_{\min}|$"); axs[0].set_ylim(0, 0.162)
    axs[1].set_ylabel(r"max. centreline error in $u$, $v$"); axs[1].set_yscale("log"); axs[1].set_ylim(3e-3, 4.0)
    for a, tt in zip(axs, ["(a)", "(b)"]):
        a.set_xlabel("$Re$"); a.set_title(tt, fontsize=8.5); a.grid(lw=0.3, alpha=0.5, which="major")
        for x0, lab in [(110, "agrees"), (620, "too weak,\nfirst order"), (5400, "wrong\nflow")]:
            a.text(x0, 0.97, lab, transform=a.get_xaxis_transform(), ha="center", va="top", fontsize=6.8, color="0.3")
    axs[0].legend(frameon=False, loc="lower left", fontsize=6.8)
    fig.tight_layout(pad=0.4, w_pad=1.2); fig.savefig(OUT_C + "cavity_summary.pdf"); plt.close(fig)
    for t in ["Q1","V1","g3","V2","V3","V4","V5","V6"]:
        print(f"  {t}: Re={int(R[t]['Re']):6d}  max centreline error {cl_err(t):.4f}")

    # vorticity fields (129x129)
    plt.rcdefaults()
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
                         "axes.linewidth": 0.6, "legend.fontsize": 7.5,
                         "savefig.bbox": "tight", "savefig.pad_inches": 0.03})
    fig, axs = plt.subplots(1, 3, figsize=(W, 2.15), constrained_layout=True)
    lev = np.array([-5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 5], float)
    for ax, (tag, Re) in zip(axs, [("V1", 100), ("V2", 400), ("V3", 1000)]):
        d = R[tag]
        X, Y = np.meshgrid(d["x"], d["y"])
        o = np.clip(d["o"], -5.5, 5.5)
        cf = ax.contourf(X, Y, o, levels=lev, cmap="RdBu_r", extend="both")
        ax.contour(X, Y, d["o"], levels=lev, colors="k", linewidths=0.3)
        ax.set_aspect("equal"); ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
        if tag != "V1": ax.set_yticklabels([])
        ax.set_title(f"$Re={Re}$, $129^2$", fontsize=8.5)
        print(f"Re={Re}: omega range {d['o'].min():.1f} .. {d['o'].max():.1f}; "
              f"interior |omega| max {np.abs(d['o'][1:-1,1:-1]).max():.1f}")
    cb = fig.colorbar(cf, ax=axs, shrink=0.9, pad=0.01, ticks=[-5, -3, -1, 0, 1, 3, 5])
    cb.set_label(r"$\omega$ (clipped at $\pm5$)")
    fig.savefig(OUT_C+"cavity_vorticity.pdf"); plt.close(fig)



    # ---- streamlines on the finest grid at every Re (black: psi<0, red: psi>0) ----------------
    from matplotlib.colors import ListedColormap, BoundaryNorm
    cmap2 = ListedColormap(["k", "#b2182b"]); norm2 = BoundaryNorm([-0.5, 0.5, 1.5], 2)
    fig, axs = plt.subplots(2, 4, figsize=(W, 3.35))
    finest = [("g3", 100, 257), ("V2", 400, 129), ("V3", 1000, 129), ("V4", 3200, 129),
              ("V5", 5000, 257), ("Q6", 7500, 101), ("V6", 10000, 257)]
    for ax, (t, Re, n) in zip(axs.ravel(), finest):
        d = R[t]; X, Y = np.meshgrid(d["x"], d["y"]); u = d["u"].copy(); u[-1, 1:-1] = 1.0
        sgn = (d["s"] > 0).astype(float)
        ax.streamplot(X, Y, u, d["v"], density=0.9, color=sgn, cmap=cmap2, norm=norm2,
                      linewidth=0.45, arrowsize=0.45)
        j, i = np.unravel_index(np.argmin(d["s"]), d["s"].shape); ax.plot(d["x"][i], d["y"][j], "k+", ms=5, mew=0.9)
        ax.set_aspect("equal"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1]); ax.tick_params(labelsize=6.5)
        ax.set_title(f"$Re={Re}$, ${n}^2$", fontsize=7.8)
    ax = axs[1, 3]; ax.axis("off")
    ax.plot([0.08, 0.28], [0.62, 0.62], "k-", lw=0.9); ax.text(0.33, 0.62, "$\\psi<0$: clockwise,\ndriven by the lid", va="center", fontsize=6.8)
    ax.plot([0.08, 0.28], [0.32, 0.32], "-", color="#b2182b", lw=0.9); ax.text(0.33, 0.32, "$\\psi>0$:\ncounter-rotating", va="center", fontsize=6.8)
    ax.plot([0.18], [0.06], "k+", ms=5, mew=0.9); ax.text(0.33, 0.06, "$\\psi_{\\min}$", va="center", fontsize=6.8)
    ax.set_xlim(0, 1); ax.set_ylim(-0.05, 0.85)
    fig.tight_layout(pad=0.3, w_pad=0.3, h_pad=0.5); fig.savefig(OUT_C + "cavity_streamlines_finest.pdf"); plt.close(fig)

    # ---- vorticity at Re >= 3200 on the finest grids (same scale as the low-Re figure) -------
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9})
    fig, axs = plt.subplots(1, 3, figsize=(W, 2.15), constrained_layout=True)
    lev = np.array([-5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 5], float)
    for ax, (tag, Re, n) in zip(axs, [("V4", 3200, 129), ("V5", 5000, 257), ("V6", 10000, 257)]):
        d = R[tag]; X, Y = np.meshgrid(d["x"], d["y"])
        cf = ax.contourf(X, Y, np.clip(d["o"], -5.5, 5.5), levels=lev, cmap="RdBu_r", extend="both")
        ax.contour(X, Y, d["o"], levels=lev, colors="k", linewidths=0.3)
        ax.set_aspect("equal"); ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
        if tag != "V4": ax.set_yticklabels([])
        ax.set_title(f"$Re={Re}$, ${n}^2$", fontsize=8.5)
    cb = fig.colorbar(cf, ax=axs, shrink=0.9, pad=0.01, ticks=[-5, -3, -1, 0, 1, 3, 5])
    cb.set_label(r"$\omega$ (clipped at $\pm5$)")
    fig.savefig(OUT_C + "cavity_vorticity_high_re.pdf"); plt.close(fig)

    # ---- Poisson residual |lap_h(psi) + omega| of the final fields, Re = 100 -------------------
    from matplotlib.colors import LogNorm
    fig, axs = plt.subplots(1, 3, figsize=(W, 2.15), constrained_layout=True)
    for ax, (tag, n, lab) in zip(axs, [("Q1", 101, "settled"), ("V1", 129, "never settled"), ("g3", 257, "never settled")]):
        d = R[tag]; s_, o_ = d["s"], d["o"]; h = d["x"][1] - d["x"][0]
        r = np.abs((s_[1:-1, 2:] - 2*s_[1:-1, 1:-1] + s_[1:-1, :-2]) / h**2
                   + (s_[2:, 1:-1] - 2*s_[1:-1, 1:-1] + s_[:-2, 1:-1]) / h**2 + o_[1:-1, 1:-1])
        X, Y = np.meshgrid(d["x"][1:-1], d["y"][1:-1])
        cf = ax.contourf(X, Y, np.clip(r, 1e-7, 0.99), levels=np.logspace(-7, 0, 15), norm=LogNorm(1e-7, 1), cmap="magma_r")
        j, i = np.unravel_index(np.argmax(r), r.shape); ax.plot(X[j, i], Y[j, i], "c^", ms=4, mec="k", mew=0.4)
        ax.set_aspect("equal"); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
        if tag != "Q1": ax.set_yticklabels([])
        ax.set_title(f"${n}^2$ ({lab})", fontsize=8.5)
        print(f"residual map {tag}: max {r.max():.2e}, median {np.median(r):.1e}")
    cb = fig.colorbar(cf, ax=axs, shrink=0.9, pad=0.01, ticks=[1e-7, 1e-5, 1e-3, 1e-1])
    cb.set_label(r"$|\nabla_h^2\psi + \omega|$")
    fig.savefig(OUT_C + "cavity_poisson_residual.pdf"); plt.close(fig)

    print("cavity figures written to", OUT_C)


if __name__ == "__main__":
    parts = sys.argv[1:] or ["schematics", "journey", "overflow", "cavity"]
    if "schematics" in parts: make_schematics()
    if "journey" in parts: make_journey()
    if "overflow" in parts: make_overflow()
    if "cavity" in parts: make_cavity()
    print("done.")
