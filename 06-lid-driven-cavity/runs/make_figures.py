"""Generate the figures committed to figures/.

Kept as a script, not a notebook, so every figure in the repository can be
regenerated from the saved .npz runs by one command.
"""
import glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from src.ghia_compare import centreline_u, centreline_v, load_ghia

os.makedirs("figures", exist_ok=True)
runs = {}
for f in sorted(glob.glob("runs/Re*_N*.npz")):
    r = dict(np.load(f))
    runs[int(r["Re"])] = r

# --- centreline validation ------------------------------------------------
yc, gu = load_ghia("u")
xc, gv = load_ghia("v")
fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
cols = {100: "tab:blue", 400: "tab:orange", 1000: "tab:red"}
for Re, r in sorted(runs.items()):
    y, u = centreline_u(r); x, v = centreline_v(r)
    ax[0].plot(u, y, "-", color=cols[Re], label=f"present, Re={Re}")
    ax[0].plot(gu[Re], yc, "o", color=cols[Re], mfc="none", ms=5)
    ax[1].plot(x, v, "-", color=cols[Re], label=f"present, Re={Re}")
    ax[1].plot(xc, gv[Re], "o", color=cols[Re], mfc="none", ms=5)
ax[0].set_xlabel("u"); ax[0].set_ylabel("y"); ax[0].set_title("u on x = 0.5")
ax[1].set_xlabel("x"); ax[1].set_ylabel("v"); ax[1].set_title("v on y = 0.5")
for a in ax:
    a.grid(alpha=0.3); a.legend(fontsize=8)
fig.suptitle("Centreline velocities: lines = present solver, circles = Ghia et al. (1982)")
fig.tight_layout(); fig.savefig("figures/centreline_validation.png", dpi=150)
print("figures/centreline_validation.png")

# --- streamfunction -------------------------------------------------------
fig, ax = plt.subplots(1, len(runs), figsize=(4.2*len(runs), 4))
for a, (Re, r) in zip(np.atleast_1d(ax), sorted(runs.items())):
    lv = np.array([-0.1,-0.09,-0.07,-0.05,-0.03,-0.01,-1e-3,-1e-5,
                   1e-8,1e-6,1e-5,5e-5,1e-4,2.5e-4])
    a.contour(r["X"], r["Y"], r["psi"], levels=lv, colors="k", linewidths=0.7)
    a.set_title(f"Re = {Re}"); a.set_aspect("equal")
    a.set_xlabel("x"); a.set_ylabel("y")
fig.suptitle("Streamfunction contours (levels chosen to expose corner vortices)")
fig.tight_layout(); fig.savefig("figures/streamfunction.png", dpi=150)
print("figures/streamfunction.png")

# --- error vs numerical viscosity ratio -----------------------------------
from src.ghia_compare import compare
from src.numerical_viscosity import numerical_viscosity
dts = {100: 1e-3, 400: 2.5e-3, 1000: 4e-3}
Res, ratios, rmses = [], [], []
for Re, r in sorted(runs.items()):
    _, _, _, m = compare(r, "u")
    nn = numerical_viscosity(1.0, float(r["dx"]), dts[Re])
    Res.append(Re); ratios.append(nn/(1.0/Re)); rmses.append(m["RMSE"])
fig, a = plt.subplots(figsize=(5.5, 4))
a.plot(ratios, rmses, "o-", color="tab:red")
for R, x, y in zip(Res, ratios, rmses):
    a.annotate(f"Re={R}", (x, y), textcoords="offset points", xytext=(6, -4))
a.axvline(1.0, ls="--", c="gray", lw=1)
a.text(1.05, max(rmses)*0.5, r"$\nu_{num}=\nu$", color="gray", fontsize=9)
a.set_xlabel(r"$\nu_{num}/\nu$"); a.set_ylabel("u-centreline RMSE vs Ghia")
a.set_title("Error tracks the added numerical viscosity"); a.grid(alpha=0.3)
fig.tight_layout(); fig.savefig("figures/error_vs_numerical_viscosity.png", dpi=150)
print("figures/error_vs_numerical_viscosity.png")
