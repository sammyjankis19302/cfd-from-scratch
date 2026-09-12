"""Regenerate the figures in figures/ from src/spectral.py."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from src.spectral import FIRST, SECOND, group_velocity

os.makedirs("figures", exist_ok=True)
K = np.linspace(1e-6, np.pi, 2000)
c = {"BD-1":"tab:blue","FD-1":"tab:green","CD-2":"tab:orange","CD-4":"tab:red"}

fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
for n, f in FIRST.items():
    s = f(K)
    ax[0].plot(K/np.pi, s.real, color=c[n], label=n)
    ax[1].plot(K/np.pi, s.imag, color=c[n], label=n)
ax[0].plot(K/np.pi, np.ones_like(K), "k--", lw=1, label="exact")
ax[1].axhline(0, color="k", ls="--", lw=1)
ax[0].set_title(r"Re$(k_{num}/k)$  —  dispersion")
ax[1].set_title(r"Im$(k_{num}/k)$  —  dissipation")
for n in ("BD-1","CD-2","CD-4"):
    ax[2].plot(K/np.pi, group_velocity(n, K), color=c[n], label=n)
ax[2].axhline(0, color="k", lw=1); ax[2].axhline(1, color="k", ls="--", lw=1)
ax[2].axvline(0.5, color="gray", ls=":", lw=1)
ax[2].text(0.52, -1.3, "4 points per\nwavelength", fontsize=8, color="gray")
ax[2].set_title(r"$V_g/c$  —  group velocity")
for a in ax:
    a.set_xlabel(r"$k\Delta x/\pi$"); a.grid(alpha=0.3); a.legend(fontsize=8)
fig.suptitle("First derivative: below 4 PPW, CD-2 carries energy upstream")
fig.tight_layout(); fig.savefig("figures/first_derivative.png", dpi=150)
print("figures/first_derivative.png")

fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.3))
c2 = {"CD-2":"tab:orange","CD-4":"tab:red","FD-2":"tab:green","BD-2":"tab:blue"}
for n, f in SECOND.items():
    s = f(K)
    ax[0].plot(K/np.pi, s.real, color=c2[n], label=n)
    ax[1].plot(K/np.pi, s.imag, color=c2[n], label=n)
ax[0].plot(K/np.pi, np.ones_like(K), "k--", lw=1, label="exact")
ax[1].axhline(0, color="k", ls="--", lw=1)
ax[0].set_title(r"Re$(k_{num}^2/k^2)$"); ax[1].set_title(r"Im$(k_{num}^2/k^2)$")
for a in ax:
    a.set_xlabel(r"$k\Delta x/\pi$"); a.grid(alpha=0.3); a.legend(fontsize=8)
fig.suptitle("Second derivative: one-sided stencils acquire an imaginary part")
fig.tight_layout(); fig.savefig("figures/second_derivative.png", dpi=150)
print("figures/second_derivative.png")
