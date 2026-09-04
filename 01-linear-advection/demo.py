"""Minimal driver: advect a Gaussian once round the domain with every scheme.

Run from the module root:

    python3 demo.py

Prints the L2 error against the exact solution. Nothing is plotted yet -- the
point of this script is to confirm the pieces talk to each other correctly
before any figures exist.
"""

import numpy as np

from src.initial_conditions import exact_solution, gaussian, make_grid
from src.schemes import FORMAL_ORDER, SCHEMES
from src.solver import l2_error, solve

np.seterr(all="ignore")  # unstable schemes will overflow; that is the point

C, T, c, L, N = 0.5, 1.0, 1.0, 1.0, 400

print(f"Advecting a Gaussian once round the domain: N={N}, C={C}, T={T}\n")
print(f"{'scheme':18s}{'L2 error':>12s}   formal order")
print("-" * 46)

for name, step in SCHEMES.items():
    x, u, t_final = solve(step, gaussian, N, C, T, c, L, sigma=0.1)
    u_exact = exact_solution(gaussian, x, L, c, t_final, sigma=0.1)
    _, dx = make_grid(N, L)
    err = l2_error(u, u_exact, dx)
    order = FORMAL_ORDER[name]
    label = "unstable" if order is None else str(order)
    shown = "diverged" if not np.isfinite(err) or err > 1e3 else f"{err:.4e}"
    print(f"{name:18s}{shown:>12s}   {label}")
