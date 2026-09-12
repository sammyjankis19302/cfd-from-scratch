"""Run the full parameter study and write .npz files for later comparison.

Usage:
    python3 runs/run_study.py           # all cases
    python3 runs/run_study.py 1000 201  # one case

Each case saves runs/Re{Re}_N{N}.npz containing psi, omega, u, v, X, Y and the
run diagnostics. The original solver saved only PNGs, so a finished run could
not be revisited without repeating it.

Timesteps are chosen near the stability limit for each case. Both constraints
must hold (module 04):

    advective:  dt <= dx / (|u| + |v|)
    diffusive:  dt <= 0.5 dx^2 / (2 nu)

At Re = 100 the diffusive limit binds; at Re = 1000 the advective one does.
Using a smaller dt than necessary is not free -- module 06's theory note shows
nu_num grows as dt shrinks.
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.cavity import save, solve_cavity

CASES = [
    # Re,    N,   dt,     T
    (100,   101, 1.0e-3, 60.0),
    (100,   201, 5.0e-4, 60.0),
    (400,   101, 2.5e-3, 80.0),
    (400,   201, 1.0e-3, 80.0),
    (1000,  101, 4.0e-3, 150.0),
    (1000,  201, 2.0e-3, 150.0),
]


def run(Re, N, dt, T):
    t0 = time.perf_counter()
    r = solve_cavity(N=N, Re=Re, dt=dt, T=T, tol_steady=1e-7, verbose=True)
    el = time.perf_counter() - t0
    os.makedirs("runs", exist_ok=True)
    save(r, f"runs/Re{Re}_N{N}.npz")
    steady = r["final_change"] < 1e-7
    print(f"Re={Re} N={N}: {el:7.1f}s  steps={r['steps']}  t={r['t_final']:.2f}  "
          f"psi_min={r['psi'].min():.6f}  poisson_res={r['poisson_residual']:.2e}  "
          f"d_omega={r['final_change']:.2e}  {'STEADY' if steady else 'NOT STEADY'}")
    if not steady:
        print("   ^ increase T: the numbers below are contaminated by "
              "incomplete convergence, not just numerical viscosity")
    return r


if __name__ == "__main__":
    if len(sys.argv) == 3:
        Re, N = int(sys.argv[1]), int(sys.argv[2])
        case = next(c for c in CASES if c[0] == Re and c[1] == N)
        run(*case)
    else:
        for c in CASES:
            run(*c)
