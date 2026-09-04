"""Verification: does each scheme converge at the order the theory claims?

This is the file that separates "my code runs" from "my code is correct".
Four independent checks, each ruling out a different class of error:

  1. measured order of accuracy      -> the discretisation is what I think
  2. exactness at C = 1              -> the upwind stencil is exactly right
  3. conservation                    -> nothing is leaking out of the domain
  4. unstable schemes really diverge -> the stability theory is not decorative

Run:
    python3 tests/test_convergence.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.initial_conditions import exact_solution, gaussian, make_grid, square
from src.schemes import FORMAL_ORDER, SCHEMES
from src.solver import l2_error, solve

np.seterr(all="ignore")  # the unstable schemes overflow on purpose


def test_order_of_accuracy():
    """Halve dx, see how much the error falls. That ratio IS the order.

    A first-order scheme halves its error; a second-order scheme quarters it.
    p = log2(e_coarse / e_fine).

    Note sigma = 0.1, not 0.05. With a narrower Gaussian the coarse grids do
    not resolve the initial condition at all, and the measured order comes out
    around 0.7 for a scheme that is genuinely first order. A convergence table
    that disagrees with theory usually means the test is under-resolved, not
    that the scheme is broken.
    """
    C, T, c, L = 0.5, 1.0, 1.0, 1.0
    grids = (200, 400, 800, 1600)
    sigma = 0.1

    print("  " + "scheme".ljust(16) + "".join(f"N={n}".rjust(12) for n in grids) + "     p  formal")
    print("  " + "-" * 76)

    failures = []
    for name, step in SCHEMES.items():
        order = FORMAL_ORDER[name]
        if order is None:
            continue

        errs = []
        for N in grids:
            x, u, t_final = solve(step, gaussian, N, C, T, c, L, sigma=sigma)
            u_exact = exact_solution(gaussian, x, L, c, t_final, sigma=sigma)
            _, dx = make_grid(N, L)
            errs.append(l2_error(u, u_exact, dx))

        p = np.log2(errs[-2] / errs[-1])
        print("  " + name.ljust(16) + "".join(f"{e:12.3e}" for e in errs)
              + f"  {p:5.2f}  {order}")

        if abs(p - order) > 0.15:
            failures.append(f"{name}: measured p={p:.2f}, expected {order}")

    for f in failures:
        print(f"  FAIL {f}")
    assert not failures, "; ".join(failures)
    print("  every scheme converges at its formal order: PASS\n")


def test_exact_at_cfl_one():
    """At C = 1 upwind is EXACT, even for a discontinuity.

    At C = 1 the update is u_new[i] = u[i-1]: the profile is shifted exactly
    one cell per step, which is precisely what the PDE does. Numerical
    viscosity nu = c*dx*(1-C)/2 vanishes identically.

    This is a much sharper test than a convergence rate, because it predicts an
    exact number rather than a trend. Run with a square wave: if there were any
    dissipation or dispersion at all, a discontinuity would expose it instantly.
    """
    N, C, T, c, L = 200, 1.0, 1.0, 1.0, 1.0

    x, u, t_final = solve(SCHEMES["FTBS"], square, N, C, T, c, L)
    u_exact = exact_solution(square, x, L, c, t_final)
    err = np.abs(u - u_exact).max()

    print(f"  FTBS, square wave, C = 1.0, max error = {err:.3e}")
    assert err < 1e-12, f"upwind is not exact at C=1 (error {err:.3e})"
    print("  exact to machine precision: PASS\n")


def test_conservation():
    """sum(u)*dx must not drift.

    The advection equation conserves the integral of u. Every conservative
    scheme here should hold it to machine precision over a long run -- errors
    can redistribute u, but they must not create or destroy it.
    """
    N, C, T, c, L = 200, 0.8, 5.0, 1.0, 1.0

    print("  " + "scheme".ljust(16) + "relative drift in sum(u)*dx")
    failures = []
    for name in ("FTBS", "Lax-Friedrichs", "Lax-Wendroff", "BD-2"):
        x, dx = make_grid(N, L)
        mass0 = np.sum(gaussian(x, L, sigma=0.1)) * dx
        _, u, _ = solve(SCHEMES[name], gaussian, N, C, T, c, L, sigma=0.1)
        drift = abs(np.sum(u) * dx - mass0) / mass0

        print("  " + name.ljust(16) + f"{drift:.3e}")
        if drift > 1e-12:
            failures.append(f"{name}: drift {drift:.3e}")

    assert not failures, "; ".join(failures)
    print("  mass conserved to machine precision: PASS\n")


def test_unstable_schemes_diverge():
    """FTCS and FTFS must blow up. If they do not, the test problem is wrong.

    A stability analysis that never gets confirmed is just algebra. This is the
    positive control: it proves the test setup is actually capable of
    detecting instability.
    """
    N, C, T, c, L = 200, 0.5, 1.0, 1.0, 1.0

    for name in ("FTCS", "FTFS"):
        _, u, _ = solve(SCHEMES[name], gaussian, N, C, T, c, L, sigma=0.1)
        blew_up = (not np.all(np.isfinite(u))) or np.abs(u).max() > 1e3
        print(f"  {name.ljust(8)} max|u| = {np.abs(u).max():.3e}  diverged = {blew_up}")
        assert blew_up, f"{name} should be unstable but stayed bounded"

    print("  unconditionally unstable schemes do diverge: PASS\n")


def test_lw_and_bd2_coincide_at_half():
    """Two different schemes give identical error at exactly C = 1/2.

    Lax-Wendroff's leading dispersive coefficient goes as (1 - C^2);
    Beam-Warming's goes as (1 - C)(2 - C). Their ratio is (1+C)/(2-C), which
    equals 1 at C = 1/2. So the schemes must coincide there and nowhere else.

    This started as a suspected bug and turned out to be a prediction.
    """
    N, T, c, L = 400, 1.0, 1.0, 1.0

    print("  " + "C".ljust(6) + "LW".rjust(12) + "BD-2".rjust(12)
          + "ratio".rjust(10) + "predicted".rjust(11))
    for C in (0.3, 0.5, 0.8):
        errs = []
        for name in ("Lax-Wendroff", "BD-2"):
            x, u, t_final = solve(SCHEMES[name], gaussian, N, C, T, c, L, sigma=0.1)
            ue = exact_solution(gaussian, x, L, c, t_final, sigma=0.1)
            _, dx = make_grid(N, L)
            errs.append(l2_error(u, ue, dx))

        ratio = errs[0] / errs[1]
        predicted = (1 + C) / (2 - C)
        print(f"  {C:<6}{errs[0]:12.4e}{errs[1]:12.4e}{ratio:10.4f}{predicted:11.4f}")
        assert abs(ratio - predicted) < 1e-3, f"ratio off at C={C}"

    print("  measured ratio matches (1+C)/(2-C): PASS\n")


if __name__ == "__main__":
    tests = [
        ("order of accuracy", test_order_of_accuracy),
        ("exactness at C = 1", test_exact_at_cfl_one),
        ("conservation", test_conservation),
        ("unstable schemes diverge", test_unstable_schemes_diverge),
        ("LW and BD-2 coincide at C = 1/2", test_lw_and_bd2_coincide_at_half),
    ]

    failed = 0
    for name, fn in tests:
        print(f"[{name}]")
        try:
            fn()
        except AssertionError as exc:
            print(f"  FAIL: {exc}\n")
            failed += 1

    if failed:
        print(f"{failed} test(s) failed.")
        sys.exit(1)
    print("All tests passed.")
