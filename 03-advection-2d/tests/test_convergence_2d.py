"""Verification for the 2D linear advection module.

Five checks:
  1. order of accuracy for upwind and Lax-Wendroff
  2. the cross-derivative term is NECESSARY, not optional
  3. the stability limit is Cx + Cy <= 1, not Cx <= 1 and Cy <= 1
  4. conservation
  5. x and y are treated symmetrically (catches transposed indexing)

Run:
    python3 tests/test_convergence_2d.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.initial_conditions_2d import (exact_solution_2d, gaussian_2d,
                                       make_grid_2d, sine_2d, square_2d)
from src.schemes_2d import (SCHEMES_2D, lax_wendroff_2d,
                            lax_wendroff_2d_no_cross, upwind_2d)
from src.solver_2d import l2_error_2d, solve_2d

np.seterr(all="ignore")

L, T = 1.0, 0.5


def _order(step, grids, Cx, Cy, cx, cy):
    """Convergence study on sine_2d.

    The initial condition must be EXACTLY periodic, or the measurement is
    meaningless. A Gaussian with sigma = 0.15 is 3.9e-03 at the domain edge;
    np.mod wraps it regardless, creating a kink no scheme can represent. That
    kink is a fixed-size error which does not shrink with dx, so the measured
    order decays as the grid refines: 1.97, 1.87, 1.63, 1.30 -- a scheme that
    is genuinely second order looking progressively worse the harder you try.

    sine_2d is periodic by construction, so there is no mismatch at all.
    """
    errs = []
    for N in grids:
        X, Y, u, tf = solve_2d(step, sine_2d, N, Cx, Cy, T, cx, cy, L)
        ue = exact_solution_2d(sine_2d, X, Y, L, cx, cy, tf)
        _, _, dx = make_grid_2d(N, L)
        errs.append(l2_error_2d(u, ue, dx))
    return errs, np.log2(errs[-2] / errs[-1])


def test_order_of_accuracy():
    """Diagonal advection, cx = cy = 1."""
    grids = (40, 80, 160, 320)
    print("  " + "scheme".ljust(20)
          + "".join(f"N={n}".rjust(12) for n in grids) + "     p")
    print("  " + "-" * 76)

    failures = []
    for name, step, expect in (("Upwind-2D", upwind_2d, 1),
                               ("Lax-Wendroff-2D", lax_wendroff_2d, 2)):
        errs, p = _order(step, grids, 0.25, 0.25, 1.0, 1.0)
        print("  " + name.ljust(20) + "".join(f"{e:12.3e}" for e in errs) + f"  {p:5.2f}")
        if abs(p - expect) > 0.15:
            failures.append(f"{name}: p={p:.2f}, expected {expect}")

    assert not failures, "; ".join(failures)
    print("  both converge at their formal order: PASS\n")


def test_cross_term_is_necessary():
    """Dropping the u_xy term degrades LW to first order on diagonals.

    In 2D the time-Taylor expansion of u_tt produces u_xx, u_yy AND a mixed
    u_xy term. Omitting the last one is a common mistake, and it is invisible
    if you only ever test flow aligned with an axis -- which is exactly why
    this test uses BOTH cases.
    """
    grids = (40, 80, 160)

    _, p_axis = _order(lax_wendroff_2d_no_cross, grids, 0.25, 0.0, 1.0, 0.0)
    _, p_diag = _order(lax_wendroff_2d_no_cross, grids, 0.25, 0.25, 1.0, 1.0)
    _, p_full = _order(lax_wendroff_2d, grids, 0.25, 0.25, 1.0, 1.0)

    print(f"  LW without cross term, flow along x : p = {p_axis:.2f}")
    print(f"  LW without cross term, flow diagonal: p = {p_diag:.2f}")
    print(f"  LW WITH    cross term, flow diagonal: p = {p_full:.2f}")

    assert p_axis > 1.8, "axis-aligned case should still be second order"
    assert p_diag < 1.5, "diagonal case should degrade without the cross term"
    assert p_full > 1.8, "full LW should be second order on diagonals"
    print("  the cross term is required for second order on diagonals: PASS\n")


def test_stability_limit_is_the_sum():
    """The 2D CFL condition is Cx + Cy <= 1, NOT Cx <= 1 and Cy <= 1.

    Tested on a checkerboard, which is pure theta_x = theta_y = pi -- the mode
    the Von Neumann analysis identifies as worst case. A smooth Gaussian
    contains almost no energy at that wavenumber, so it stays bounded well past
    the true limit and would hide the instability entirely. The initial
    condition has to be chosen to excite the mode being tested.
    """
    N, n_steps = 64, 200
    idx = np.arange(N)
    u0 = ((idx[None, :] + idx[:, None]) % 2).astype(float) * 2.0 - 1.0

    print("  " + "Cx=Cy".ljust(8) + "Cx+Cy".rjust(8)
          + "max|u| after 200 steps".rjust(26) + "   verdict")
    for C in (0.45, 0.49, 0.50, 0.51, 0.55):
        u = u0.copy()
        for _ in range(n_steps):
            u = upwind_2d(u, C, C)
        m = np.abs(u).max()
        verdict = "stable" if m <= 1.0 + 1e-9 else "UNSTABLE"
        print("  " + f"{C:.2f}".ljust(8) + f"{2*C:.2f}".rjust(8)
              + f"{m:26.4e}" + f"   {verdict}")

        if 2 * C <= 1.0:
            assert m <= 1.0 + 1e-9, f"should be stable at Cx+Cy={2*C}"
        else:
            assert m > 10.0, f"should be unstable at Cx+Cy={2*C}"

    print("  boundary is exactly at Cx + Cy = 1: PASS\n")


def test_conservation():
    N, C = 64, 0.25
    X, Y, dx = make_grid_2d(N, L)
    mass0 = np.sum(gaussian_2d(X, Y, L, sigma=0.15)) * dx**2

    for name in ("Upwind-2D", "Lax-Wendroff-2D"):
        _, _, u, _ = solve_2d(SCHEMES_2D[name], gaussian_2d, N, C, C, 2.0,
                              1.0, 1.0, L, sigma=0.15)
        drift = abs(np.sum(u) * dx**2 - mass0) / mass0
        print(f"  {name.ljust(20)} relative drift = {drift:.3e}")
        assert drift < 1e-12, f"{name} does not conserve mass"

    print("  mass conserved to machine precision: PASS\n")


def test_x_y_symmetry():
    """Advecting in +x then transposing must equal transposing then +y.

    A transposed index (u[i,j] where u[j,i] was meant) is the classic 2D bug.
    It survives every symmetric test, so this one is deliberately asymmetric.
    """
    N = 64
    X, Y, _ = make_grid_2d(N, L)
    u0 = square_2d(X, Y, L)

    ux = u0.copy()
    for _ in range(40):
        ux = upwind_2d(ux, 0.3, 0.0)

    uy = u0.T.copy()
    for _ in range(40):
        uy = upwind_2d(uy, 0.0, 0.3)

    diff = np.abs(ux.T - uy).max()
    print(f"  max difference between x-advection and transposed y-advection = {diff:.3e}")
    assert diff < 1e-14, "x and y directions are not treated symmetrically"
    print("  x and y handled identically: PASS\n")


if __name__ == "__main__":
    tests = [
        ("order of accuracy", test_order_of_accuracy),
        ("the cross term is necessary", test_cross_term_is_necessary),
        ("stability limit is Cx + Cy <= 1", test_stability_limit_is_the_sum),
        ("conservation", test_conservation),
        ("x-y symmetry", test_x_y_symmetry),
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
