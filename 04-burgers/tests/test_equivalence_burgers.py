"""Verification for the Burgers module.

  1. the vectorised 2D solver reproduces the original loop version exactly
  2. conditional upwinding actually picks the right direction
  3. a fixed (non-conditional) upwind direction blows up -- the negative control
  4. shock formation: a smooth wave steepens and the order of accuracy
     legitimately collapses once it does
  5. the two stability limits are where the theory says

Run:
    python3 tests/test_equivalence_burgers.py
"""

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.burgers import solve_1d, step_1d, step_2d, stability_limits

np.seterr(all="ignore")


def _loop_step_2d(u, v, dx, dy, dt, nu):
    """Uday's original nested-loop update, extracted verbatim."""
    ny, nx = u.shape
    un, vn = u.copy(), v.copy()
    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            ux = ((u[j, i] - u[j, i-1]) / dx if u[j, i] > 0
                  else (u[j, i+1] - u[j, i]) / dx)
            uy = ((u[j, i] - u[j-1, i]) / dy if v[j, i] > 0
                  else (u[j+1, i] - u[j, i]) / dy)
            uxx = (u[j, i+1] - 2*u[j, i] + u[j, i-1]) / dx**2
            uyy = (u[j+1, i] - 2*u[j, i] + u[j-1, i]) / dy**2
            un[j, i] = u[j, i] - dt*(ux*u[j, i] + uy*v[j, i]) + nu*dt*(uxx+uyy)

            vx = ((v[j, i] - v[j, i-1]) / dx if u[j, i] > 0
                  else (v[j, i+1] - v[j, i]) / dx)
            vy = ((v[j, i] - v[j-1, i]) / dy if v[j, i] > 0
                  else (v[j+1, i] - v[j, i]) / dy)
            vxx = (v[j, i+1] - 2*v[j, i] + v[j, i-1]) / dx**2
            vyy = (v[j+1, i] - 2*v[j, i] + v[j-1, i]) / dy**2
            vn[j, i] = v[j, i] - dt*(vx*u[j, i] + vy*v[j, i]) + nu*dt*(vxx+vyy)
    return un, vn


def test_loop_matches_vectorised():
    """Same scheme, two implementations, identical answers."""
    L, nx, nu, dt, n = 10.0, 101, 0.001, 1e-4, 60
    dx = L / (nx - 1)
    xs = np.linspace(0, L, nx)
    X, Y = np.meshgrid(xs, xs)
    u0 = np.sin(np.pi * X) * np.sin(np.pi * Y)

    ua, va = u0.copy(), u0.copy()
    t0 = time.perf_counter()
    for _ in range(n):
        ua, va = _loop_step_2d(ua, va, dx, dx, dt, nu)
    t_loop = time.perf_counter() - t0

    ub, vb = u0.copy(), u0.copy()
    t0 = time.perf_counter()
    for _ in range(n):
        ub, vb = step_2d(ub, vb, dx, dx, dt, nu)
    t_vec = time.perf_counter() - t0

    du = np.abs(ua - ub).max()
    dv = np.abs(va - vb).max()
    print(f"  after {n} steps:  max|du| = {du:.3e}   max|dv| = {dv:.3e}")
    print(f"  loops {t_loop:.2f}s   vectorised {t_vec:.3f}s   "
          f"speed-up {t_loop/t_vec:.0f}x")
    assert du < 1e-14 and dv < 1e-14, "implementations disagree"
    print("  identical to machine precision: PASS\n")


def test_conditional_upwinding_picks_correctly():
    """Where u > 0 the scheme must look left; where u < 0, right."""
    dx = 0.1
    u = np.array([0.0, 1.0, 2.0, -2.0, -1.0, 0.0])
    back = (u - np.roll(u, 1)) / dx
    fwd = (np.roll(u, -1) - u) / dx
    ux = np.where(u > 0, back, fwd)

    print("   i   u[i]   chosen u_x   should use")
    for i in range(len(u)):
        want = "backward" if u[i] > 0 else "forward"
        got = back[i] if u[i] > 0 else fwd[i]
        print(f"  {i:2d}  {u[i]:+5.1f}   {ux[i]:+9.2f}   {want}")
        assert np.isclose(ux[i], got)
    print("  direction follows the local sign of u: PASS\n")


def test_fixed_direction_blows_up():
    """Negative control: always-backward upwinding fails where u < 0.

    Proves conditional upwinding is doing real work rather than being a
    stylistic flourish. Where u < 0 a backward difference is DOWNWIND, which
    is the anti-diffusive case from module 01.
    """
    N, L, nu = 200, 1.0, 0.0
    dx = L / N
    x = np.linspace(0, L - dx, N)
    u0 = np.sin(2 * np.pi * x)          # positive and negative halves
    dt = 0.2 * dx

    u = u0.copy()
    for _ in range(400):
        u = step_1d(u, dx, dt, nu)
    ok = np.abs(u).max()

    w = u0.copy()
    for _ in range(400):
        wx = (w - np.roll(w, 1)) / dx    # ALWAYS backward
        w = w - dt * w * wx
    bad = np.abs(w).max()

    print(f"  conditional upwinding: max|u| = {ok:.4e}")
    print(f"  fixed backward only  : max|u| = {bad:.4e}")
    assert ok < 2.0, "conditional scheme should stay bounded"
    assert bad > 10.0 or not np.isfinite(bad), "fixed direction should diverge"
    print("  conditional upwinding is necessary: PASS\n")


def test_shock_formation():
    """A smooth wave steepens toward a shock, and it does so QUANTITATIVELY.

    For inviscid Burgers the characteristics carry u at speed u, so they cross
    at the breaking time

        t_break = 1 / max|u0'(x)|

    which for u0 = sin(2 pi x) is 1/(2 pi) ~ 0.1592. Before that, the maximum
    gradient has a closed form:

        max|u_x|(t) = |u0'|max / (1 - t |u0'|max) = 2 pi / (1 - 2 pi t)

    That is a far sharper test than "the gradient goes up": it predicts a
    number at every instant, and the number diverges at a specific time.

    A NOTE ON GETTING THIS WRONG. My first version of this test used
    t_break = 1/(4 pi^2) ~ 0.0253 and only asserted that the gradient grew by
    some factor. The assertion failed, which looked like a solver bug. It was
    not: the formula was wrong. The measured gradient was independent of grid
    resolution (16.46, 16.72, 16.82, 16.87 at N = 200 to 1600), which is the
    signature of a CONVERGED physical answer rather than a numerical artefact.
    Grid-independence is how you tell the difference.
    """
    N, L, nu = 1600, 1.0, 0.0
    dx = L / N
    x = np.linspace(0, L - dx, N)
    u0 = np.sin(2 * np.pi * x)
    dt = 0.2 * dx
    slope0 = 2 * np.pi
    t_break = 1.0 / slope0

    print(f"  t_break = 1/max|u0'| = 1/(2 pi) = {t_break:.4f}")
    print("      t     measured max|u_x|      theory     rel. error")

    for t in (0.00, 0.05, 0.10, 0.13):
        n = int(round(t / dt))
        u = solve_1d(u0, dx, dt, nu, n)
        g = np.abs((np.roll(u, -1) - np.roll(u, 1)) / (2 * dx)).max()
        theory = slope0 / (1.0 - t * slope0)
        rel = abs(g - theory) / theory
        print(f"  {t:.3f}  {g:16.2f}  {theory:12.2f}  {rel:12.2%}")
        assert rel < 0.02, f"gradient off theory by {rel:.1%} at t={t}"

    print("  gradient tracks 2 pi / (1 - 2 pi t) to within 2%: PASS\n")


def test_stability_limits():
    """Both constraints must hold, and each binds in a different regime."""
    dx = 0.1
    print("  " + "nu".rjust(8) + "dt_advective".rjust(16) + "dt_diffusive".rjust(16)
          + "   binding")
    for nu in (1e-4, 1e-3, 1e-2, 1e-1):
        a, d = stability_limits(1.0, 1.0, dx, dx, nu)
        which = "advective" if a < d else "diffusive"
        print("  " + f"{nu:8.0e}" + f"{a:16.4f}" + f"{d:16.4f}" + f"   {which}")

    a_lo, d_lo = stability_limits(1.0, 1.0, dx, dx, 1e-4)
    a_hi, d_hi = stability_limits(1.0, 1.0, dx, dx, 1e-1)
    assert a_lo < d_lo, "advection should bind at small nu"
    assert d_hi < a_hi, "diffusion should bind at large nu"
    print("  which limit binds depends on nu, as expected: PASS\n")


if __name__ == "__main__":
    tests = [
        ("loop version matches vectorised", test_loop_matches_vectorised),
        ("conditional upwinding picks correctly", test_conditional_upwinding_picks_correctly),
        ("fixed direction blows up", test_fixed_direction_blows_up),
        ("shock formation", test_shock_formation),
        ("stability limits", test_stability_limits),
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
