"""Verification for the Poisson module.

  1. method of manufactured solutions -- the solver gets the RIGHT answer
  2. second-order convergence of the discretisation
  3. all three methods converge to the same solution
  4. iteration-count scaling: Jacobi O(N^2), SOR O(N)
  5. residual vs change-per-sweep -- why the stopping test matters

Run:
    python3 tests/test_poisson.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.poisson import (SOLVERS, gauss_seidel, jacobi, optimal_omega,
                         residual, sor)


def manufactured(N, L=1.0):
    """s = sin(pi x) sin(pi y), so laplacian(s) = -2 pi^2 s exactly.

    METHOD OF MANUFACTURED SOLUTIONS. Rather than finding a physical problem
    with a known answer, PICK the answer first, substitute it into the
    equation, and whatever comes out is the source term f that makes it true.
    Then the solver must reproduce the function you started from.

    This is the strongest verification available: it tests the discretisation,
    the boundary handling and the iteration together, against an exact answer,
    for a problem of your choosing. Chosen here to be zero on all four
    boundaries, so homogeneous Dirichlet conditions are exact.
    """
    dx = L / (N - 1)
    x = np.linspace(0, L, N)
    X, Y = np.meshgrid(x, x)
    exact = np.sin(np.pi * X) * np.sin(np.pi * Y)
    f = -2.0 * np.pi**2 * exact
    return np.zeros((N, N)), f, exact, dx


def test_manufactured_solution():
    N = 41
    s0, f, exact, dx = manufactured(N)
    print("  " + "solver".ljust(16) + "iters".rjust(8) + "max|s - exact|".rjust(18))
    for name, fn in SOLVERS.items():
        s, it, ok = fn(s0, f, dx, tol=1e-8)
        err = np.abs(s - exact).max()
        print("  " + name.ljust(16) + f"{it:8d}" + f"{err:18.3e}")
        assert ok, f"{name} did not converge"
        assert err < 1e-3, f"{name} converged to the wrong answer"
    print("  every solver reproduces the manufactured solution: PASS\n")


def test_second_order_convergence():
    """Discretisation error must fall as dx^2.

    Distinct from iterative convergence. Even a perfectly converged iteration
    leaves the error of the five-point STENCIL, which module 01 derived as
    -dx^2/12 * u'''' -- second order. Iterating harder cannot beat it; only a
    finer grid can.
    """
    print("  " + "N".rjust(5) + "max error".rjust(14) + "ratio".rjust(9) + "   order")
    prev = None
    for N in (11, 21, 41, 81):
        s0, f, exact, dx = manufactured(N)
        s, _, _ = jacobi(s0, f, dx, tol=1e-10)
        err = np.abs(s - exact).max()
        if prev is None:
            print("  " + f"{N:5d}" + f"{err:14.3e}")
        else:
            p = np.log2(prev / err)
            print("  " + f"{N:5d}" + f"{err:14.3e}" + f"{prev/err:9.2f}" + f"   {p:.2f}")
            assert abs(p - 2.0) < 0.1, f"order {p:.2f}, expected 2"
        prev = err
    print("  discretisation is second order: PASS\n")


def test_all_methods_agree():
    """Three different iterations, one linear system, one answer."""
    N = 31
    s0, f, _, dx = manufactured(N)
    sols = {n: fn(s0, f, dx, tol=1e-10)[0] for n, fn in SOLVERS.items()}
    ref = sols["Jacobi"]
    for name, s in sols.items():
        d = np.abs(s - ref).max()
        print(f"  {name.ljust(16)} max difference from Jacobi = {d:.3e}")
        assert d < 1e-8, f"{name} disagrees with Jacobi"
    print("  all three solve the same system: PASS\n")


def test_iteration_scaling():
    """Jacobi costs O(N^2) sweeps; SOR with optimal omega costs O(N).

    Not a constant-factor speed-up -- a different complexity class. Doubling
    the grid roughly quadruples Jacobi's work and merely doubles SOR's.
    """
    print("  " + "N".rjust(5) + "omega_opt".rjust(11)
          + "Jacobi".rjust(10) + "Gauss-Seidel".rjust(14) + "SOR".rjust(8))
    counts = {"Jacobi": [], "SOR": []}
    grids = (11, 21, 41)
    for N in grids:
        s0, f, _, dx = manufactured(N)
        row = []
        for name, fn in (("Jacobi", jacobi), ("Gauss-Seidel", gauss_seidel), ("SOR", sor)):
            _, it, _ = fn(s0, f, dx, tol=1e-6)
            row.append(it)
            if name in counts:
                counts[name].append(it)
        print("  " + f"{N:5d}" + f"{optimal_omega(N):11.4f}"
              + f"{row[0]:10d}" + f"{row[1]:14d}" + f"{row[2]:8d}")

    jac = np.log2(counts["Jacobi"][-1] / counts["Jacobi"][-2])
    s_r = np.log2(counts["SOR"][-1] / counts["SOR"][-2])
    print(f"\n  scaling exponent per grid doubling:  Jacobi {jac:.2f}   SOR {s_r:.2f}")
    assert jac > 1.7, "Jacobi should scale close to N^2"
    assert s_r < 1.4, "SOR should scale close to N"
    print("  Jacobi is O(N^2), SOR is O(N): PASS\n")


def test_residual_vs_change():
    """Stopping on the change per sweep badly overestimates convergence.

    For Jacobi, change = residual * dx^2 / 4 exactly. So the residual is
    4/dx^2 times larger than the change -- and that factor GROWS as the grid
    refines, meaning the error is worst precisely when you refine to get a
    better answer.
    """
    print("  " + "N".rjust(5) + "dx".rjust(9) + "change".rjust(12)
          + "residual".rjust(12) + "ratio".rjust(10) + "4/dx^2".rjust(10))
    for N in (11, 21, 41, 81):
        L = 1.0
        dx = L / (N - 1)
        f = np.zeros((N, N))
        s = np.zeros((N, N))
        s[-1, :] = 1.0                      # lid
        for _ in range(100000):
            so = s.copy()
            s[1:-1, 1:-1] = 0.25 * (so[1:-1, 2:] + so[1:-1, :-2]
                                    + so[2:, 1:-1] + so[:-2, 1:-1])
            ch = np.abs(s - so).max()
            if ch < 1e-5:
                break
        res = np.abs(residual(s, f, dx)).max()
        print("  " + f"{N:5d}" + f"{dx:9.4f}" + f"{ch:12.3e}"
              + f"{res:12.3e}" + f"{res/ch:10.1f}" + f"{4/dx**2:10.1f}")
        assert abs(res / ch - 4 / dx**2) / (4 / dx**2) < 0.01, "ratio off theory"
    print("  the gap is exactly 4/dx^2, and grows as the grid refines: PASS\n")


if __name__ == "__main__":
    tests = [
        ("manufactured solution", test_manufactured_solution),
        ("second-order convergence", test_second_order_convergence),
        ("all methods agree", test_all_methods_agree),
        ("iteration scaling", test_iteration_scaling),
        ("residual vs change per sweep", test_residual_vs_change),
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
