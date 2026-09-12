"""Verification for the centreline extraction and comparison machinery.

The error tables are only trustworthy if the tooling that produces them is.

  1. interpolation reproduces an exactly-known function
  2. the centreline extractors pick the right row/column
  3. relative error excludes near-zero reference values (and why)
  4. a run compared against itself gives exactly zero error

Run:
    python3 tests/test_ghia_compare.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.ghia_compare import (centreline_u, centreline_v, error_metrics,
                              interpolate_to, load_ghia)


def test_interpolation_is_accurate():
    """Interpolating a known function onto Ghia's coordinates.

    Ghia's locations are 129-grid points (0.9766, 0.9688, ...), which are NOT
    grid points of a 101 grid: 0.9766/0.01 = 97.66. Nearest-point sampling
    would misplace the sample by up to dx/2, and near the lid -- where du/dy is
    steepest -- that is a large error in value, silently blamed on the scheme.
    """
    have = np.linspace(0, 1, 101)
    f = lambda y: np.sin(3 * np.pi * y) * y**2
    targets, _ = load_ghia("u")

    got = interpolate_to(targets, have, f(have))
    err = np.abs(got - f(targets)).max()
    print(f"  max interpolation error on a smooth test function = {err:.2e}")
    assert err < 2e-3, "linear interpolation is losing too much"

    # nearest-point sampling, for contrast
    idx = np.abs(have[None, :] - targets[:, None]).argmin(axis=1)
    err_near = np.abs(f(have)[idx] - f(targets)).max()
    print(f"  same, using nearest grid point instead        = {err_near:.2e}")
    print(f"  interpolation is {err_near/err:.1f}x better: PASS\n")


def test_centrelines_pick_the_right_line():
    """x = 0.5 and y = 0.5 are exact grid points on any odd-N grid."""
    N = 101
    xs = np.linspace(0, 1, N)
    X, Y = np.meshgrid(xs, xs)
    r = {"N": N, "X": X, "Y": Y, "u": Y.copy(), "v": X.copy()}

    y, u = centreline_u(r)
    x, v = centreline_v(r)
    print(f"  vertical centreline at x = {X[0, (N-1)//2]:.4f}")
    print(f"  horizontal centreline at y = {Y[(N-1)//2, 0]:.4f}")
    assert abs(X[0, (N - 1) // 2] - 0.5) < 1e-12
    assert abs(Y[(N - 1) // 2, 0] - 0.5) < 1e-12
    assert np.allclose(u, y) and np.allclose(v, x)
    print("  both centrelines land exactly on 0.5: PASS\n")


def test_relative_error_excludes_near_zero():
    """Why relative error alone is misleading on this benchmark.

    Ghia's Table I contains u = 0.00332 at y = 0.7344 -- the profile crosses
    zero there. A physically negligible absolute discrepancy of 0.003 at that
    point reads as ~90% relative error and would dominate any average.
    """
    ref = np.array([1.0, 0.5, 0.00332, -0.5])
    mine = ref + 0.003

    no_floor = error_metrics(mine, ref, rel_floor=0.0)
    with_floor = error_metrics(mine, ref, rel_floor=0.05)
    print(f"  no floor  : mean rel = {no_floor['rel_pct']:7.2f}%  "
          f"({no_floor['n_rel']}/{no_floor['n_total']} pts)")
    print(f"  floor 0.05: mean rel = {with_floor['rel_pct']:7.2f}%  "
          f"({with_floor['n_rel']}/{with_floor['n_total']} pts)")
    assert no_floor["rel_pct"] > 5 * with_floor["rel_pct"]
    assert with_floor["n_rel"] == 3
    print("  the near-zero point is excluded, as it must be: PASS\n")


def test_self_comparison_is_zero():
    """Comparing anything against itself must give identically zero."""
    ref = np.array([1.0, 0.3, -0.2, -0.05])
    m = error_metrics(ref, ref)
    print(f"  L_inf = {m['L_inf']:.1e}   RMSE = {m['RMSE']:.1e}   "
          f"rel = {m['rel_pct']:.1e}%")
    assert m["L_inf"] == 0.0 and m["RMSE"] == 0.0
    print("  identical inputs give zero error: PASS\n")


if __name__ == "__main__":
    tests = [
        ("interpolation accuracy", test_interpolation_is_accurate),
        ("centreline extraction", test_centrelines_pick_the_right_line),
        ("relative error floor", test_relative_error_excludes_near_zero),
        ("self-comparison", test_self_comparison_is_zero),
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
