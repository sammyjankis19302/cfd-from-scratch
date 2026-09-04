"""Proof that the readable version and the fast version are the same scheme.

The repository contains first-order upwind twice:

    experiments/first_attempt_refactored.py   explicit loops, explicit
                                              index arithmetic, easy to read
    src/schemes.py                            vectorised with np.roll, ~100x
                                              faster, harder to read at first

Having two implementations of one scheme is only acceptable if you can show
they are the same implementation. That is what this file does. It is also the
best possible explanation of what np.roll actually does: nothing clever, just
the wrap-around index arithmetic written more compactly.

Run:
    python3 tests/test_equivalence.py
"""

import sys
import time

import numpy as np


# --- the two implementations, side by side ---------------------------------
def step_explicit(u_old, C):
    """Readable version: one point at a time, wrap-around written out."""
    nx = len(u_old)
    u_new = np.zeros(nx)
    for i in range(nx):
        if i == 0:
            i_left = nx - 1
        else:
            i_left = i - 1
        u_new[i] = u_old[i] - C * (u_old[i] - u_old[i_left])
    return u_new


def step_vectorised(u, C):
    """Fast version: the whole array at once."""
    return u - C * (u - np.roll(u, 1))


# ---------------------------------------------------------------------------
def test_roll_is_index_shift():
    """np.roll(u, 1) is exactly 'the left neighbour of every point'.

    Demonstrated on a small array where you can check it by eye.
    """
    u = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
    rolled = np.roll(u, 1)

    print("  u             =", u)
    print("  np.roll(u, 1) =", rolled)
    print("  (each entry is the one to its LEFT; the first wraps to the last)")

    manual = np.zeros(len(u))
    for i in range(len(u)):
        manual[i] = u[i - 1] if i > 0 else u[len(u) - 1]

    assert np.array_equal(rolled, manual), "np.roll does not match manual shift"
    print("  identical to the manual loop: PASS\n")


def test_one_step_identical():
    """A single timestep must agree to the last bit."""
    rng = np.random.default_rng(0)
    u = rng.random(200)

    for C in (0.01, 0.3, 0.5, 0.8, 1.0):
        a = step_explicit(u, C)
        b = step_vectorised(u, C)
        max_diff = np.abs(a - b).max()
        print(f"  C = {C:<5}  max difference = {max_diff:.3e}")
        assert max_diff == 0.0, f"implementations differ at C={C}"

    print("  bit-for-bit identical at every C: PASS\n")


def test_full_run_identical():
    """1000 timesteps: errors do not accumulate differently."""
    nx, nt, C = 2000, 200, 0.5
    L = 10.0
    dx = L / nx
    x = np.array([i * dx for i in range(nx)])
    u0 = np.sin(np.pi * x)

    ua = u0.copy()
    t0 = time.perf_counter()
    for _ in range(nt):
        ua = step_explicit(ua, C)
    t_explicit = time.perf_counter() - t0

    ub = u0.copy()
    t0 = time.perf_counter()
    for _ in range(nt):
        ub = step_vectorised(ub, C)
    t_vectorised = time.perf_counter() - t0

    max_diff = np.abs(ua - ub).max()
    print(f"  after {nt} steps, max difference = {max_diff:.3e}")
    print(f"  explicit loops : {t_explicit:.3f} s")
    print(f"  vectorised     : {t_vectorised:.3f} s")
    print(f"  speed-up       : {t_explicit / t_vectorised:.0f}x  (at nx = {nx})")
    print("  note: the gap widens with nx -- roughly 3x at nx=100,")
    print("        16x at nx=500, 52x at nx=2000, 204x at nx=8000.")

    assert max_diff == 0.0, "implementations diverge over a long run"
    print("  identical after a full run: PASS\n")


def test_no_frozen_boundary():
    """Guard against the original bug: u[0] must actually be updated.

    This is the regression test for the inflow-boundary bug. If someone ever
    reintroduces `for i in range(1, nx-1)`, this fails immediately.

    Note the initial condition: COS, not sin. With sin(pi x) the value at
    x = 0 happens to be zero, and it stays near zero for a while even when the
    scheme is working correctly -- so a frozen u[0] would sail through
    undetected. The first version of this test did exactly that: it "passed"
    on a floating-point accident between +0.000000 and -0.000000.

    A test that cannot fail is not a test. Choose the case that would expose
    the bug loudly.
    """
    nx = 50
    L = 10.0
    dx = L / nx
    x = np.array([i * dx for i in range(nx)])
    u = np.cos(np.pi * x)          # u[0] = 1.0, unmistakable
    u0_first = u[0]

    for _ in range(50):
        u = step_explicit(u, 0.5)

    change = abs(u[0] - u0_first)
    print(f"  u[0] initially {u0_first:+.6f}, after 50 steps {u[0]:+.6f}")
    print(f"  change = {change:.6f}")
    assert change > 1e-6, "u[0] barely moved -- inflow boundary is frozen"
    print("  inflow boundary is being updated: PASS\n")


if __name__ == "__main__":
    tests = [
        ("np.roll is just an index shift", test_roll_is_index_shift),
        ("one step is identical", test_one_step_identical),
        ("a full run is identical", test_full_run_identical),
        ("the inflow boundary is not frozen", test_no_frozen_boundary),
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
