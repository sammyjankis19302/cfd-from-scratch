"""Quantifying what first-order upwind costs the cavity solver.

  1. nu_num matches the module 01 modified equation exactly
  2. nu_num is LARGEST as dt -> 0 (shrinking dt makes diffusion worse)
  3. the table of effective Reynolds numbers for every case actually run
  4. the Ghia reference data is well formed and self-consistent

Run:
    python3 tests/test_numerical_viscosity.py
"""

import csv
import io
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.numerical_viscosity import (cell_reynolds, effective_reynolds,
                                     numerical_viscosity, report)

HERE = os.path.dirname(__file__)
REF = os.path.join(HERE, "..", "reference")


def test_matches_modified_equation():
    """nu_num = |c| dx (1-C)/2, and it vanishes at C = 1 exactly.

    Same expression module 01 verified against the Gaussian width-growth law
    sigma(t)^2 = sigma0^2 + 2 nu t. Reused here, not re-derived.
    """
    dx, u = 0.01, 1.0
    print("  " + "C".rjust(6) + "nu_num".rjust(12) + "   expected")
    for C in (0.0, 0.2, 0.5, 0.9, 1.0):
        dt = C * dx / u
        nn = numerical_viscosity(u, dx, dt)
        exp = u * dx * (1 - C) / 2
        print("  " + f"{C:6.2f}" + f"{nn:12.6f}" + f"   {exp:.6f}")
        assert abs(nn - exp) < 1e-15
    assert numerical_viscosity(u, dx, dx / u) == 0.0, "must vanish at C = 1"
    print("  matches the modified equation; zero at C = 1: PASS\n")


def test_smaller_dt_is_worse():
    """Counterintuitive but important: shrinking dt INCREASES numerical diffusion.

    The (1 - C) factor means nu_num -> |u| dx / 2 as dt -> 0. Taking a smaller
    timestep 'to be safe' buys stability margin and costs accuracy. Only
    refining dx reduces nu_num.
    """
    dx, u = 0.01, 1.0
    prev = None
    print("  " + "dt".rjust(10) + "C".rjust(8) + "nu_num".rjust(12))
    for dt in (5e-3, 2e-3, 1e-3, 1e-4):
        nn = numerical_viscosity(u, dx, dt)
        print("  " + f"{dt:10.0e}" + f"{u*dt/dx:8.2f}" + f"{nn:12.6f}")
        if prev is not None:
            assert nn > prev, "nu_num should grow as dt shrinks"
        prev = nn
    print("  numerical diffusion grows as the timestep shrinks: PASS\n")


def test_effective_reynolds_table():
    """Every case actually run, and what Reynolds number it really solved."""
    cases = [("Re=100", 0.01, 51), ("Re=100", 0.01, 101), ("Re=100", 0.01, 201),
             ("Re=400", 0.0025, 101), ("Re=1000", 0.001, 101),
             ("Re=1000", 0.001, 201)]

    print("  Near the lid, |u| = 1.0 (worst case):")
    rows = report(cases, u_mag=1.0)
    print()
    print("  Representative interior speed, |u| = 0.2:")
    report(cases, u_mag=0.2)

    by = {(l, N): r for l, _, N, _, _, r, _, _ in
          [(a, b, c, d, e, f, g, h) for a, b, c, d, e, f, g, h in rows]}
    assert by[("Re=100", 101)] < 0.5, "Re=100 on 101^2 should be mildly affected"
    assert by[("Re=1000", 101)] > 3.0, "Re=1000 on 101^2 should be dominated"
    assert by[("Re=1000", 201)] < by[("Re=1000", 101)], \
        "refining the grid must reduce the ratio"
    print("\n  numerical viscosity dominates at high Re, eases on finer grids: PASS\n")


def _load(path):
    rows = [r for r in open(path) if not r.startswith("#")]
    return list(csv.DictReader(io.StringIO("".join(rows))))


def test_reference_data_is_sane():
    """Structural checks on the transcribed Ghia tables.

    These cannot prove the digits are right -- only your own copy of the paper
    can do that. They do catch a transcription that broke the structure:
    wrong row count, a coordinate out of [0,1], a wall value that is not zero,
    or a sign flipped on a whole column.
    """
    t1 = _load(os.path.join(REF, "ghia1982_table1_u.csv"))
    t2 = _load(os.path.join(REF, "ghia1982_table2_v.csv"))

    print(f"  Table I  : {len(t1)} rows")
    print(f"  Table II : {len(t2)} rows")
    assert len(t1) == 17 and len(t2) == 17

    for name, tab, coord in (("I", t1, "y"), ("II", t2, "x")):
        cs = [float(r[coord]) for r in tab]
        assert cs == sorted(cs, reverse=True), f"Table {name} coord not descending"
        assert 0.0 <= min(cs) and max(cs) <= 1.0

    # Wall values: u = 1 at the lid, everything else zero at the walls.
    assert float(t1[0]["Re100"]) == 1.0, "u must equal the lid speed at y = 1"
    assert float(t1[-1]["Re100"]) == 0.0, "u must vanish at y = 0"
    assert float(t2[0]["Re100"]) == 0.0 and float(t2[-1]["Re100"]) == 0.0

    # u on the vertical centreline must change sign (recirculation).
    for re in ("Re100", "Re400", "Re1000"):
        col = np.array([float(r[re]) for r in t1])
        assert col.max() > 0 and col.min() < 0, f"{re}: no sign change in u"
        print(f"  Table I {re}: min {col.min():+.5f}  max {col.max():+.5f}")

    print("  reference tables are structurally sound: PASS")
    print("  (digits still need checking against your own copy)\n")


if __name__ == "__main__":
    tests = [
        ("nu_num matches the modified equation", test_matches_modified_equation),
        ("smaller dt is worse", test_smaller_dt_is_worse),
        ("effective Reynolds number table", test_effective_reynolds_table),
        ("Ghia reference data", test_reference_data_is_sane),
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
