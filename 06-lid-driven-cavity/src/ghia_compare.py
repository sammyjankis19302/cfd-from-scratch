"""Extract centreline profiles and compare them against the Ghia benchmark.

WHY INTERPOLATION IS NECESSARY

Ghia's tabulated locations are the grid points of a 129x129 mesh: 0.9766,
0.9688, 0.9609, 0.9531, ... On a 101x101 grid, 0.9766/0.01 = 97.66 -- there is
no such grid point. Taking the nearest point instead introduces an error of up
to dx/2 in POSITION, which near the lid (where du/dy is largest) is a large
error in VALUE, and it would be silently attributed to the scheme.

So the solution is interpolated to Ghia's coordinates. Linear interpolation is
used deliberately rather than a spline: the underlying discretisation is
second-order, and a higher-order interpolant would invent smoothness the
solution does not have.

The centrelines themselves, x = 0.5 and y = 0.5, ARE exact grid points on any
odd-N grid, so no interpolation is needed in that direction.

WHY MORE THAN ONE ERROR METRIC

Relative error is the natural first choice and it is treacherous here. Ghia's
Table I contains u = 0.00332 at y = 0.7344 -- a value that passes through zero.
An absolute discrepancy of 0.003 there is negligible physically and reads as
90% relative error. So report:

    L_inf   max absolute error       -- worst single point
    RMSE    root mean square error   -- overall agreement
    relative error, EXCLUDING points where |Ghia| is below a floor

and never quote relative error alone.
"""

import csv
import io
import os

import numpy as np

_REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "reference")


def load_ghia(table):
    """Load Table I ('u') or Table II ('v'). Returns (coords, {Re: values})."""
    fn = {"u": "ghia1982_table1_u.csv", "v": "ghia1982_table2_v.csv"}[table]
    rows = [r for r in open(os.path.join(_REF, fn)) if not r.startswith("#")]
    rd = list(csv.DictReader(io.StringIO("".join(rows))))
    coord_key = "y" if table == "u" else "x"
    coords = np.array([float(r[coord_key]) for r in rd])
    data = {int(k[2:]): np.array([float(r[k]) for r in rd])
            for k in rd[0] if k.startswith("Re")}
    return coords, data


def centreline_u(result):
    """u along the vertical centreline x = 0.5. Returns (y, u)."""
    N = int(result["N"])
    i_mid = (N - 1) // 2
    y = result["Y"][:, i_mid]
    return y, result["u"][:, i_mid]


def centreline_v(result):
    """v along the horizontal centreline y = 0.5. Returns (x, v)."""
    N = int(result["N"])
    j_mid = (N - 1) // 2
    x = result["X"][j_mid, :]
    return x, result["v"][j_mid, :]


def interpolate_to(coords_target, coords_have, values_have):
    """Linear interpolation onto Ghia's coordinates.

    np.interp needs increasing x, and Ghia's tables are listed descending.
    """
    order = np.argsort(coords_have)
    return np.interp(coords_target, coords_have[order], values_have[order])


def error_metrics(mine, reference, rel_floor=0.05):
    """L_inf, RMSE, and mean relative error over points above the floor.

    rel_floor excludes reference values too close to zero from the relative
    error, because a small absolute discrepancy there produces an enormous and
    meaningless percentage. The count of points used is returned so the reader
    can see how many were excluded.
    """
    d = mine - reference
    keep = np.abs(reference) > rel_floor
    rel = (np.abs(d[keep] / reference[keep]).mean() * 100.0
           if keep.any() else float("nan"))
    return {
        "L_inf": np.abs(d).max(),
        "RMSE": float(np.sqrt((d**2).mean())),
        "rel_pct": rel,
        "n_rel": int(keep.sum()),
        "n_total": len(reference),
    }


def compare(result, table="u", rel_floor=0.05):
    """Compare one finished run against the matching Ghia column."""
    Re = int(round(float(result["Re"])))
    coords, data = load_ghia(table)
    if Re not in data:
        raise KeyError(f"no Ghia column for Re={Re}; have {sorted(data)}")

    have_c, have_v = (centreline_u(result) if table == "u"
                      else centreline_v(result))
    mine = interpolate_to(coords, have_c, have_v)
    return coords, mine, data[Re], error_metrics(mine, data[Re], rel_floor)


def print_table(result, table="u"):
    coords, mine, ref, m = compare(result, table)
    label = "y" if table == "u" else "x"
    Re, N = int(round(float(result["Re"]))), int(result["N"])
    print(f"  {table}-velocity,  Re = {Re},  grid {N}x{N}")
    print(f"  {label:>7s}{'Ghia':>11s}{'mine':>11s}{'diff':>11s}")
    for c, a, b in zip(coords, ref, mine):
        print(f"  {c:7.4f}{a:11.5f}{b:11.5f}{b-a:+11.5f}")
    print(f"\n  L_inf = {m['L_inf']:.5f}   RMSE = {m['RMSE']:.5f}   "
          f"mean rel = {m['rel_pct']:.2f}% ({m['n_rel']}/{m['n_total']} pts)")
    return m
