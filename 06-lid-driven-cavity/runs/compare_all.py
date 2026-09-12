"""Compare every finished run against the Ghia benchmark.

Usage:
    python3 runs/compare_all.py           # summary table
    python3 runs/compare_all.py --full    # full point-by-point tables
"""

import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.ghia_compare import compare, print_table
from src.numerical_viscosity import numerical_viscosity

full = "--full" in sys.argv
files = sorted(glob.glob("runs/Re*_N*.npz"))
if not files:
    sys.exit("no runs found -- run runs/run_study.py first")

print(f"{'Re':>6}{'grid':>7}{'u L_inf':>10}{'u RMSE':>10}{'u rel%':>9}"
      f"{'v L_inf':>10}{'v RMSE':>10}{'v rel%':>9}{'psi_min':>11}{'steady':>9}")
print("-" * 91)

for f in files:
    r = dict(np.load(f))
    try:
        _, _, _, mu = compare(r, "u")
        _, _, _, mv = compare(r, "v")
    except KeyError as e:
        print(f"  skipping {f}: {e}")
        continue
    steady = "yes" if float(r["final_change"]) < 1e-7 else "NO"
    print(f"{int(r['Re']):6d}{int(r['N']):5d}^2"
          f"{mu['L_inf']:10.5f}{mu['RMSE']:10.5f}{mu['rel_pct']:9.2f}"
          f"{mv['L_inf']:10.5f}{mv['RMSE']:10.5f}{mv['rel_pct']:9.2f}"
          f"{float(r['psi'].min()):11.6f}{steady:>9}")
    if full:
        print()
        print_table(r, "u")
        print()
        print_table(r, "v")
        print()

print()
print("Any row marked NO under 'steady' mixes incomplete convergence with")
print("numerical viscosity. Increase T for that case before drawing conclusions.")
