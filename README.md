# CFD from Scratch

Finite-difference solvers for computational fluid dynamics, built from first
principles — from the 1D linear advection equation up to a lid-driven cavity
solver in vorticity–streamfunction form. No CFD libraries: every scheme, every
solver and every verification test is written from the discretisation up.

Each module keeps three things together that are usually kept apart:

- **theory/** — the derivation, done by hand, before any code exists
- **src/** and **tests/** — the implementation, and evidence that it is
  *correct* rather than merely runnable
- **experiments/** — numerical experiments in the form
  *question → prediction from theory → result → verdict*

The derivations are written to be self-contained. A reader needs derivatives
and power series; every algebraic step is shown.

---

## Progress

| # | Module | Status |
|---|---|---|
| 01 | [Linear advection](01-linear-advection/) — six schemes, Taylor-series stencils, verification | **complete** |
| 02 | Spectral analysis — numerical wavenumber, dispersion, group velocity, spurious modes | planned |
| 03 | [2D advection](03-advection-2d/) — the 2D CFL condition, Lax–Wendroff cross term | **complete** |
| 04 | Burgers — 1D and 2D nonlinear | planned |
| 05 | Poisson — Jacobi, Gauss–Seidel, SOR | planned |
| 06 | Lid-driven cavity — incompressible Navier–Stokes, Ghia validation | planned |

Modules are numbered by topic, not by the order they were built. 02 is
deliberately out of sequence: the spectral analysis is the most interesting part
of this repository, but modules 03–05 are what the lid-driven cavity actually
depends on, so they came first.

---

## Module 01 — verification summary

Gaussian advected once round a periodic domain at `C = 0.5`. Order `p` is
`log₂(e_coarse / e_fine)` between the two finest grids.

| Scheme | N=200 | N=1600 | measured `p` | formal |
|---|---|---|---|---|
| FTBS (upwind) | 3.952e-02 | 5.588e-03 | **0.97** | 1 |
| Lax–Friedrichs | 9.448e-02 | 1.615e-02 | **0.92** | 1 |
| Lax–Wendroff | 1.800e-03 | 2.815e-05 | **2.00** | 2 |
| BD-2 (Beam–Warming) | 1.800e-03 | 2.815e-05 | **2.00** | 2 |
| FTCS, FTFS | — | — | — | unconditionally unstable |

Also verified: upwind is **exact at `C = 1`** — max error `0.000e+00` on a
*square wave*, where any dissipation or dispersion would show immediately; mass
conserved to machine precision over five domain transits; and FTCS and FTFS do
diverge, reaching `10¹²` and `10¹¹²`, which confirms the test rig can actually
detect instability.

Finite-difference stencils (FD-1, BD-1, CD-2, BD-2, CD-4) are derived from
Taylor series and checked numerically — measured orders 1.01, 0.99, 2.00, 2.02,
4.00, with leading truncation coefficients matching theory to four significant
figures.

**One result worth singling out.** Lax–Wendroff and Beam–Warming produce
identical errors at exactly `C = 1/2`. Their leading dispersive coefficients
scale as `(1 − C²)` and `(1 − C)(2 − C)`, so the ratio is `(1+C)/(2−C)`, which
equals 1 only there. Measured `0.7648 / 1.0000 / 1.4999` against predicted
`0.7647 / 1.0000 / 1.5000`.

---

## Module 03 — verification summary

Diagonal advection on an exactly periodic initial condition.

| Scheme | N=40 | N=320 | measured `p` | formal |
|---|---|---|---|---|
| Upwind-2D | 1.581e-01 | 2.372e-02 | **0.96** | 1 |
| Lax–Wendroff-2D | 6.857e-03 | 1.071e-04 | **2.00** | 2 |

Two results with no 1D analogue, both invisible to a test with flow aligned to
an axis:

**The CFL limit is the sum.** Not `Cx ≤ 1` and `Cy ≤ 1` separately, but
`Cx + Cy ≤ 1` — so at `Cx = Cy` the limit per direction is 0.5, not 1. Measured
on a checkerboard (the worst Von Neumann mode): neutrally stable at exactly
`Cx + Cy = 1.00`, diverging at 1.02.

**Lax–Wendroff needs a cross-derivative term.** In 2D `u_tt` grows a
`2 cx cy u_xy` term, which vanishes when `cy = 0`. A scheme missing it measures
a flawless 2.00 along an axis and 1.00 on a diagonal. The broken version is kept
in the repository as a control.

---

## Running

```bash
pip install numpy matplotlib
cd 01-linear-advection
python3 advection_1d.py          # standalone upwind solver, explicit loops
python3 demo.py                  # all six schemes compared
python3 tests/test_convergence.py
```

---

## Author

**Akula Uday Kiran** — B.Tech Aerospace Engineering, RV College of
Engineering, Bengaluru. Interested in numerical methods for CFD: high-accuracy
schemes, compact discretisations, and spectral analysis of numerical error.
