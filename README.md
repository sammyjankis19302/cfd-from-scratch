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
| 02 | [Spectral analysis](02-spectral-analysis/) — numerical wavenumber, dispersion, group velocity | **complete** |
| 03 | [2D advection](03-advection-2d/) — the 2D CFL condition, Lax–Wendroff cross term | **complete** |
| 04 | [Burgers](04-burgers/) — 1D and 2D nonlinear, conditional upwinding, shock formation | **complete** |
| 05 | [Poisson](05-poisson/) — Jacobi, Gauss–Seidel, SOR, manufactured solutions | **complete** |
| 06 | [Lid-driven cavity](06-lid-driven-cavity/) — vorticity–streamfunction, Ghia validation | in progress |

Modules are numbered by topic, not by the order they were built. 02 came last
of the completed modules: the spectral analysis is the most interesting part of
this repository, but 03–05 are what the lid-driven cavity depends on, so they
came first.

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

## Module 02 — resolving power is not order of accuracy

Feed `exp(ikx)` into a stencil and compare what comes back against the exact
derivative. The **numerical wavenumber** `k_num` is complex: its real part
carries dispersion, its imaginary part dissipation. All eight closed forms
verified against the stencils to `1e-14` or better.

**Points per wavelength for 1% phase error:**

| stencil | order | PPW |
|---|---|---|
| BD-1 | 1 | **25.6** |
| CD-2 | 2 | **25.6** |
| CD-4 | 4 | **8.3** |

BD-1 and CD-2 need *identical* resolution despite differing in formal order —
they share the real part `sin(K)/K`, and phase error lives entirely there. A
convergence study cannot distinguish them on this axis at all.

**Group velocity changes sign at 4 points per wavelength.** For CD-2,
`V_g/c = cos(K)`, so under-resolved modes carry energy *upstream* at up to full
speed — the spurious q-waves. Invisible to convergence studies (they vanish as
`Δx → 0`) and to stability analysis (`|G| ≤ 1` throughout).

---

## Module 06 — what first-order upwind costs

The cavity solver advects vorticity with first-order upwind, which module 01
showed adds a real diffusion `ν_num = |u|Δx(1−C)/2`. In the vorticity transport
equation that sits directly alongside the physical `ν`:

| case | grid | `ν` | `ν_num` | ratio | `Re_eff` |
|---|---|---|---|---|---|
| Re=100 | 101² | 0.0100 | 0.00450 | 0.45 | 69 |
| Re=400 | 101² | 0.0025 | 0.00450 | **1.80** | 143 |
| Re=1000 | 101² | 0.0010 | 0.00450 | **4.50** | 182 |

(near the lid; interior speeds give `Re_eff` of 91, 287, 505)

Since `ν_num` depends on `Δx` and not on `ν`, raising the target Reynolds number
at fixed grid eventually lets it dominate. This predicts — before any comparison
is run — that agreement with the Ghia benchmark should be good at Re = 100,
degrade sharply between 400 and 1000, and improve under grid refinement by *more*
at high Re than at low Re.

Corollary: `ν_num → |u|Δx/2` as `Δt → 0`, so **shrinking the timestep makes
numerical diffusion worse**. Only refining `Δx` reduces it.

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
