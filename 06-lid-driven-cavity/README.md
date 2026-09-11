# 06 — Lid-driven cavity (vorticity–streamfunction)

$$\omega_t + u\,\omega_x + v\,\omega_y = \nu\nabla^2\omega, \qquad \nabla^2\psi = -\omega, \qquad u = \psi_y,\;\; v = -\psi_x$$

Everything in modules 01–05 meets here. Poisson solve for `ψ` (module 05),
velocities by central differences (module 01), wall vorticity from ghost points,
vorticity advected by conditional upwinding (module 04) with central diffusion,
marched explicitly.

**Status: solver complete and running; Ghia comparison in progress.**

---

## The central result, stated before any comparison

The solver uses first-order upwind for vorticity advection. Module 01's modified
equation says that scheme adds a real diffusion:

```
nu_num = |u| dx (1 - C) / 2
```

In the vorticity transport equation this sits **directly alongside the physical
`ν`**. So the solver runs at `Re_eff ≈ UL/(ν + ν_num)`, not at `UL/ν`.

Near the lid (`|u| ≈ 1`), with `dt = 1e-3`:

| case | grid | `ν` | `ν_num` | ratio | `Re_eff` | `Re_cell` |
|---|---|---|---|---|---|---|
| Re=100 | 51² | 0.0100 | 0.00950 | 0.95 | 51 | 2.00 |
| Re=100 | 101² | 0.0100 | 0.00450 | 0.45 | 69 | 1.00 |
| Re=100 | 201² | 0.0100 | 0.00200 | 0.20 | 83 | 0.50 |
| Re=400 | 101² | 0.0025 | 0.00450 | **1.80** | 143 | 4.00 |
| Re=1000 | 101² | 0.0010 | 0.00450 | **4.50** | 182 | 10.00 |
| Re=1000 | 201² | 0.0010 | 0.00200 | **2.00** | 333 | 5.00 |

With an interior speed `|u| ≈ 0.2` the same cases give `Re_eff` of 91, 287 and
505 — so the truth lies between, worst near the lid.

**`ν_num` depends on `Δx`, not on `ν`.** Raise the target Re at fixed grid and
`ν` falls while `ν_num` stays put, until it dominates. At Re = 1000 on 101² the
numerical viscosity is 4.5× the physical.

### Three predictions this makes

1. Re = 100 should agree well with Ghia.
2. Agreement should degrade monotonically with Re, sharply between 400 and 1000,
   where the ratio crosses 1.
3. Refining 101² → 201² at Re = 1000 should move the solution measurably toward
   Ghia — and by **more** than the same refinement achieves at Re = 100, where
   `ν_num` was never the limiting factor.

Prediction 3 is the interesting one: grid refinement should do something
*qualitatively different* at high Re. "Finer is better" does not predict that.

### A corollary that catches people out

`ν_num → |u|Δx/2` as `Δt → 0`. Measured at `Δx = 0.01`, `|u| = 1`:

| `Δt` | `C` | `ν_num` |
|---|---|---|
| 5e-03 | 0.50 | 0.002500 |
| 1e-03 | 0.10 | 0.004500 |
| 1e-04 | 0.01 | 0.004950 |

**Shrinking the timestep makes numerical diffusion worse.** Only refining `Δx`
reduces it.

---

## Reference data

`reference/ghia1982_table1_u.csv` and `ghia1982_table2_v.csv` hold the Re = 100,
400 and 1000 columns of Ghia, Ghia & Shin (1982) Tables I and II — `u` on the
vertical centreline and `v` on the horizontal centreline.

Transcribed from the paper by OCR, then cleaned: the scan renders minus signs as
`a.`, `AI.`, `XJ.` and splits digits (`0.7887 1`). Every value was
cross-checked. **Verify against your own copy before publishing error metrics** —
`tests/` checks structure (row counts, coordinate ordering, wall values, sign
changes) but cannot prove individual digits.

Note: the paper's printed heading on Table II reads "u-Velocity"; the quantity
tabulated is `v`.

---

## Two fixes needed before extracting centrelines

Both present in `experiments/first_attempt_cavity.py`:

**`u` and `v` are never set on the boundaries.** Both loops run `range(1,n-1)`,
so `u[-1,:]` stays 0 where it should equal the lid speed 1. Harmless for the
interior march, but Ghia's Table I lists `y = 1.0000, u = 1.0000` — that row
would compare against zero. It is also why the streamplot shows a stationary lid.

**The Poisson stopping test measures the change per sweep, not the residual.**
Module 05 established the two differ by exactly `4/Δx²`; at `nx = 101` that is
**40,000×**. A tolerance of `1e-6` on the change corresponds to a residual near
`0.04`. Worth measuring the true residual once before publishing error tables.

---

## Theory

[**01 — What first-order upwind costs**](theory/01-numerical-viscosity.md) —
the full argument, the `Re_cell` analysis, the cost of reaching `Re_cell ≤ 2`
at Re = 1000 (a 501² grid, ~25× the cells), and the connection to module 02.

Module 02 reaches the same conclusion from wavenumber space: BD-1 needs 25.6
points per wavelength for 1% phase error, while the Re = 1000 boundary layer
(`~Re^-1/2 ≈ 0.03`) spans about **3 cells** on a 101² grid. Two analyses, one
answer.

---

## Code

| File | Purpose |
|---|---|
| `experiments/first_attempt_cavity.py` | The original hand-written solver, unchanged |
| `src/numerical_viscosity.py` | `ν_num`, `Re_cell`, `Re_eff`, and the case table |
| `reference/` | Ghia 1982 Tables I and II, Re = 100 / 400 / 1000 |
| `tests/test_numerical_viscosity.py` | Four checks |

The original solver is correct: Jacobi Poisson, ghost-point wall vorticity on all
four walls, conditional upwinding, central diffusion, explicit Euler, the
`v = −ψ_x` sign, corner exclusion, and the `o_old`/`o_new` discipline.

Runs completed: grids 51/101/201, Re 100/400/1000, `t* = 50/100/200`.

## Running

```bash
python3 tests/test_numerical_viscosity.py
```

---

## Why the comparison is stronger for this

A disagreement with Ghia at Re = 1000 is not a failed validation. It is a
measured consequence of a documented property of the discretisation, predicted
from the modified equation before the comparison was run.

"My solver disagrees with the benchmark" is weak. "First-order upwind adds
`ν_num = 4.5ν` at this grid and Reynolds number, so it should behave like
`Re ≈ 180–500`, and here is the measurement confirming it" is a different kind
of statement.
