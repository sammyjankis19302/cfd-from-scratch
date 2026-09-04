# 03 — 2D linear advection

$$\frac{\partial u}{\partial t} + c_x\frac{\partial u}{\partial x} + c_y\frac{\partial u}{\partial y} = 0$$

One dimension added. The exact solution is still known — `f(x − c_x t, y − c_y t)`
— so every deviation remains pure numerical error. What changes are two things
that have no 1D analogue, and both are invisible to a test with flow aligned to
an axis.

---

## Verification

Diagonal advection, `c_x = c_y = 1`, `C_x = C_y = 0.25`, on `sine_2d`.

| Scheme | N=40 | N=80 | N=160 | N=320 | measured `p` | formal |
|---|---|---|---|---|---|---|
| Upwind-2D | 1.581e-01 | 8.754e-02 | 4.617e-02 | 2.372e-02 | **0.96** | 1 |
| Lax–Wendroff-2D | 6.857e-03 | 1.713e-03 | 4.282e-04 | 1.071e-04 | **2.00** | 2 |

Plus: mass conserved to `2.0e-16`, and advection in `x` matches transposed
advection in `y` to `0.000e+00` — a check specifically for transposed indexing,
the classic 2D bug, which survives every symmetric test.

---

## The two results that matter

### The CFL condition is the sum, not each direction

The natural guess — `C_x ≤ 1` and `C_y ≤ 1` independently — is wrong. Von
Neumann analysis gives an amplification factor in which the two directions'
dissipation **adds**, so the worst mode (`θx = θy = π`, a checkerboard) requires

```
Cx + Cy <= 1
```

At `C_x = C_y` the limit per direction is **0.5, not 1**. In 3D it is 1/3 each.

| `Cx = Cy` | `Cx + Cy` | max\|u\| after 200 steps |
|---|---|---|
| 0.45 | 0.90 | 4.15e-20 |
| 0.49 | 0.98 | 2.85e-04 |
| 0.50 | **1.00** | **1.0000** — neutrally stable |
| 0.51 | 1.02 | 2.55e+03 |
| 0.55 | 1.10 | 6.86e+15 |

The test uses a checkerboard because a smooth Gaussian stayed bounded even at
`Cx + Cy = 1.02` — it contains almost no energy at that wavenumber. The
instability was real and simply invisible. **A stability test proves nothing
unless the initial condition excites the mode being tested.**

### Lax–Wendroff needs a cross-derivative term

In 2D, `u_tt` expands to `cx² u_xx + 2 cx cy u_xy + cy² u_yy`. The mixed term
is not optional — but it vanishes identically when `cy = 0`, so a scheme
missing it looks perfectly second order in any axis-aligned test.

| case | measured order |
|---|---|
| no cross term, flow along x | **2.00** |
| no cross term, flow diagonal | **1.00** |
| with cross term, flow diagonal | **2.00** |

A whole order of accuracy lost silently, in exactly the case an axis-aligned
test cannot detect. The broken version is kept in `src/schemes_2d.py` as
`lax_wendroff_2d_no_cross`, so this is demonstrated rather than asserted.

---

## Theory

[**01 — 2D linear advection**](theory/01-2d-advection.md): the characteristics
argument in 2D, the full Von Neumann derivation of `Cx + Cy ≤ 1`, the
cross-term expansion, and a write-up of a trap encountered in the verification
itself.

That last one is worth reading. The first convergence study gave Lax–Wendroff
orders of **1.97, 1.87, 1.63, 1.30** — getting *worse* as the grid refined. The
scheme was fine: the Gaussian used was `3.9e-03` at the domain edge, so it was
not periodic, and wrapping it created a kink that no scheme can represent — a
fixed-size error floor that stops shrinking with `Δx`. Module 01's theory notes
had already warned that the initial condition must genuinely be periodic. The
warning was written, then not followed.

**An order of convergence that degrades under refinement is almost always the
test, not the scheme.**

---

## Code

| File | Purpose |
|---|---|
| `src/schemes_2d.py` | Upwind-2D, Lax–Wendroff-2D, and the no-cross-term control |
| `src/initial_conditions_2d.py` | Grid, Gaussian, square patch, sine mode, exact solution |
| `src/solver_2d.py` | Time-marching driver and 2D L2 error |
| `tests/test_convergence_2d.py` | Five verification checks |

Arrays are indexed `u[j, i]` — row `j` is `y`, column `i` is `x`. Fixed once
and used everywhere; it is the main source of transposed-plot confusion in 2D.

The interface follows module 01: `step(u, Cx, Cy) → u_new`, knowing nothing
about grids or final times.

## Running

```bash
python3 tests/test_convergence_2d.py
```
