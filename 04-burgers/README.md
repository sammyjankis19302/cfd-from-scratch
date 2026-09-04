# 04 — Burgers equation

$$\text{1D:}\quad u_t + u\,u_x = \nu\,u_{xx} \qquad\qquad \text{2D:}\quad u_t + u\,u_x + v\,u_y = \nu(u_{xx}+u_{yy})$$

The first nonlinear equation in this repository. The advection speed is now the
solution itself, which changes three things that modules 01 and 03 could not
show.

---

## What is genuinely new

**1. The upwind direction is no longer known in advance.** With `c > 0` fixed,
a backward difference was always correct. Here the sign of the velocity varies
across the domain, so the direction must be chosen per grid point, per timestep:

```
u > 0  ->  information arrives from the LEFT   ->  backward difference
u < 0  ->  information arrives from the RIGHT  ->  forward difference
```

This is **conditional upwinding**, and it is not a stylistic choice.
`tests/` includes a negative control that uses a fixed backward direction
everywhere: on `sin(2πx)`, which is negative over half the domain, it reaches
`NaN` while the conditional version stays bounded at `0.841`. Where `u < 0`, a
backward difference is *downwind* — the anti-diffusive case from module 01.

In 2D, note which velocity selects which direction: the x-derivative of **both**
`u` and `v` is upwinded on the sign of `u`, because `u` is the advecting
velocity in `x`. Using `v` to upwind a `v_x` term is a natural-looking mistake
and is wrong.

**2. Shocks.** A perfectly smooth initial condition develops a discontinuity in
finite time. Characteristics carry `u` at speed `u`, so they cross at

$$t_{\text{break}} = \frac{1}{\max|u_0'(x)|}$$

and before that the maximum gradient has a closed form,
`max|u_x|(t) = 2π/(1 − 2πt)` for `u₀ = sin(2πx)`. Measured against it:

| `t` | measured `max\|u_x\|` | theory | error |
|---|---|---|---|
| 0.000 | 6.28 | 6.28 | 0.00% |
| 0.050 | 9.16 | 9.16 | 0.04% |
| 0.100 | 16.87 | 16.90 | 0.23% |
| 0.130 | 33.99 | 34.30 | 0.91% |

with `t_break = 1/(2π) = 0.1592`. A prediction that gives a *number at every
instant* and diverges at a *specific time* is far stronger evidence than "the
gradient goes up".

**3. Two independent stability limits.** The viscous term adds a constraint
unrelated to the CFL condition, and both must hold:

```
advective:  dt <= 1 / (|u|/dx + |v|/dy)
diffusive:  dt <= 0.5 / (nu (1/dx^2 + 1/dy^2))
```

| `ν` | `dt` advective | `dt` diffusive | binding |
|---|---|---|---|
| 1e-04 | 0.0500 | 25.0000 | advective |
| 1e-03 | 0.0500 | 2.5000 | advective |
| 1e-02 | 0.0500 | 0.2500 | advective |
| 1e-01 | 0.0500 | 0.0250 | **diffusive** |

Refining the grid tightens the diffusive limit as `dx²` but the advective limit
only as `dx`, so on fine grids the viscous term usually decides the timestep.

---

## The original version

`experiments/first_attempt_2d.py` is my hand-written 2D Burgers solver, kept
unchanged. It is **correct** — the conditional upwinding, the choice of which
velocity selects each direction, and the two-array update are all right, and the
Dirichlet boundary is legitimate here because `sin(πx)sin(πy)` is exactly zero
on all four edges of `[0,10]²`.

Two things about it are worth recording:

- It ran at `dt = 0.0001`, which is **CFL 0.002 — about 500× smaller than
  stability requires**. At `dt = 0.01` the answer changes by 0.5%
  (`max u = 0.5199` versus `0.5224`) for one hundredth of the work.
- The nested Python loops take **~10.6 minutes** for the full 10 000 steps. The
  vectorised version in `src/burgers.py` does it in 0.06 s and agrees to
  `1.2e-15` (`tests/test_equivalence_burgers.py`).

---

## A mistake worth reading

The shock test failed on its first run, which looked like a solver bug. It was
not — my predicted breaking time was wrong (`1/(4π²)` instead of `1/(2π)`).

What settled it: the measured gradient at `t = 0.1` was **16.46, 16.72, 16.82,
16.87** at `N = 200, 400, 800, 1600`. Essentially independent of the grid.
A numerical artefact changes when you change the grid; a converged physical
answer does not. **Grid-independence is how you tell a real result from a
numerical one**, and here it pointed straight at the theory rather than the code.

---

## Code

| File | Purpose |
|---|---|
| `experiments/first_attempt_2d.py` | The original hand-written solver, unchanged |
| `src/burgers.py` | Vectorised 1D and 2D, plus the stability-limit helper |
| `tests/test_equivalence_burgers.py` | Five checks |

## Running

```bash
python3 tests/test_equivalence_burgers.py
```
