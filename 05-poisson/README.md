# 05 — Poisson and Laplace

$$\nabla^2 s = f \qquad\text{with Dirichlet boundaries}$$

`f = 0` gives Laplace. For the lid-driven cavity the streamfunction satisfies
`∇²ψ = −ω`, so this module is the engine that runs inside every timestep of
module 06.

Three iterative methods on the five-point stencil derived in module 01, all
sharing one update and differing only in **which values they read on the right**.

---

## Verification

Method of manufactured solutions: pick `s = sin(πx)sin(πy)`, substitute, and
whatever falls out (`f = −2π²s`) is the source term that makes it true. The
solver must then reproduce the function you started from.

| Solver | iterations (N=41) | `max\|s − exact\|` |
|---|---|---|
| Jacobi | 6933 | 5.142e-04 |
| Gauss–Seidel | 3468 | 5.142e-04 |
| SOR | 172 | 5.142e-04 |

Same answer from all three, agreeing to `4.5e-12` — three different iterations,
one linear system.

**Discretisation is second order**, exactly as the stencil derivation predicts:

| N | max error | ratio | order |
|---|---|---|---|
| 11 | 8.265e-03 | | |
| 21 | 2.059e-03 | 4.01 | **2.01** |
| 41 | 5.142e-04 | 4.00 | **2.00** |
| 81 | 1.285e-04 | 4.00 | **2.00** |

This is separate from iterative convergence. Even a perfectly converged
iteration leaves the stencil's own `−dx²/12 · u''''` error. Iterating harder
cannot beat it; only a finer grid can.

---

## Iteration count is a complexity question, not a constant factor

| N | `ω_opt` | Jacobi | Gauss–Seidel | SOR |
|---|---|---|---|---|
| 11 | 1.5279 | 335 | 169 | 36 |
| 21 | 1.7295 | 1356 | 680 | 72 |
| 41 | 1.8545 | 5441 | 2722 | 148 |

**Scaling exponent per grid doubling: Jacobi 2.00, SOR 1.04.** Jacobi is
`O(N²)`, SOR is `O(N)`. Gauss–Seidel is reliably half of Jacobi — new values
are used the moment they exist rather than waiting a full sweep.

SOR is Gauss–Seidel that deliberately overshoots:

```
s_new = (1 - omega) * s_old + omega * s_gauss_seidel
```

`ω = 1` is exactly Gauss–Seidel; `ω > 1` overshoots the correction. It sounds
reckless and it is the single largest speed-up available here. For a square
Dirichlet domain the optimum is known analytically,
`ω_opt = 2/(1 + sin(π/(N−1)))`, tending to 2 as the grid refines.

Inside a Navier–Stokes timestep loop this is the difference between a solver
that finishes and one that does not.

---

## The stopping test — the thing that would have bitten me

`experiments/first_attempt_laplace.py` stops when `max|s_new − s_old| < tol`:
how much the solution *changed* this sweep. That is not how well it *satisfies
the equation*. For Jacobi the two are related exactly:

```
change = residual * dx^2 / 4
```

so the residual is `4/dx²` times larger than the quantity being tested:

| N | dx | change | residual | ratio | `4/dx²` |
|---|---|---|---|---|---|
| 11 | 0.1000 | 9.74e-06 | 3.894e-03 | 400.0 | 400.0 |
| 21 | 0.0500 | 9.96e-06 | 1.593e-02 | 1599.9 | 1600.0 |
| 41 | 0.0250 | 1.00e-05 | 6.397e-02 | 6398.6 | 6400.0 |
| 81 | 0.0125 | 9.99e-06 | 2.556e-01 | 25583.4 | 25600.0 |

Because Jacobi converges slowly, the change per sweep goes small long before
the answer is right. **And the gap grows as `dx²`, so it is worst exactly when
you refine the grid to get a better answer.** At N=81 a solution reported as
converged to `1e-5` has a residual of `0.26`.

This is the kind of defect that produces a Ghia comparison that is subtly wrong
with no visible cause. Every solver in `src/poisson.py` tests the residual.

---

## The original version

`experiments/first_attempt_laplace.py`, unchanged. The Jacobi update, the
Dirichlet boundaries and the two-array pattern are all correct. Besides the
stopping test above, one reporting bug: `print("Ittt =", max_iterations + 1)`
always prints 501 regardless of what happened. It converged at iteration
**152**. Should be `ti + 1` — which matters, because iteration count is the
entire point of comparing these three methods.

---

## Code

| File | Purpose |
|---|---|
| `experiments/first_attempt_laplace.py` | The original hand-written Jacobi solver |
| `src/poisson.py` | Jacobi, Gauss–Seidel, SOR, residual, optimal ω |
| `tests/test_poisson.py` | Five verification checks |

## Running

```bash
python3 tests/test_poisson.py
```

## Next

Module 06 couples this to vorticity transport: solve `∇²ψ = −ω` for the
streamfunction, differentiate to get velocities, advect vorticity, repeat.
Having verified the Poisson solve in isolation means that when the cavity
misbehaves, this piece is already ruled out.
