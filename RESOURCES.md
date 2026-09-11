# Resources

Everything I used while building this repository, in case it helps someone
learning the same material. Free and openly available sources are marked
**[free]**.

> **A note on what is here.** These are links and my own notes. No copyrighted
> PDFs are hosted in this repository — please support authors by obtaining
> books through proper channels. Where a free, author-sanctioned version exists,
> it is linked.

---

## Books

<!-- TODO Uday: keep only what you ACTUALLY used, and add the rest.
     A short honest list is worth more than a long aspirational one.
     For each, one line on what it was useful FOR. -->

| Book | Used for |
|---|---|
| *(fill in — keep only what you actually used)* | |

---

## Papers

Linked, not hosted. All are behind publishers; use your institutional access.

| Paper | Used for |
|---|---|
| Ghia, Ghia & Shin (1982), *J. Comput. Phys.* **48**, 387–411, [doi:10.1016/0021-9991(82)90058-4](https://doi.org/10.1016/0021-9991(82)90058-4) | The lid-driven cavity benchmark. Tables I and II are the centreline velocities module 06 validates against; Table IV gives wall vorticity. Their CSI-multigrid method is the contrast to this repository's Jacobi + explicit Euler. |
| Nallasamy & Krishna Prasad (1977), *J. Fluid Mech.* **79**(2), 391–414 | Cavity flow at high Reynolds number, Re up to 50 000. Secondary eddy behaviour. |
| Tokuda (1968), *J. Fluid Mech.* **33**(4), 657–672 | Impulsive motion of a flat plate — the unsteady startup problem. |
| Agarwal (2012), third-order upwind scheme for Navier–Stokes at high Re | Directly relevant to module 06: the cost of first-order upwind is the reason higher-order upwinding exists. |
| Jeong & Hussain, on the definition of a vortex (λ₂ criterion) | Vortex identification, for locating cavity vortex centres objectively. |

## Courses and lectures

<!-- TODO: NPTEL courses, YouTube lecture series, university notes. -->

| Source | Used for |
|---|---|
| *(fill in)* | |

---

## Software and tools

- **NumPy / Matplotlib** — everything numerical here
- **Python 3** on Ubuntu
- **Git** — version control; this repository was also how I learned it

---

## My own notes

<!-- TODO: if you scan your handwritten derivations, link them here.
     theory/ has the typed versions; the handwritten originals show the
     work was done by hand first. Scan as PDF, keep file sizes small. -->

---

## If you are starting from zero

Suggested order for someone with calculus and basic Python but no CFD:

1. **Understand why `u(x,t) = f(x − ct)`** before writing any code —
   [`01-linear-advection/theory/01-analytical-solution.md`](01-linear-advection/theory/01-analytical-solution.md).
   The method of characteristics also explains *why* upwinding works, which
   turns a rule you memorise into something obvious.
2. **Derive the stencils yourself** from Taylor series —
   [`01-linear-advection/theory/02-finite-difference-stencils.md`](01-linear-advection/theory/02-finite-difference-stencils.md).
   Everything else is built on these six formulas.
3. **Write the simplest possible solver** and get it wrong. Read
   [`01-linear-advection/BUGS.md`](01-linear-advection/BUGS.md) for the four
   mistakes I made and what each one cost.
4. **Verify before you trust.** Measure the order of accuracy. A code that runs
   is not a code that is correct.
5. **Learn what a convergence study cannot see** —
   [`02-spectral-analysis/`](02-spectral-analysis/). Order of accuracy describes
   the limit `Δx → 0`; it says nothing about the grid you can afford. This is
   where dispersion, dissipation and group velocity come from.
6. Then 2D, then nonlinear, then elliptic solvers, then Navier–Stokes.

Three things I learned the hard way, all documented in the modules:

- *The bugs you can find are limited by the diagnostic you choose.* A frozen
  boundary destroyed 10% of my domain and was invisible because I was only
  printing `max(u)`.
- *An order of convergence that gets worse under refinement is almost always
  the test, not the scheme.*
- *Grid-independence is how you tell a real result from a numerical artefact.*
- *A test that cannot fail is not a test.* Write the assertion, then ask what a
  broken version would print. If it prints the same thing, the test proves
  nothing.
- *Order of accuracy is not resolving power.* A first-order and a second-order
  stencil can need exactly the same points per wavelength.
