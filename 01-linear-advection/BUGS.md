# Bugs in my first advection script

Four things were wrong with `first_attempt.py`. Keeping the record because
finding them taught me more than writing the original did.

---

### 1. The inflow boundary was frozen — this one was serious

```python
for i in range(1, nx - 1):        # u[0] is never updated by anything
```

For `c > 0`, information enters the domain at `x = 0`. Leaving `u[0]` stuck at
`sin(0) = 0` forever injected a permanent zero that marched inward at speed `c`
and destroyed the first unit of the domain.

At `T = 1`:

```
x=0.3   numerical +0.0010   exact -0.8090
x=0.5   numerical +0.0123   exact -1.0000     <- 100% wrong
x=0.9   numerical +0.2018   exact -0.3090
```

Mean absolute error was **0.675 for `x < 1`** against 0.246 elsewhere.

**Fix:** loop over every point, `range(nx)`, and give point 0 its left
neighbour from the end of the grid. Max error dropped from 1.01 to 0.385.

**What I actually learned:** I was only printing `max(u)` and `min(u)`. That
gave 0.641 against an exact 1.0 — which looks exactly like ordinary upwind
dissipation, so I moved on. The bug was in a region my diagnostic never looked
at. *The bugs you can find are limited by the diagnostic you choose.*

---

### 2. Duplicate grid point

`np.linspace(0, L, nx)` includes both `x = 0` and `x = L`. On a periodic
domain those are the same physical point, so one degree of freedom was stored
twice.

**Fix:** `nx` points, endpoint excluded, `dx = L/nx` instead of `L/(nx-1)`.

---

### 3. Signed error reported as maximum

```python
e = s - u
print("Error max is ", max(e))
```

A large *negative* error shows up as a small positive number, hiding exactly
the failures worth seeing.

**Fix:** `np.abs(s - u)`.

---

### 4. The exact solution written twice

`np.sin(np.pi*(x - c*T))` appeared in two places, one commented out.
Duplicated formulas drift apart. Now it lives in `exact()` and cannot disagree
with itself.

---

### Not a bug, but worth knowing

`CFL = 0.01` is the *worst* setting for accuracy. Numerical viscosity for FTBS
is `ν = c·Δx·(1−C)/2`, which is largest as `C → 0` and vanishes at `C = 1`. So
1000 timesteps bought maximum smearing. Compare `C = 0.8` and `C = 1.0`.

Also: `sin(πx)` has period 2, and `L = 10` fits exactly five wavelengths, so
periodic boundaries are legitimate here. With `L = 3` they would not be, and
the domain edge would quietly pollute every result.
