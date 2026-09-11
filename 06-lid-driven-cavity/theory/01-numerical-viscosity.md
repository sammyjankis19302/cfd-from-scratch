# What first-order upwind costs the cavity solver

The solver uses first-order upwind for vorticity advection. This document
quantifies what that means for the Reynolds number it actually solves — and it
does so **before** any comparison with Ghia, because the answer predicts the
comparison.

---

## 1. The argument

Module 01 derived the modified equation for FTBS. The scheme does not solve

$$u_t + c\,u_x = 0$$

It solves

$$u_t + c\,u_x = \nu_{\text{num}}\,u_{xx}, \qquad \nu_{\text{num}} = \frac{|c|\,\Delta x\,(1-C)}{2}$$

$\nu_{\text{num}}$ is not a metaphor. It has the units of viscosity and the
effect of viscosity. In the vorticity transport equation

$$\omega_t + u\,\omega_x + v\,\omega_y = \nu\,\nabla^2\omega$$

it sits **directly alongside the physical $\nu$**. The solver therefore runs at

$$\boxed{\;Re_{\text{eff}} \approx \frac{U L}{\nu + \nu_{\text{num}}}\;}$$

---

## 2. The consequence

$\nu_{\text{num}}$ depends on $\Delta x$, **not on $\nu$**. So as you target a
higher Reynolds number at fixed grid, $\nu$ falls and $\nu_{\text{num}}$ stays
where it is — until it dominates.

Near the lid, where $|u| \approx 1$, with $\Delta t = 10^{-3}$:

| case | grid | $\Delta x$ | $\nu$ | $\nu_{\text{num}}$ | ratio | $Re_{\text{eff}}$ | $Re_{\text{cell}}$ |
|---|---|---|---|---|---|---|---|
| Re=100 | $51^2$ | 0.0200 | 0.0100 | 0.00950 | 0.95 | 51 | 2.00 |
| Re=100 | $101^2$ | 0.0100 | 0.0100 | 0.00450 | 0.45 | 69 | 1.00 |
| Re=100 | $201^2$ | 0.0050 | 0.0100 | 0.00200 | 0.20 | 83 | 0.50 |
| Re=400 | $101^2$ | 0.0100 | 0.0025 | 0.00450 | **1.80** | 143 | 4.00 |
| Re=1000 | $101^2$ | 0.0100 | 0.0010 | 0.00450 | **4.50** | 182 | 10.00 |
| Re=1000 | $201^2$ | 0.0050 | 0.0010 | 0.00200 | **2.00** | 333 | 5.00 |

With a representative interior speed $|u| \approx 0.2$ the picture is milder —
$Re_{\text{eff}}$ of 91, 287 and 505 for the three $101^2$ cases — so the truth
sits between the two tables, worst near the lid and mildest in the core.

**At Re = 1000 on a $101^2$ grid the numerical viscosity is roughly 4.5 times
the physical viscosity near the lid.** The run is not solving Re = 1000.

---

## 3. What this predicts

Three falsifiable statements, all made before looking at any Ghia comparison:

1. **Re = 100 should agree well.** $\nu_{\text{num}}/\nu \le 0.45$, and less in
   the interior.
2. **Agreement should degrade monotonically with Re**, sharply between 400 and
   1000, because the ratio crosses 1 in that interval.
3. **Refining $101^2 \to 201^2$ at Re = 1000 should move the solution
   measurably toward Ghia**, because it halves $\nu_{\text{num}}$ — and the
   improvement should be *larger* than the same refinement produces at Re = 100,
   where $\nu_{\text{num}}$ was never the limiting factor.

Point 3 is the interesting one: it predicts that grid refinement does something
*qualitatively different* at high Re than at low Re. A generic "finer is better"
expectation does not predict that.

---

## 4. The cell Reynolds number

The standard shorthand for the same idea:

$$Re_{\text{cell}} = \frac{|u|\,\Delta x}{\nu}$$

Above roughly $Re_{\text{cell}} \approx 2$, the upwind scheme's own diffusion is
comparable to the physical diffusion it is meant to be resolving. The table
above gives $Re_{\text{cell}} = 10$ for Re = 1000 on $101^2$.

Reaching $Re_{\text{cell}} \le 2$ at Re = 1000 needs $\Delta x \le 0.002$, i.e.
a $501^2$ grid — about 25 times the cells of $101^2$, and a proportionally
smaller timestep. That is the real cost of first-order upwind at high Re, and it
is why Ghia used a different approach.

---

## 5. A counterintuitive corollary

$$\nu_{\text{num}} = \frac{|u|\,\Delta x\,(1-C)}{2}, \qquad C = \frac{|u|\Delta t}{\Delta x}$$

The $(1-C)$ factor means $\nu_{\text{num}} \to |u|\Delta x/2$ as
$\Delta t \to 0$. Measured at $\Delta x = 0.01$, $|u| = 1$:

| $\Delta t$ | $C$ | $\nu_{\text{num}}$ |
|---|---|---|
| 5e-03 | 0.50 | 0.002500 |
| 2e-03 | 0.20 | 0.004000 |
| 1e-03 | 0.10 | 0.004500 |
| 1e-04 | 0.01 | 0.004950 |

**Shrinking the timestep makes the numerical diffusion worse.** Taking a smaller
$\Delta t$ "to be safe" buys stability margin and pays for it in accuracy. Only
refining $\Delta x$ reduces $\nu_{\text{num}}$.

---

## 6. Connection to module 02

Global spectral analysis says the same thing in wavenumber space. BD-1's
imaginary part $-\,(1-\cos K)/K$ is exactly this dissipation, wavenumber by
wavenumber, and it reaches $-0.6366$ at $K = \pi$. Module 02 also shows that
BD-1 needs **25.6 points per wavelength** for 1% phase error — so a cavity
feature spanning fewer than ~26 cells is not faithfully represented regardless
of how converged the iteration is.

At Re = 1000 the boundary layer thickness scales roughly as $Re^{-1/2} \approx
0.03$, which is **3 cells** on a $101^2$ grid. Module 02 says that is far below
the resolution needed. The two analyses agree, from different directions.

---

## 7. Why this makes the comparison stronger, not weaker

A disagreement with Ghia at Re = 1000 is not a failed validation. It is a
**measured, predicted consequence of a documented property of the
discretisation** — and the prediction was made from the modified equation before
the comparison was run.

"My solver disagrees with the benchmark" is a weak result. "My solver's
first-order upwind adds $\nu_{\text{num}} = 4.5\nu$ at this grid and Reynolds
number, so it should behave like $Re \approx 180$–500, and here is the
measurement confirming that" is a different kind of statement entirely.
