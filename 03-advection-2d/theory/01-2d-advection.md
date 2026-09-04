# 2D linear advection: what changes and what doesn't

$$\frac{\partial u}{\partial t} + c_x\frac{\partial u}{\partial x} + c_y\frac{\partial u}{\partial y} = 0$$

Exact solution, by the same characteristics argument as in module 01 — now the
characteristics are straight lines in 3D $(x,y,t)$ space with direction
$(c_x, c_y, 1)$:

$$u(x,y,t) = f(x - c_x t,\; y - c_y t)$$

The profile translates diagonally, shape unchanged. Every deviation is again
pure numerical error.

Two things genuinely change in 2D. Both are easy to get wrong, and both are
invisible if you only test flow aligned with an axis.

---

## 1. The stability limit is the SUM

The natural guess is that the 1D condition applies per direction: $C_x \le 1$
and $C_y \le 1$. **That is wrong.**

Von Neumann analysis for 2D upwind. Substitute
$u_{j,i}^n = G^n e^{\mathrm{i}(k_x x_i + k_y y_j)}$ and write
$\theta_x = k_x\Delta x$, $\theta_y = k_y\Delta y$:

$$G = 1 - C_x\left(1 - e^{-\mathrm{i}\theta_x}\right) - C_y\left(1 - e^{-\mathrm{i}\theta_y}\right)$$

Each direction contributes its own dissipation, and the contributions **add**.
The worst case is the mode where both terms are maximally negative,
$\theta_x = \theta_y = \pi$, where $e^{-\mathrm{i}\pi} = -1$:

$$G = 1 - 2C_x - 2C_y$$

Requiring $|G| \le 1$ gives $-1 \le 1 - 2(C_x + C_y)$, hence

$$\boxed{\;C_x + C_y \le 1\;}$$

So for $C_x = C_y$ the limit per direction is **0.5, not 1**. Add a third
dimension and it becomes $C_x + C_y + C_z \le 1$, i.e. 1/3 each. Every
dimension you add tightens the timestep further.

**Verified numerically.** The worst mode $\theta_x = \theta_y = \pi$ is a
checkerboard, so that is what the test uses:

| $C_x = C_y$ | $C_x + C_y$ | $\max|u|$ after 200 steps |
|---|---|---|
| 0.45 | 0.90 | 4.15e-20 |
| 0.49 | 0.98 | 2.85e-04 |
| 0.50 | **1.00** | **1.0000** — neutrally stable |
| 0.51 | 1.02 | 2.55e+03 |
| 0.55 | 1.10 | 6.86e+15 |

The boundary sits exactly where the algebra says.

**Why the test uses a checkerboard.** A smooth Gaussian stayed bounded even at
$C_x + C_y = 1.02$, because it contains almost no energy at
$\theta_x = \theta_y = \pi$. The instability was real and simply invisible.
*The initial condition has to excite the mode you are testing*, or a stability
test proves nothing.

---

## 2. Lax–Wendroff grows a cross term

In 1D, Lax–Wendroff comes from a Taylor expansion in time with
$u_{tt} = c^2 u_{xx}$. In 2D:

$$u_t = -(c_x u_x + c_y u_y)$$

$$u_{tt} = c_x^2 u_{xx} + 2c_x c_y\,u_{xy} + c_y^2 u_{yy}$$

The mixed term $u_{xy}$ appears and is **not optional**. Substituting into
$u^{n+1} = u^n + \Delta t\,u_t + \tfrac{\Delta t^2}{2}u_{tt}$ and discretising
$u_{xy}$ on the four-corner stencil

$$u_{xy} \approx \frac{u_{j+1,i+1} - u_{j-1,i+1} - u_{j+1,i-1} + u_{j-1,i-1}}{4\Delta x\Delta y}$$

gives the extra contribution $\tfrac{1}{4}C_x C_y(u_{++} - u_{+-} - u_{-+} + u_{--})$.

**Why it is easy to miss.** If $c_y = 0$ the cross term vanishes identically.
So a scheme that omits it looks perfectly second order in any axis-aligned
test, and only fails when the flow has a diagonal component. Measured:

| case | order |
|---|---|
| no cross term, flow along $x$ | **2.00** |
| no cross term, flow diagonal | **1.00** |
| with cross term, flow diagonal | **2.00** |

An entire order of accuracy, lost silently, in the one case an axis-aligned
test cannot see. `src/schemes_2d.py` keeps the broken version as
`lax_wendroff_2d_no_cross` so this is demonstrated rather than asserted.

---

## A trap in the verification itself

The first convergence study here used a Gaussian with $\sigma = 0.15$ and gave
measured orders of **1.97, 1.87, 1.63, 1.30** for Lax–Wendroff — getting
*worse* as the grid refined.

The scheme was fine. That Gaussian is $3.9\times10^{-3}$ at the domain edge, so
it is **not periodic**. The periodic exact solution wraps it anyway, creating a
kink no scheme can represent — a fixed-size error that does not shrink with
$\Delta x$. As $\Delta x$ falls, that floor comes to dominate and the apparent
order collapses.

Switching to `sine_2d`, which is exactly periodic, gives a clean
**2.00, 2.00, 2.00**.

Module 01's theory notes already warned that the initial condition must
genuinely be periodic on the domain. The warning was written and then not
followed. *An order of convergence that gets worse under refinement is almost
always the test, not the scheme.*
