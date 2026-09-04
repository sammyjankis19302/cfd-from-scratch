# Why the solution is `u(x,t) = f(x − ct)`

**The equation**

$$\frac{\partial u}{\partial t} + c\,\frac{\partial u}{\partial x} = 0, \qquad u(x,0) = f(x), \qquad c = \text{constant} > 0$$

Everything in this repository is verified against the exact solution of this
equation, so the exact solution had better be right. Below are three
independent derivations of it. They are all worth seeing, because each one
reveals something different: the first tells you *why*, the second is the
fastest, and the third is the one that leads directly into spectral analysis.

---

## Derivation 1 — Method of characteristics

**The idea.** Instead of asking "what is `u` at a fixed point in space", ask
"is there a path through the `(x,t)` plane along which `u` does not change?"
If such paths exist, the PDE collapses to an ODE along them.

Suppose we move along some path `x(t)`. The value of `u` seen by an observer
riding that path is `u(x(t), t)`. How fast does it change? By the chain rule:

$$\frac{d}{dt}\,u\big(x(t),t\big) = \frac{\partial u}{\partial t} + \frac{dx}{dt}\,\frac{\partial u}{\partial x}$$

Now compare this with the PDE itself:

$$\frac{\partial u}{\partial t} + c\,\frac{\partial u}{\partial x} = 0$$

The two expressions are identical **provided we choose the path such that**

$$\boxed{\;\frac{dx}{dt} = c\;}$$

Make that choice, and the chain-rule expression becomes exactly the left-hand
side of the PDE, which is zero:

$$\frac{du}{dt} = 0 \qquad \text{along } \frac{dx}{dt} = c$$

**This is the whole content of the equation.** `u` is constant along paths
that move at speed `c`. Those paths are the *characteristics*.

Solve the path equation — it is a trivial ODE:

$$x(t) = x_0 + ct \quad\Longrightarrow\quad x_0 = x - ct$$

So the characteristic passing through the point `(x,t)` started at
`x_0 = x − ct` when `t = 0`. And since `u` never changed along that path, the
value at `(x,t)` is the value it had at the start:

$$u(x,t) = u(x_0, 0) = f(x_0) = \boxed{f(x - ct)}$$

**Read what this says.** The initial profile is translated rightward at speed
`c`, with its shape completely unaltered. No spreading, no decay, no change in
amplitude — forever. Any spreading or decay you see in a numerical solution is
therefore *entirely* an artefact of the discretisation. That is what makes
this equation the ideal test problem: every error is your error.

**Why the characteristics also tell you which stencil to use.** Information
arrives at `x` from the left, from `x − ct`. A backward difference
`(u_i − u_{i−1})/Δx` looks left, in the direction the information actually
came from — that is *upwinding*, and it is why FTBS is stable for `c > 0`
while FTFS is not. Stencil direction is a physical statement, not a taste.

---

## Derivation 2 — Change of variables

Same result, three lines, less insight. Define

$$\xi = x - ct, \qquad \eta = t$$

By the chain rule, with `u(x,t) = U(ξ,η)`:

$$\frac{\partial u}{\partial t} = \frac{\partial U}{\partial \xi}\frac{\partial \xi}{\partial t} + \frac{\partial U}{\partial \eta}\frac{\partial \eta}{\partial t} = -c\,U_\xi + U_\eta$$

$$\frac{\partial u}{\partial x} = \frac{\partial U}{\partial \xi}\frac{\partial \xi}{\partial x} + \frac{\partial U}{\partial \eta}\frac{\partial \eta}{\partial x} = U_\xi$$

Substitute into the PDE:

$$(-c\,U_\xi + U_\eta) + c\,U_\xi = 0 \quad\Longrightarrow\quad U_\eta = 0$$

`U` does not depend on `η` at all. So `U` is a function of `ξ` alone:

$$u(x,t) = U(\xi) = F(x - ct)$$

and applying `u(x,0) = f(x)` gives `F = f`. Same answer.

---

## Derivation 3 — Separation of variables and Fourier modes

This is the route that leads into the spectral analysis in module 02, so it
is worth doing carefully.

**Step 1 — assume a separable solution.** Try

$$u(x,t) = X(x)\,T(t)$$

Substituting into the PDE:

$$X T' + c\,X' T = 0$$

Divide through by `XT`:

$$\frac{T'}{T} = -c\,\frac{X'}{X}$$

The left side depends only on `t`; the right side only on `x`. Two functions
of independent variables can only be equal everywhere if both equal the same
constant. Call it `−ickc` (the factor is chosen to make the algebra come out
cleanly; any constant works and you would rediscover this one):

$$\frac{T'}{T} = -ick, \qquad \frac{X'}{X} = ik$$

**Step 2 — solve the two ODEs.** Both are first order and linear:

$$X(x) = A\,e^{ikx}, \qquad T(t) = B\,e^{-ickt}$$

**Step 3 — recombine.**

$$u_k(x,t) = C\,e^{ikx}\,e^{-ickt} = C\,e^{ik(x - ct)}$$

So every single Fourier mode is transported at speed `c`, unchanged in
amplitude. Note carefully what that means:

$$|e^{ik(x-ct)}| = 1 \quad \text{for every } k, \text{ for all time}$$

**No amplitude decay, and every wavenumber travels at exactly the same speed
`c`.** The exact equation is neither dissipative nor dispersive.

**Step 4 — superpose.** The equation is linear, so any sum of solutions is a
solution. Build the initial condition out of Fourier modes:

$$f(x) = \sum_k \hat{f}_k\,e^{ikx}$$

Advect each one:

$$u(x,t) = \sum_k \hat{f}_k\,e^{ik(x-ct)} = f(x-ct)$$

Same answer again, now via an infinite superposition.

---

## Why derivation 3 matters more than the other two

Derivation 3 sets up the entire framework of numerical error analysis.

The exact solution has, for every wavenumber `k`:

- **amplification factor** `|G| = 1` — no growth, no decay
- **phase speed** `c_p = c` — the same for all `k`
- **group velocity** `c_g = c` — the same for all `k`

A numerical scheme applied to the same problem gives some `G(kΔx, C)` which is
**not** equal to 1, and an effective wavenumber `k_num ≠ k`. The two failures
have names and separate physical consequences:

| | exact | numerical | consequence |
|---|---|---|---|
| `\|G\|` | 1 | `≠ 1` | **dissipation** — amplitude decays (or blows up) |
| `c_p(k)` | `c` for all `k` | depends on `k` | **dispersion** — modes separate, wave packet falls apart |

Numerical **dissipation** is why your sine wave lost 38% of its amplitude with
FTBS at `C = 0.01`. Numerical **dispersion** is why Lax–Wendroff grows
oscillations behind a square wave. Both are measured against the baseline this
derivation established: `|G| = 1`, `c_p = c`, for every `k`.

That comparison, done properly, is global spectral analysis — module 02.

---

## The periodic domain

All the code here uses periodic boundaries: `u(0,t) = u(L,t)`. Two reasons.

1. **It is exact for this problem.** A profile leaving the right edge re-enters
   at the left, which is precisely what `f(x − ct)` does if `f` is periodic
   with period `L`. So the boundary contributes *zero* error, and everything
   left over is scheme error. That is what makes a clean convergence study
   possible.
2. **It permits arbitrarily long runs.** You can advect a hundred times round
   the domain and watch error accumulate, without ever needing an inflow
   condition.

The corresponding exact solution wraps:

$$u(x,t) = f\big(\bmod(x - ct,\; L)\big)$$

which is exactly `np.mod(x - c*t, L)` in the code.

**The initial condition must genuinely be periodic on `[0,L)`.** In the code
here, `f(x) = sin(πx)` has period 2, and `L = 10` fits exactly five
wavelengths — so it is legitimate. `sin(πx)` on `L = 3` would *not* be, and
the domain edge would introduce a discontinuity that pollutes everything.
Check this before trusting any result.
