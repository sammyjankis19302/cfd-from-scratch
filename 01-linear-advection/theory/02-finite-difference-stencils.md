# Where finite-difference stencils come from

**Assumed background:** derivatives, and the idea that a smooth function can be
written as a power series. Nothing else. Every step is shown.

Everything in this document comes from **one** tool — the Taylor series — used
in exactly **one** way. Once you see the pattern, you can derive any stencil
you will ever need, including ones that appear in no textbook.

---

## 0. The one tool

For a function `u(x)` that is smooth near `x_i`, the value a distance `Δx`
away is

$$u(x_i + \Delta x) = u(x_i) + \Delta x\,u'(x_i) + \frac{\Delta x^2}{2!}u''(x_i) + \frac{\Delta x^3}{3!}u'''(x_i) + \frac{\Delta x^4}{4!}u''''(x_i) + \dots$$

Write `u_i = u(x_i)`, `u_{i+1} = u(x_i + Δx)`, `u_{i-1} = u(x_i − Δx)`. The two
series we will use constantly are:

$$u_{i+1} = u_i + \Delta x\,u' + \frac{\Delta x^2}{2}u'' + \frac{\Delta x^3}{6}u''' + \frac{\Delta x^4}{24}u'''' + \frac{\Delta x^5}{120}u^{(5)} + \dots \tag{A}$$

$$u_{i-1} = u_i - \Delta x\,u' + \frac{\Delta x^2}{2}u'' - \frac{\Delta x^3}{6}u''' + \frac{\Delta x^4}{24}u'''' - \frac{\Delta x^5}{120}u^{(5)} + \dots \tag{B}$$

**(B) is (A) with `Δx → −Δx`.** That single fact is where all the sign
patterns come from: odd derivatives flip sign, even ones do not.

We will also need the two-step versions, obtained by putting `2Δx` in place of
`Δx`:

$$u_{i+2} = u_i + 2\Delta x\,u' + 2\Delta x^2 u'' + \frac{8\Delta x^3}{6}u''' + \frac{16\Delta x^4}{24}u'''' + \frac{32\Delta x^5}{120}u^{(5)} + \dots \tag{C}$$

$$u_{i-2} = u_i - 2\Delta x\,u' + 2\Delta x^2 u'' - \frac{8\Delta x^3}{6}u''' + \frac{16\Delta x^4}{24}u'''' - \frac{32\Delta x^5}{120}u^{(5)} + \dots \tag{D}$$

**The strategy, every single time:** take a combination of these series,
arrange for the terms you do not want to cancel, and divide by whatever is
left in front of the term you do want.

---

## 1. Forward difference, FD-1

Take (A) and subtract `u_i`:

$$u_{i+1} - u_i = \Delta x\,u' + \frac{\Delta x^2}{2}u'' + \frac{\Delta x^3}{6}u''' + \dots$$

Divide by `Δx`:

$$\frac{u_{i+1} - u_i}{\Delta x} = u' + \frac{\Delta x}{2}u'' + \frac{\Delta x^2}{6}u''' + \dots$$

Rearranging for the thing we actually want:

$$\boxed{u' = \frac{u_{i+1} - u_i}{\Delta x} - \frac{\Delta x}{2}u'' + O(\Delta x^2)}$$

The stencil is `(u_{i+1} − u_i)/Δx`. Everything after it is the error we make
by throwing those terms away — the **truncation error**. Its leading term is
`−(Δx/2)u''`, proportional to `Δx¹`, so this is **first-order accurate**.

*What first-order means concretely:* halve `Δx`, and the error halves.

---

## 2. Backward difference, BD-1

Take (B), and form `u_i − u_{i−1}`:

$$u_i - u_{i-1} = \Delta x\,u' - \frac{\Delta x^2}{2}u'' + \frac{\Delta x^3}{6}u''' - \dots$$

Divide by `Δx`:

$$\boxed{u' = \frac{u_i - u_{i-1}}{\Delta x} + \frac{\Delta x}{2}u'' + O(\Delta x^2)}$$

Also first order. Note the leading error is `+(Δx/2)u''`, the **opposite sign**
to FD-1 — same magnitude, opposite direction.

**This sign difference is not cosmetic.** Substituted into the advection
equation, one produces an effective `+ν u_xx` (physical-looking diffusion,
which damps and is stable) and the other `−ν u_xx` (anti-diffusion, which
amplifies every mode and blows up). That is the entire reason FTBS works for
`c > 0` and FTFS does not.

---

## 3. Central difference, CD-2

Now subtract (B) from (A). Watch what happens to each term:

$$u_{i+1} - u_{i-1} = \underbrace{(u_i - u_i)}_{0} + 2\Delta x\,u' + \underbrace{\left(\frac{\Delta x^2}{2} - \frac{\Delta x^2}{2}\right)}_{0}u'' + \frac{2\Delta x^3}{6}u''' + \underbrace{0}_{u''''} + \dots$$

**Every even-derivative term cancels**, because it has the same sign in both
series. That is the payoff of symmetry. What survives:

$$u_{i+1} - u_{i-1} = 2\Delta x\,u' + \frac{\Delta x^3}{3}u''' + \dots$$

Divide by `2Δx`:

$$\boxed{u' = \frac{u_{i+1} - u_{i-1}}{2\Delta x} - \frac{\Delta x^2}{6}u''' + O(\Delta x^4)}$$

**Second-order accurate**, using the same number of points as FD-1 or BD-1.
Symmetry bought an entire order of accuracy for free.

But look at what the leading error term now *is*: `u'''`, an **odd** derivative.
Odd-derivative errors are **dispersive** — they make different wavelengths
travel at different speeds. Even-derivative errors are **dissipative** — they
damp amplitude. So the two families of stencil fail in qualitatively different
ways:

| stencil | leading error | type | how it looks |
|---|---|---|---|
| FD-1 / BD-1 | `Δx u''` | dissipative | peaks flatten, corners round off |
| CD-2 | `Δx² u'''` | dispersive | wiggles appear, packets fall apart |

A one-sided stencil smears your solution. A central stencil keeps the
amplitude but scrambles the phase. Neither is simply "better".

---

## 4. Second-order backward, BD-2

Can we be upwind *and* second order? We need a third point on the upwind side.
Look for coefficients `a, b, d` such that

$$\frac{a\,u_i + b\,u_{i-1} + d\,u_{i-2}}{\Delta x} \approx u'$$

Substitute (B) and (D) and collect by derivative:

- **`u_i` terms:** `a + b + d` — must be `0` (we want a derivative, not a value)
- **`u'` terms:** `−b Δx − 2d Δx` — must equal `Δx` (this is the term we keep)
- **`u''` terms:** `b Δx²/2 + 2d Δx²` — must be `0` (kill it for second order)

Three equations:

$$a + b + d = 0, \qquad -b - 2d = 1, \qquad \frac{b}{2} + 2d = 0$$

From the third, `b = −4d`. Substituting into the second: `4d − 2d = 1`, so
`d = 1/2` and `b = −2`. From the first, `a = 3/2`. Therefore:

$$\boxed{u' = \frac{3u_i - 4u_{i-1} + u_{i-2}}{2\Delta x} + \frac{\Delta x^2}{3}u''' + O(\Delta x^3)}$$

(The leading error follows by keeping the `u'''` terms: `−4(−Δx³/6) + (−8Δx³/6)
= −4Δx³/6`, divided by `2Δx`, giving `−Δx²/3` on the stencil side and `+Δx²/3`
when moved across.)

**Second order, and entirely upwind.** This is the spatial part of the
Beam–Warming scheme.

**A warning that cost me a debugging session.** A second-order *spatial*
stencil marched with plain forward Euler in time is **unconditionally
unstable** — just like FTCS. Spatial accuracy buys nothing if the time
discretisation stays first order. Beam–Warming needs its `C²` correction term
to be second order in time as well. Order of accuracy in space and stability
are separate properties, and you must check both.

---

## 5. Fourth-order central, CD-4

CD-2 killed the `u''` term. To kill `u'''` as well we need four points. Take
the combination `α(u_{i+1} − u_{i−1}) + β(u_{i+2} − u_{i−2})`.

From §3 we already know `u_{i+1} − u_{i−1} = 2Δx u' + (Δx³/3)u''' + (Δx⁵/60)u⁽⁵⁾ + …`

By the same cancellation applied to (C) and (D):

$$u_{i+2} - u_{i-2} = 4\Delta x\,u' + \frac{16\Delta x^3}{6}u''' + \frac{64\Delta x^5}{60}u^{(5)} + \dots$$

Two conditions — keep `u'`, kill `u'''`:

$$2\alpha + 4\beta = 1 \qquad\text{(coefficient of } \Delta x\,u')$$
$$\frac{\alpha}{3} + \frac{16\beta}{6} = 0 \qquad\text{(coefficient of } \Delta x^3 u''')$$

The second gives `α = −8β`. Substituting: `−16β + 4β = 1`, so `β = −1/12` and
`α = 8/12 = 2/3`. Assembling:

$$\boxed{u' = \frac{-u_{i+2} + 8u_{i+1} - 8u_{i-1} + u_{i-2}}{12\Delta x} + \frac{\Delta x^4}{30}u^{(5)} + O(\Delta x^6)}$$

**Fourth-order accurate.** Halve `Δx` and the error drops by a factor of 16.

*Why this matters far more than it sounds.* To reach a fixed error level, a
fourth-order scheme needs dramatically fewer points per direction than a
second-order one. In 3D that saving is cubed — the cell-count advantage can be
a factor of roughly 30. This is the entire economic argument for high-order
methods in DNS and LES.

The price: a wider stencil, which is awkward at boundaries, and no guarantee of
better behaviour for non-smooth solutions — the Taylor series assumed
smoothness, and a shock has none.

---

## 6. Second derivative, CD-2

For the diffusion and Poisson equations we need `u''`. **Add** (A) and (B)
instead of subtracting:

$$u_{i+1} + u_{i-1} = 2u_i + \Delta x^2 u'' + \frac{2\Delta x^4}{24}u'''' + \dots$$

Now the **odd** terms cancel. Rearranging:

$$\boxed{u'' = \frac{u_{i+1} - 2u_i + u_{i-1}}{\Delta x^2} - \frac{\Delta x^2}{12}u'''' + O(\Delta x^4)}$$

Second-order accurate. This one stencil is the backbone of the Jacobi and
Gauss–Seidel solvers in module 05 and the whole cavity solver in module 06.

---

## 7. Summary

| stencil | formula | order | leading error | type |
|---|---|---|---|---|
| FD-1 | `(u_{i+1} − u_i)/Δx` | 1 | `−(Δx/2)u''` | dissipative |
| BD-1 | `(u_i − u_{i−1})/Δx` | 1 | `+(Δx/2)u''` | dissipative |
| CD-2 | `(u_{i+1} − u_{i−1})/2Δx` | 2 | `−(Δx²/6)u'''` | dispersive |
| BD-2 | `(3u_i − 4u_{i−1} + u_{i−2})/2Δx` | 2 | `+(Δx²/3)u'''` | dispersive |
| CD-4 | `(−u_{i+2} + 8u_{i+1} − 8u_{i−1} + u_{i−2})/12Δx` | 4 | `+(Δx⁴/30)u⁽⁵⁾` | dispersive |
| CD-2 (2nd deriv) | `(u_{i+1} − 2u_i + u_{i−1})/Δx²` | 2 | `−(Δx²/12)u''''` | — |

**Verified numerically** (`u = sin x` at `x = 0.7`, measured convergence rate):

```
stencil        h=0.1      h=0.05     h=0.025  measured order
FD-1       3.346e-02   1.642e-02   8.132e-03      1.01
BD-1       3.091e-02   1.578e-02   7.973e-03      0.99
CD-2       1.274e-03   3.186e-04   7.967e-05      2.00
BD-2       2.701e-03   6.569e-04   1.618e-04      2.02
CD-4       2.546e-06   1.593e-07   9.958e-09      4.00
```

The leading truncation-error *coefficients* were also checked against the
predictions above and agree to four significant figures. Derivations are
worthless until they have been tested against a machine that does not care
what you believe.

---

## 8. The general recipe

You now have everything needed to construct any stencil:

1. Decide which points you will use (which side, how many).
2. Write the Taylor series for each of them about `x_i`.
3. Take a general linear combination with unknown coefficients.
4. Set the coefficient of every derivative you *don't* want to zero, and the
   one you *do* want to one.
5. Solve the linear system. `n` points gives `n` equations.
6. Whatever term survives first is your truncation error, and its power of
   `Δx` is your order of accuracy.

**Two rules to carry forward, both derived above rather than asserted:**

- Symmetric stencils cancel one whole family of terms, so they get an extra
  order of accuracy free — but their leading error is odd-derivative, hence
  dispersive.
- One-sided stencils have even-derivative leading error, hence dissipative,
  and they encode a direction — which is exactly what you need when the
  physics has one.

**And one caution.** Order of accuracy is a statement about the limit
`Δx → 0`. It tells you how fast error falls as you refine. It does **not** tell
you which scheme is more accurate on the grid you can actually afford. Those
are different questions, and answering the second one is what global spectral
analysis in module 02 is for.
