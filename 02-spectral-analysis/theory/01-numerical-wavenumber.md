# The numerical wavenumber

Module 01 derived every stencil from Taylor series and measured its order of
accuracy. Order of accuracy is a statement about the limit $\Delta x \to 0$: it
says how fast error falls as you refine. **It does not say which scheme is more
accurate on the grid you can actually afford.**

That second question is the one that matters in practice, and this is how you
answer it.

---

## 1. The idea

Feed a single Fourier mode into a stencil and see what comes back.

The exact first derivative of $e^{\mathrm{i}kx}$ is $\mathrm{i}k\,e^{\mathrm{i}kx}$.
A stencil returns $\mathrm{i}k_{\text{num}}\,e^{\mathrm{i}kx}$, where
$k_{\text{num}}$ is the **numerical wavenumber** — in general complex, and in
general not equal to $k$.

Everything follows from comparing $k_{\text{num}}/k$ against 1:

| | meaning | consequence |
|---|---|---|
| $\mathrm{Re}(k_{\text{num}}/k)\ne 1$ | mode travels at the wrong speed | **dispersion** |
| $\mathrm{Im}(k_{\text{num}}/k)\ne 0$ | mode grows or decays | **dissipation** |

Module 01's third derivation established the baseline: for the exact advection
equation, $|G| = 1$ and $c_p = c$ for *every* wavenumber. This is the comparison
against it.

Write $K = k\,\Delta x$. Since the shortest representable wave spans two grid
points, $K \in (0,\pi]$, and points per wavelength is $\text{PPW} = 2\pi/K$.

---

## 2. Deriving one, then reading off the rest

Take BD-1. Apply it to $u_j = e^{\mathrm{i}kx_j}$:

$$\frac{u_j - u_{j-1}}{\Delta x} = \frac{e^{\mathrm{i}kx_j} - e^{\mathrm{i}k(x_j-\Delta x)}}{\Delta x} = \frac{1 - e^{-\mathrm{i}K}}{\Delta x}\,e^{\mathrm{i}kx_j}$$

Set this equal to $\mathrm{i}k_{\text{num}}e^{\mathrm{i}kx_j}$:

$$\mathrm{i}k_{\text{num}}\Delta x = 1 - e^{-\mathrm{i}K} = 1 - \cos K + \mathrm{i}\sin K$$

Divide by $\mathrm{i}$ and normalise by $K$:

$$\boxed{\;\frac{k_{\text{num}}}{k}\Big|_{\text{BD-1}} = \frac{\sin K}{K} - \mathrm{i}\,\frac{1-\cos K}{K}\;}$$

The same substitution for the others:

| stencil | $k_{\text{num}}/k$ |
|---|---|
| BD-1 | $\dfrac{\sin K}{K} - \mathrm{i}\dfrac{1-\cos K}{K}$ |
| FD-1 | $\dfrac{\sin K}{K} + \mathrm{i}\dfrac{1-\cos K}{K}$ |
| CD-2 | $\dfrac{\sin K}{K}$ |
| CD-4 | $\dfrac{8\sin K - \sin 2K}{6K}$ |

And for the second derivative, $k_{\text{num}}^2/k^2$:

| stencil | $k_{\text{num}}^2/k^2$ |
|---|---|
| CD-2 | $\left(\dfrac{\sin(K/2)}{K/2}\right)^2$ |
| CD-4 | $\dfrac{30 - 32\cos K + 2\cos 2K}{12K^2}$ |
| FD-2 | $\left(\dfrac{\sin(K/2)}{K/2}\right)^2(\cos K + \mathrm{i}\sin K)$ |
| BD-2 | $\left(\dfrac{\sin(K/2)}{K/2}\right)^2(\cos K - \mathrm{i}\sin K)$ |

**All eight verified** by applying the stencils numerically to $e^{\mathrm{i}kx}$
and comparing — agreement to $10^{-14}$ or better. Derivations are worthless
until tested against a machine that does not care what you believe.

---

## 3. What the algebra immediately tells you

**BD-1 and FD-1 have identical real parts and opposite imaginary parts.** Same
order, same phase behaviour, opposite stability. At $K=\pi$:

$$\mathrm{Im}(\text{BD-1}) = -0.6366, \qquad \mathrm{Im}(\text{FD-1}) = +0.6366$$

The entire difference between a scheme that works and one that explodes is that
sign. Module 01 showed this through the modified equation's $\pm\nu u_{xx}$;
here it is visible directly.

**Central stencils are purely real.** $\mathrm{Im} = 0$ at *every* wavenumber,
exactly. So CD-2 has no mechanism whatsoever to control the two-point wave — which
is why it cannot be used alone for advection, and why damping must come from the
time discretisation or from added artificial viscosity instead.

---

## 4. The result that changes how you pick a scheme

Ask: how many points per wavelength does each stencil need to keep the phase
error below 1%?

| stencil | order | PPW for 1% | PPW for 0.1% |
|---|---|---|---|
| BD-1 | 1 | **25.6** | 81.1 |
| FD-1 | 1 | **25.6** | 81.1 |
| CD-2 | 2 | **25.6** | 81.1 |
| CD-4 | 4 | **8.3** | 15.0 |

**BD-1 and CD-2 need exactly the same resolution** — because they share the real
part $\sin K/K$, and phase error lives entirely in the real part. One is first
order and the other second, and for *resolving power* they are indistinguishable.

Order of accuracy and resolving power are different properties. A convergence
study measures the first and says nothing about the second.

CD-4 needs about a third as many points per direction. In 3D that is roughly
$3^3 \approx 27\times$ fewer cells, which is the whole economic argument for
high-order methods in DNS and LES.

---

## 5. Group velocity, and waves that travel the wrong way

Phase speed says how fast a crest moves. **Group velocity** says how fast energy
moves — a wave packet's envelope:

$$\frac{V_g}{c} = \frac{\mathrm{d}\big(\mathrm{Re}(k_{\text{num}})\Delta x\big)}{\mathrm{d}K}$$

For the exact equation both equal $c$ at every $k$. For CD-2,
$k_{\text{num}}\Delta x = \sin K$, so

$$\boxed{\;\frac{V_g}{c}\bigg|_{\text{CD-2}} = \cos K\;}$$

| $K/\pi$ | CD-2 | CD-4 |
|---|---|---|
| 0.10 | 0.9511 | 0.9984 |
| 0.25 | 0.7071 | 0.9428 |
| 0.50 | **0.0000** | 0.3333 |
| 0.75 | **−0.7071** | −0.9428 |
| 1.00 | **−1.0000** | −1.6667 |

**$V_g$ changes sign at $K = \pi/2$, i.e. at 4 points per wavelength.** Below
that, energy travels *upstream* while the exact solution carries it downstream.
At $K=\pi$ it moves at full speed in the wrong direction.

These are the **spurious upstream q-waves**. They are not instability — the
amplitude can be perfectly bounded. They are not smearing. They are signal going
backwards, and:

- a convergence study will not see them (they vanish as $\Delta x \to 0$)
- a stability analysis will not see them ($|G| \le 1$ throughout)
- they appear as oscillations *ahead of* a feature, where nothing should be yet

CD-4's group velocity is better at low $K$ but *worse* at high $K$ — it reaches
$-1.67$, overshooting the exact speed in the wrong direction. Higher order is not
uniformly better; it is better where you resolve and worse where you do not.

---

## 6. Why this module exists

Everything in modules 01 and 03–05 was verified by refining the grid and watching
error fall. That is necessary and it is not sufficient. It cannot distinguish a
scheme that resolves a 6-point wave from one that needs 26 points, and it cannot
see energy propagating the wrong way at all.

Global spectral analysis answers the question a convergence study cannot:
**what does this scheme do to the waves actually present on the grid I can
afford?**
