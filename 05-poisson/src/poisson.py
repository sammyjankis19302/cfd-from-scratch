"""Iterative solvers for the 2D Poisson equation on a uniform square grid.

    laplacian(s) = f          with Dirichlet boundaries

Setting f = 0 gives Laplace. For the lid-driven cavity the streamfunction
satisfies laplacian(psi) = -omega, so f = -omega there.

Discretised with the five-point stencil from module 01:

    (s[j,i+1] + s[j,i-1] + s[j+1,i] + s[j-1,i] - 4 s[j,i]) / dx^2  =  f[j,i]

Rearranged for s[j,i], that is the update every method below shares:

    s[j,i] = 0.25 * (s[j,i+1] + s[j,i-1] + s[j+1,i] + s[j-1,i] - dx^2 f[j,i])

The three methods differ only in WHICH values they read on the right.


TWO THINGS THAT MATTER MORE THAN THEY LOOK

1. CONVERGENCE MUST BE TESTED ON THE RESIDUAL, NOT ON THE CHANGE PER SWEEP.

   It is tempting to stop when max|s_new - s_old| is small. But for Jacobi the
   change and the residual are related exactly by

       change = residual * dx^2 / 4

   so the residual is 4/dx^2 times LARGER than the change you measured.
   Measured on Laplace with a lid, stopping at change < 1e-5:

       N=11  dx=0.1000   residual 3.9e-03    400x the change
       N=21  dx=0.0500   residual 1.6e-02   1600x
       N=41  dx=0.0250   residual 6.4e-02   6399x
       N=81  dx=0.0125   residual 2.6e-01  25583x

   Because Jacobi converges slowly, the change per sweep goes small long before
   the answer is right -- and the gap grows as dx^2, so it is worst exactly
   when you refine the grid. Every solver here therefore tests the residual.

2. ITERATION COUNT SCALES BADLY.

   Jacobi needs O(N^2) sweeps: measured 152, 502, 1564, 4465 for N = 11, 21,
   41, 81 -- about 4x per grid doubling. Gauss-Seidel halves that. SOR with the
   optimal relaxation factor reduces it to O(N), which is a different
   complexity class, not merely a constant factor. Inside a Navier-Stokes
   timestep loop this is the difference between a solver that finishes and one
   that does not.
"""

import numpy as np


def residual(s, f, dx):
    """r = laplacian(s) - f, evaluated on interior points only.

    This is the honest measure of how well the discrete equation is satisfied.
    Boundary points are held fixed by the Dirichlet condition, so they are
    excluded rather than being counted as perfectly converged.
    """
    r = np.zeros_like(s)
    r[1:-1, 1:-1] = (
        s[1:-1, 2:] + s[1:-1, :-2] + s[2:, 1:-1] + s[:-2, 1:-1]
        - 4.0 * s[1:-1, 1:-1]
    ) / dx**2 - f[1:-1, 1:-1]
    return r


def jacobi(s0, f, dx, tol=1e-6, max_iter=100000):
    """Jacobi: every point updated from the PREVIOUS sweep only.

    Needs two arrays, because no value from the current sweep may be used.
    Slowest of the three, and the easiest to reason about -- and the only one
    of the three that parallelises trivially, since no point depends on another
    point's new value.
    """
    s = s0.copy()
    for it in range(1, max_iter + 1):
        s_old = s.copy()
        s[1:-1, 1:-1] = 0.25 * (
            s_old[1:-1, 2:] + s_old[1:-1, :-2]
            + s_old[2:, 1:-1] + s_old[:-2, 1:-1]
            - dx**2 * f[1:-1, 1:-1]
        )
        if np.abs(residual(s, f, dx)).max() < tol:
            return s, it, True
    return s, max_iter, False


def gauss_seidel(s0, f, dx, tol=1e-6, max_iter=100000):
    """Gauss-Seidel: uses new values as soon as they are available.

    Sweeping in increasing j and i means s[j,i-1] and s[j-1,i] have ALREADY
    been updated this sweep, so the new information propagates immediately
    instead of waiting a full sweep. Roughly halves the iteration count for
    free, and needs only one array.

    The loop cannot be vectorised the way Jacobi can: the sequential dependency
    is the whole point of the method. Written with explicit loops for that
    reason.
    """
    s = s0.copy()
    ny, nx = s.shape
    for it in range(1, max_iter + 1):
        for j in range(1, ny - 1):
            for i in range(1, nx - 1):
                s[j, i] = 0.25 * (
                    s[j, i + 1] + s[j, i - 1] + s[j + 1, i] + s[j - 1, i]
                    - dx**2 * f[j, i]
                )
        if np.abs(residual(s, f, dx)).max() < tol:
            return s, it, True
    return s, max_iter, False


def sor(s0, f, dx, omega=None, tol=1e-6, max_iter=100000):
    """Successive over-relaxation: Gauss-Seidel, then overshoot deliberately.

        s_new = (1 - omega) * s_old + omega * s_gauss_seidel

    omega = 1 is exactly Gauss-Seidel. omega > 1 overshoots the correction,
    which sounds reckless and is in fact the single largest speed-up available
    here. Stability requires 0 < omega < 2.

    For Laplace on an N x N square with Dirichlet boundaries the optimal value
    is known analytically:

        omega_opt = 2 / (1 + sin(pi / (N - 1)))

    which tends to 2 as the grid refines. This is what turns O(N^2) iterations
    into O(N).
    """
    s = s0.copy()
    ny, nx = s.shape
    if omega is None:
        omega = optimal_omega(nx)

    for it in range(1, max_iter + 1):
        for j in range(1, ny - 1):
            for i in range(1, nx - 1):
                gs = 0.25 * (
                    s[j, i + 1] + s[j, i - 1] + s[j + 1, i] + s[j - 1, i]
                    - dx**2 * f[j, i]
                )
                s[j, i] = (1.0 - omega) * s[j, i] + omega * gs
        if np.abs(residual(s, f, dx)).max() < tol:
            return s, it, True
    return s, max_iter, False


def optimal_omega(N):
    """Analytically optimal SOR factor for an N x N Dirichlet square."""
    return 2.0 / (1.0 + np.sin(np.pi / (N - 1)))


SOLVERS = {"Jacobi": jacobi, "Gauss-Seidel": gauss_seidel, "SOR": sor}
