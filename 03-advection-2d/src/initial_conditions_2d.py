"""Grid and initial conditions for the 2D linear advection equation.

    du/dt + cx du/dx + cy du/dy = 0

The exact solution is u(x, y, t) = f(x - cx t, y - cy t): the initial profile is
translated diagonally, unchanged in shape. Same structure as 1D, one dimension
added, so every error is still pure numerical error.

Arrays are indexed u[j, i] -- ROW index j is y, COLUMN index i is x. This is the
NumPy convention (first index is the slowest-varying) and it is the single most
common source of transposed-plot confusion in 2D CFD. Fixed here once, used
everywhere.
"""

import numpy as np


def make_grid_2d(N, L=1.0):
    """Uniform periodic N x N grid on [0, L) x [0, L).

    Endpoints excluded for the same reason as in 1D: on a periodic domain
    x = L is the same point as x = 0.

    Returns
    -------
    X, Y : ndarray, shape (N, N) -- meshgrid coordinates, indexed [j, i]
    dx   : float (dy = dx here; the grid is square)
    """
    dx = L / N
    x = np.linspace(0.0, L - dx, N)
    X, Y = np.meshgrid(x, x, indexing="xy")   # X[j,i] = x[i], Y[j,i] = x[j]
    return X, Y, dx


def gaussian_2d(X, Y, L, x0=0.5, y0=0.5, sigma=0.1):
    """Circular Gaussian bump centred at (x0*L, y0*L)."""
    return np.exp(-((X - x0 * L) ** 2 + (Y - y0 * L) ** 2) / (2.0 * sigma**2))


def square_2d(X, Y, L, lo=0.2, hi=0.4):
    """Square patch equal to 1 on [lo*L, hi*L)^2.

    Four corners and four edges: exposes whether a scheme treats the two
    directions differently, which a smooth blob would hide.
    """
    return np.where(
        (X >= lo * L) & (X < hi * L) & (Y >= lo * L) & (Y < hi * L), 1.0, 0.0
    )


def sine_2d(X, Y, L, kx=1, ky=1):
    """Single 2D Fourier mode, kx wavelengths in x and ky in y."""
    return np.sin(2.0 * np.pi * kx * X / L) * np.sin(2.0 * np.pi * ky * Y / L)


def exact_solution_2d(f, X, Y, L, cx, cy, t, **kwargs):
    """u(x, y, t) = f(x - cx t, y - cy t), wrapped periodically."""
    return f(np.mod(X - cx * t, L), np.mod(Y - cy * t, L), L, **kwargs)
