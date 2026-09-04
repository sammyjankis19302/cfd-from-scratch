"""Initial conditions for the 1D linear advection equation.

Each function returns u(x, 0) on the grid `x`. Every one of them is periodic
on [0, L), so it can be advected round the domain indefinitely without the
boundary condition introducing an error of its own.

The choice of initial condition is a diagnostic decision, not a physical one:
the equation transports anything, so f is picked to expose a specific defect
of the scheme.

    gaussian     smooth  -> Taylor-based error analysis applies; use for
                            measuring order of accuracy and for the
                            width-growth test of numerical viscosity
    square       C^0     -> Taylor analysis breaks down; exposes dissipative
                            corner-rounding and dispersive oscillations
    sine         single k-> one wavenumber only; use for spectral analysis
    wave_packet  narrow  -> envelope and carrier travel together in the exact
                            solution, so any splitting is group-velocity error
"""

import numpy as np


def gaussian(x, L, x0=0.5, sigma=0.05):
    """Gaussian bump of width sigma centred at x0 * L."""
    return np.exp(-((x - x0 * L) ** 2) / (2.0 * sigma**2))


def square(x, L, x0=0.3, x1=0.5):
    """Top-hat equal to 1 on [x0*L, x1*L) and 0 elsewhere."""
    return np.where((x >= x0 * L) & (x < x1 * L), 1.0, 0.0)


def sine(x, L, k_modes=1):
    """Single Fourier mode with `k_modes` full wavelengths in the domain."""
    return np.sin(2.0 * np.pi * k_modes * x / L)


def wave_packet(x, L, x0=0.5, sigma=0.08, k_modes=12):
    """Sine carrier under a Gaussian envelope."""
    envelope = np.exp(-((x - x0 * L) ** 2) / (2.0 * sigma**2))
    carrier = np.sin(2.0 * np.pi * k_modes * x / L)
    return envelope * carrier


def make_grid(N, L=1.0):
    """Uniform periodic grid: N cells on [0, L).

    The final point L is deliberately excluded. On a periodic domain x = L is
    the same point as x = 0, so including both would duplicate a degree of
    freedom and quietly break conservation.

    Returns
    -------
    x  : ndarray, shape (N,)
    dx : float
    """
    dx = L / N
    x = np.linspace(0.0, L - dx, N)
    return x, dx


def exact_solution(f, x, L, c, t, **kwargs):
    """Exact solution u(x, t) = f(x - c t) with periodic wrap-around.

    This exists because the whole verification strategy rests on it: for the
    linear advection equation the true answer is known for *any* f, so every
    deviation the solver produces is pure numerical error.
    """
    return f(np.mod(x - c * t, L), L, **kwargs)
