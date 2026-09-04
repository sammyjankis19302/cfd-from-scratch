"""Time-marching driver for the 1D linear advection equation.

Separating this from schemes.py is deliberate. The scheme knows how to take one
step; the driver owns the grid, the initial condition, the timestep and the
stopping criterion. Neither needs to know anything about the other beyond the
step(u, C) contract.
"""

import numpy as np

from .initial_conditions import make_grid


def solve(step, f, N, C, T, c=1.0, L=1.0, **f_kwargs):
    """March u_t + c u_x = 0 from t = 0 to t = T.

    Parameters
    ----------
    step : callable
        One of the schemes from schemes.py, with signature step(u, C).
    f : callable
        Initial condition, called as f(x, L, **f_kwargs).
    N : int
        Number of grid cells.
    C : float
        Courant number. dt is derived from it, not the other way round --
        C is the parameter the stability analysis is written in terms of.
    T : float
        Final time.

    Returns
    -------
    x      : ndarray, the grid
    u      : ndarray, the solution at t = t_final
    t_final: float, the time actually reached

    Note on t_final: an integer number of steps of size dt will not in general
    land exactly on T. Rather than shortening the last step -- which would
    change C for that step and quietly invalidate the error analysis -- the
    driver stops at the last whole step and reports where it stopped. The
    caller compares against the exact solution at t_final, not at T.
    """
    x, dx = make_grid(N, L)
    dt = C * dx / abs(c)
    n_steps = int(round(T / dt))
    t_final = n_steps * dt

    u = f(x, L, **f_kwargs)
    for _ in range(n_steps):
        u = step(u, C)

    return x, u, t_final


def l2_error(u_num, u_exact, dx):
    """Grid-normalised L2 norm of the error.

    The sqrt(dx) factor is what makes this a norm of a *function* rather than
    of a vector of numbers. Without it the value would drift as N changes for
    reasons that have nothing to do with accuracy, and the measured order of
    convergence would be wrong.
    """
    return np.sqrt(dx * np.sum((u_num - u_exact) ** 2))
