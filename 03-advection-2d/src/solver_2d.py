"""Time-marching driver for 2D linear advection.

Same separation as the 1D module: the scheme knows how to take one step, the
driver owns the grid, the initial condition, the timestep and the stopping
criterion.
"""

import numpy as np

from .initial_conditions_2d import make_grid_2d


def solve_2d(step, f, N, Cx, Cy, T, cx=1.0, cy=1.0, L=1.0, **f_kwargs):
    """March du/dt + cx du/dx + cy du/dy = 0 to time T.

    dt is derived from Cx (the grid is square, so dx = dy and the two Courant
    numbers are consistent by construction as long as cy/cx = Cy/Cx).

    Returns X, Y, u, t_final. As in 1D, an integer number of steps will not
    generally land exactly on T; rather than shortening the last step and
    invalidating the error analysis, the driver reports where it stopped.
    """
    X, Y, dx = make_grid_2d(N, L)
    dt = Cx * dx / abs(cx)
    n_steps = int(round(T / dt))
    t_final = n_steps * dt

    u = f(X, Y, L, **f_kwargs)
    for _ in range(n_steps):
        u = step(u, Cx, Cy)

    return X, Y, u, t_final


def l2_error_2d(u_num, u_exact, dx):
    """Grid-normalised L2 norm in 2D.

    The measure is dx*dy = dx^2 here, so the normalisation factor is dx, not
    sqrt(dx) as in 1D. Getting this wrong does not change the measured ORDER
    of convergence, but it does change the reported error magnitudes.
    """
    return dx * np.sqrt(np.sum((u_num - u_exact) ** 2))
