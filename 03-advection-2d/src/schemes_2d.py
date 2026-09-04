"""Finite-difference schemes for 2D linear advection, cx, cy > 0.

    du/dt + cx du/dx + cy du/dy = 0

Uniform signature, same idea as the 1D module:

    step(u, Cx, Cy) -> u_new

where Cx = cx dt / dx and Cy = cy dt / dy.

Indexing is u[j, i]: j indexes y (rows), i indexes x (columns). So

    np.roll(u, 1, axis=1)   is u[j, i-1]   -- neighbour in -x
    np.roll(u, 1, axis=0)   is u[j-1, i]   -- neighbour in -y

THE RESULT WORTH KNOWING. In 1D, upwind is stable for C <= 1. The natural guess
in 2D is that Cx <= 1 and Cy <= 1 independently. That guess is WRONG. The true
condition is

    Cx + Cy <= 1

so at Cx = Cy the limit per direction is 0.5, not 1. Von Neumann analysis in
theory/01-2d-advection.md shows why: the amplification factor accumulates
dissipation from both directions at once, and it is the SUM that must stay
bounded. tests/test_convergence_2d.py demonstrates this numerically -- the
scheme is stable at Cx = Cy = 0.45 and diverges at Cx = Cy = 0.55, either side
of the predicted 0.5.
"""

import numpy as np


def upwind_2d(u, Cx, Cy):
    """First-order upwind in both directions. Stable for Cx + Cy <= 1.

    Dimension-by-dimension: each direction contributes its own backward
    difference, and they simply add. Leading error is dissipative in both x
    and y.
    """
    return (
        u
        - Cx * (u - np.roll(u, 1, axis=1))
        - Cy * (u - np.roll(u, 1, axis=0))
    )


def lax_wendroff_2d(u, Cx, Cy):
    """Two-dimensional Lax-Wendroff, including the cross-derivative term.

    Obtained from the 2D Taylor expansion in time: u_tt now generates
    u_xx, u_yy AND a mixed u_xy term. Dropping the cross term is a common
    mistake -- it leaves the scheme first order along diagonals even though it
    looks second order along each axis. Since the exact solution here travels
    diagonally, this module would expose that immediately.

    Stable roughly for Cx^2 + Cy^2 <= 1.
    """
    ux_p, ux_m = np.roll(u, -1, axis=1), np.roll(u, 1, axis=1)
    uy_p, uy_m = np.roll(u, -1, axis=0), np.roll(u, 1, axis=0)

    # mixed second difference u_xy, four-corner stencil
    u_pp = np.roll(ux_p, -1, axis=0)
    u_pm = np.roll(ux_p, 1, axis=0)
    u_mp = np.roll(ux_m, -1, axis=0)
    u_mm = np.roll(ux_m, 1, axis=0)

    return (
        u
        - 0.5 * Cx * (ux_p - ux_m)
        - 0.5 * Cy * (uy_p - uy_m)
        + 0.5 * Cx**2 * (ux_p - 2.0 * u + ux_m)
        + 0.5 * Cy**2 * (uy_p - 2.0 * u + uy_m)
        + 0.25 * Cx * Cy * (u_pp - u_pm - u_mp + u_mm)
    )


def lax_wendroff_2d_no_cross(u, Cx, Cy):
    """Lax-Wendroff with the cross term deliberately omitted.

    Kept as a negative control. It looks correct and converges at second order
    when the flow is aligned with an axis, but drops toward first order for
    diagonal advection. Included so the cross term's necessity is demonstrated
    rather than asserted.
    """
    ux_p, ux_m = np.roll(u, -1, axis=1), np.roll(u, 1, axis=1)
    uy_p, uy_m = np.roll(u, -1, axis=0), np.roll(u, 1, axis=0)
    return (
        u
        - 0.5 * Cx * (ux_p - ux_m)
        - 0.5 * Cy * (uy_p - uy_m)
        + 0.5 * Cx**2 * (ux_p - 2.0 * u + ux_m)
        + 0.5 * Cy**2 * (uy_p - 2.0 * u + uy_m)
    )


SCHEMES_2D = {
    "Upwind-2D": upwind_2d,
    "Lax-Wendroff-2D": lax_wendroff_2d,
    "LW-2D (no cross)": lax_wendroff_2d_no_cross,
}

FORMAL_ORDER_2D = {
    "Upwind-2D": 1,
    "Lax-Wendroff-2D": 2,
    "LW-2D (no cross)": None,   # second order on axes, degrades on diagonals
}
