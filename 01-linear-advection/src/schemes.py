"""Finite-difference schemes for u_t + c u_x = 0, c > 0.

Every scheme has the same signature:

    step(u, C) -> u_new

where `u` is the solution at time level n and `C = c dt / dx` is the Courant
number. One step advances the solution by one dt. Nothing here knows about the
grid, the initial condition or the final time -- that is the driver's job.

The uniform signature is the load-bearing design decision in this module.
Because every scheme looks identical from outside, the convergence study, the
spectral analysis and every experiment can loop over a dictionary of schemes
without a single special case.

Periodic boundaries are handled by np.roll, which wraps the array end-to-end:

    np.roll(u,  1)[i] == u[i-1]     (shift right)
    np.roll(u, -1)[i] == u[i+1]     (shift left)

so the boundary condition costs nothing and is applied identically everywhere.
Nothing below writes into `u`; each scheme builds a new array. That is the
two-array discipline -- updating in place would silently feed already-updated
values back into the stencil.
"""

import numpy as np


def ftbs(u, C):
    """First-order upwind (forward time, backward space).

    Stable for 0 <= C <= 1. Modified equation carries a +u_xx term with
    nu_num = c dx (1 - C) / 2, so the error is dissipative: peaks flatten,
    corners round off. Exact at C = 1, where nu_num vanishes.
    """
    return u - C * (u - np.roll(u, 1))


def ftfs(u, C):
    """Downwind (forward time, forward space). Unconditionally unstable for c > 0.

    Kept deliberately. It is the cleanest demonstration that the direction of
    the stencil is a physical statement about where information comes from,
    not a matter of taste.
    """
    return u - C * (np.roll(u, -1) - u)


def ftcs(u, C):
    """Forward time, centred space. Unconditionally unstable.

    Second-order accurate in space and useless on its own: the Von Neumann
    amplification factor has modulus sqrt(1 + C^2 sin^2 theta) > 1 for every
    non-zero wavenumber.
    """
    return u - 0.5 * C * (np.roll(u, -1) - np.roll(u, 1))


def lax_friedrichs(u, C):
    """Lax-Friedrichs: FTCS with u[i] replaced by the neighbour average.

    That replacement is an added dissipation of exactly the size needed to
    stabilise FTCS -- which is also why it is far more diffusive than upwind
    at the same C. Stable for |C| <= 1.
    """
    return 0.5 * (np.roll(u, -1) + np.roll(u, 1)) - 0.5 * C * (
        np.roll(u, -1) - np.roll(u, 1)
    )


def lax_wendroff(u, C):
    """Lax-Wendroff. Second-order in space and time, stable for |C| <= 1.

    Leading error term is dispersive (u_xxx) rather than dissipative, with a
    coefficient proportional to (1 - C^2). So it keeps a square wave's corners
    sharp but grows trailing oscillations -- and both the oscillations and the
    error vanish at C = 1.
    """
    return (
        u
        - 0.5 * C * (np.roll(u, -1) - np.roll(u, 1))
        + 0.5 * C**2 * (np.roll(u, -1) - 2.0 * u + np.roll(u, 1))
    )


def bd2(u, C):
    """Second-order upwind (Beam-Warming). Stable for 0 <= C <= 2.

    Uses u[i], u[i-1], u[i-2] -- still entirely upwind, but second order.

    The second term is not optional decoration. A second-order upwind stencil
    marched with plain forward Euler is *unconditionally unstable*, exactly like
    FTCS: the spatial order buys nothing if the time discretisation is left at
    first order. The u^2 term is the Lax-Wendroff-style time correction that
    makes the scheme second order in time as well, and stable.

    Formally more accurate than FTBS, yet by global spectral analysis it needs
    *more* points per wavelength than FTBS for the same phase accuracy. Order
    of accuracy and resolving power are not the same property.
    """
    return (
        u
        - 0.5 * C * (3.0 * u - 4.0 * np.roll(u, 1) + np.roll(u, 2))
        + 0.5 * C**2 * (u - 2.0 * np.roll(u, 1) + np.roll(u, 2))
    )


#: Registry used by the driver, the tests and every experiment.
SCHEMES = {
    "FTBS": ftbs,
    "FTFS": ftfs,
    "FTCS": ftcs,
    "Lax-Friedrichs": lax_friedrichs,
    "Lax-Wendroff": lax_wendroff,
    "BD-2": bd2,
}

#: Formal order of accuracy in space, for the convergence study to check against.
#: None means the scheme does not converge at all.
FORMAL_ORDER = {
    "FTBS": 1,
    "FTFS": None,
    "FTCS": None,
    "Lax-Friedrichs": 1,
    "Lax-Wendroff": 2,
    "BD-2": 2,
}
