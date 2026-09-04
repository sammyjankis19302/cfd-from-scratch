"""Burgers equation, 1D and 2D, with conditional upwinding.

    1D:  u_t + u u_x = nu u_xx
    2D:  u_t + u u_x + v u_y = nu (u_xx + u_yy)
         v_t + u v_x + v v_y = nu (v_xx + v_yy)

WHAT IS NEW COMPARED WITH LINEAR ADVECTION

1. The equation is NONLINEAR. The advection speed is the solution itself, so
   fast parts of the wave catch up with slow parts. A perfectly smooth initial
   condition can develop a discontinuity in finite time -- a shock. Nothing in
   modules 01 or 03 could do that.

2. The advection speed CHANGES SIGN across the domain. In linear advection with
   c > 0 the upwind direction was fixed and known in advance. Here it must be
   decided per grid point, per timestep, from the local sign of the velocity:

       u > 0  ->  information arrives from the LEFT   ->  backward difference
       u < 0  ->  information arrives from the RIGHT  ->  forward difference

   This is CONDITIONAL UPWINDING. Getting it wrong -- using a fixed direction --
   gives a scheme that is upwind on half the domain and downwind on the other
   half, and the downwind half blows up.

3. The viscous term nu*u_xx imposes a SECOND stability constraint, unrelated to
   the CFL condition. See theory/01-burgers.md; both must be satisfied.

In 2D note carefully WHICH velocity selects the direction: the x-derivative of
BOTH u and v is upwinded on the sign of u, because u is the advecting velocity
in x. Likewise the y-derivative of both uses the sign of v. Using v to upwind a
v_x term is a natural-looking mistake that is wrong.
"""

import numpy as np


# ---------------------------------------------------------------------------
# 1D
# ---------------------------------------------------------------------------
def step_1d(u, dx, dt, nu):
    """One step of conditionally-upwinded Burgers in 1D. Periodic."""
    back = (u - np.roll(u, 1)) / dx          # u_x looking left
    fwd = (np.roll(u, -1) - u) / dx          # u_x looking right
    ux = np.where(u > 0, back, fwd)
    uxx = (np.roll(u, -1) - 2.0 * u + np.roll(u, 1)) / dx**2
    return u - dt * u * ux + nu * dt * uxx


def solve_1d(u0, dx, dt, nu, n_steps):
    u = u0.copy()
    for _ in range(n_steps):
        u = step_1d(u, dx, dt, nu)
    return u


# ---------------------------------------------------------------------------
# 2D
# ---------------------------------------------------------------------------
def step_2d(u, v, dx, dy, dt, nu, dirichlet=True):
    """One step of 2D Burgers. Vectorised form of experiments/first_attempt_2d.py.

    Arrays are indexed [j, i]: row j is y, column i is x.

    dirichlet=True zeroes all four edges, matching the loop version. With
    dirichlet=False the boundaries are periodic (np.roll handles them already).
    """
    ux = np.where(u > 0,
                  (u - np.roll(u, 1, axis=1)) / dx,
                  (np.roll(u, -1, axis=1) - u) / dx)
    uy = np.where(v > 0,
                  (u - np.roll(u, 1, axis=0)) / dy,
                  (np.roll(u, -1, axis=0) - u) / dy)
    vx = np.where(u > 0,
                  (v - np.roll(v, 1, axis=1)) / dx,
                  (np.roll(v, -1, axis=1) - v) / dx)
    vy = np.where(v > 0,
                  (v - np.roll(v, 1, axis=0)) / dy,
                  (np.roll(v, -1, axis=0) - v) / dy)

    uxx = (np.roll(u, -1, axis=1) - 2.0 * u + np.roll(u, 1, axis=1)) / dx**2
    uyy = (np.roll(u, -1, axis=0) - 2.0 * u + np.roll(u, 1, axis=0)) / dy**2
    vxx = (np.roll(v, -1, axis=1) - 2.0 * v + np.roll(v, 1, axis=1)) / dx**2
    vyy = (np.roll(v, -1, axis=0) - 2.0 * v + np.roll(v, 1, axis=0)) / dy**2

    un = u - dt * (u * ux + v * uy) + nu * dt * (uxx + uyy)
    vn = v - dt * (u * vx + v * vy) + nu * dt * (vxx + vyy)

    if dirichlet:
        for a in (un, vn):
            a[0, :] = a[-1, :] = a[:, 0] = a[:, -1] = 0.0

    return un, vn


def solve_2d(u0, v0, dx, dy, dt, nu, n_steps, dirichlet=True):
    u, v = u0.copy(), v0.copy()
    for _ in range(n_steps):
        u, v = step_2d(u, v, dx, dy, dt, nu, dirichlet)
    return u, v


# ---------------------------------------------------------------------------
# Stability limits -- see theory/01-burgers.md
# ---------------------------------------------------------------------------
def stability_limits(u_max, v_max, dx, dy, nu):
    """Return (advective dt limit, diffusive dt limit).

    Advective (CFL):   dt <= 1 / (|u|/dx + |v|/dy)
    Diffusive:         dt <= 0.5 / (nu (1/dx^2 + 1/dy^2))

    The two are independent and BOTH must hold. Which one binds depends on nu
    and dx: refining the grid tightens the diffusive limit as dx^2 but the
    advective limit only as dx, so on fine grids the viscous term is usually
    what forces a small timestep.
    """
    dt_adv = 1.0 / (abs(u_max) / dx + abs(v_max) / dy) if (u_max or v_max) else np.inf
    dt_dif = 0.5 / (nu * (1.0 / dx**2 + 1.0 / dy**2)) if nu > 0 else np.inf
    return dt_adv, dt_dif
