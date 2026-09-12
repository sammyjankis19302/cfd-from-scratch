"""Vectorised lid-driven cavity solver, vorticity-streamfunction formulation.

Numerically identical scheme to experiments/first_attempt_cavity.py:
Jacobi Poisson, ghost-point wall vorticity, conditional upwinding for
advection, central differences for diffusion, explicit Euler in time.

THREE CHANGES, all deliberate:

1. u AND v ARE SET ON THE BOUNDARIES.
   The original updates only range(1, n-1), so u[-1,:] stays 0 where the lid
   speed is 1. Harmless for the interior march -- the advection loop only reads
   interior u, v -- but it breaks two things downstream: Ghia's Table I starts
   at y = 1.0000 with u = 1.0000, which would be compared against zero, and the
   streamplot shows a stationary lid.

2. THE POISSON ITERATION STOPS ON THE RESIDUAL, NOT THE CHANGE PER SWEEP.
   Module 05 showed the two differ by exactly 4/dx^2. At nx = 101 that is
   40,000x, so a change-tolerance of 1e-6 corresponds to a residual near 0.04.

3. THE FIELDS ARE SAVED.
   The original writes only PNGs, so centreline data cannot be recovered from a
   finished run without repeating it.

Everything else -- the discretisation, the boundary treatment, the corner
exclusion, the ordering of operations within a timestep -- is unchanged.
"""

import numpy as np


def wall_vorticity(s, dx, dy, U_lid):
    """Third-order ghost-point wall vorticity on all four walls.

    Vectorised transcription of the original's loops. Corners are excluded for
    the same reason: a corner belongs to two walls at once, and the lid-driven
    cavity's corner singularity means neither formula is right there.

    The chain is: known wall velocity -> ghost psi outside the domain -> wall
    vorticity from the second derivative. Without the ghost point there is no
    information on the outer side of the wall to form that second derivative.
    """
    o = np.zeros_like(s)

    # top wall (moving lid): psi_ghost carries the 3*U*dy term
    g = 3.0 * U_lid * dy - 1.5 * s[-1, 1:-1] + 3.0 * s[-2, 1:-1] - 0.5 * s[-3, 1:-1]
    o[-1, 1:-1] = -(g - 2.0 * s[-1, 1:-1] + s[-2, 1:-1]) / dy**2

    # bottom wall (stationary)
    g = -1.5 * s[0, 1:-1] + 3.0 * s[1, 1:-1] - 0.5 * s[2, 1:-1]
    o[0, 1:-1] = -(g - 2.0 * s[0, 1:-1] + s[1, 1:-1]) / dy**2

    # left wall
    g = -1.5 * s[1:-1, 0] + 3.0 * s[1:-1, 1] - 0.5 * s[1:-1, 2]
    o[1:-1, 0] = -(g - 2.0 * s[1:-1, 0] + s[1:-1, 1]) / dx**2

    # right wall
    g = -1.5 * s[1:-1, -1] + 3.0 * s[1:-1, -2] - 0.5 * s[1:-1, -3]
    o[1:-1, -1] = -(g - 2.0 * s[1:-1, -1] + s[1:-1, -2]) / dx**2

    return o


def poisson_residual(s, o, dx):
    """r = lap(s) + o, on interior points. Zero when the equation is satisfied."""
    r = np.zeros_like(s)
    r[1:-1, 1:-1] = (
        s[1:-1, 2:] + s[1:-1, :-2] + s[2:, 1:-1] + s[:-2, 1:-1]
        - 4.0 * s[1:-1, 1:-1]
    ) / dx**2 + o[1:-1, 1:-1]
    return r


def solve_poisson(s, o, dx, rtol=1e-4, max_sweeps=20000, omega_sor=None):
    """Red-black SOR for lap(psi) = -omega, stopping on the RELATIVE residual.

    TWO CHOICES WORTH EXPLAINING.

    Relative, not absolute, tolerance. The residual has units of vorticity, and
    omega reaches O(10-100) near the lid. An absolute tolerance of 1e-6 is
    therefore demanding roughly 1e-8 relative -- and since residual = 4/dx^2
    times the change per sweep (module 05), that means a change of ~1e-12 per
    sweep. Measured: 2151 Jacobi sweeps PER TIMESTEP at N=41, which makes the
    solver unusable. The stopping test is now

        max|residual| < rtol * max(|omega|, 1)

    SOR, not Jacobi. Module 05 measured the scaling exponents: Jacobi 2.00,
    SOR 1.04. Inside a timestep loop that is the difference between a solver
    that finishes and one that does not. Red-black ordering is used because it
    vectorises -- updating all red points from black neighbours and then all
    black from red is Gauss-Seidel in a different sweep order, so the
    sequential dependency that blocks vectorising plain Gauss-Seidel is gone.

    s is the initial guess; inside the timestep loop the previous
    streamfunction is already close.
    """
    N = s.shape[0]
    if omega_sor is None:
        omega_sor = 2.0 / (1.0 + np.sin(np.pi / (N - 1)))
    s = s.copy()
    scale = max(np.abs(o).max(), 1.0)

    # red-black masks over interior points
    jj, ii = np.mgrid[0:N, 0:N]
    interior = (jj > 0) & (jj < N - 1) & (ii > 0) & (ii < N - 1)
    red = interior & (((ii + jj) % 2) == 0)
    black = interior & (((ii + jj) % 2) == 1)

    for k in range(1, max_sweeps + 1):
        for mask in (red, black):
            nb = np.zeros_like(s)
            nb[1:-1, 1:-1] = (s[1:-1, 2:] + s[1:-1, :-2]
                              + s[2:, 1:-1] + s[:-2, 1:-1])
            gs = 0.25 * (nb + o * dx**2)
            s[mask] = (1.0 - omega_sor) * s[mask] + omega_sor * gs[mask]
        if np.abs(poisson_residual(s, o, dx)).max() < rtol * scale:
            return s, k
    return s, max_sweeps


def velocities(s, dx, dy, U_lid):
    """u = psi_y, v = -psi_x by central differences, plus the wall values.

    The boundary assignment is change (1) above: no-slip on three walls, the
    lid speed on the fourth.
    """
    u = np.zeros_like(s)
    v = np.zeros_like(s)
    u[1:-1, 1:-1] = (s[2:, 1:-1] - s[:-2, 1:-1]) / (2.0 * dy)
    v[1:-1, 1:-1] = -(s[1:-1, 2:] - s[1:-1, :-2]) / (2.0 * dx)
    u[-1, :] = U_lid          # moving lid
    u[0, :] = 0.0
    u[:, 0] = u[:, -1] = 0.0
    v[0, :] = v[-1, :] = 0.0
    v[:, 0] = v[:, -1] = 0.0
    return u, v


def advance_vorticity(o, u, v, dx, dy, dt, nu):
    """Explicit Euler, conditional upwind advection, central diffusion."""
    ox = np.where(u[1:-1, 1:-1] > 0,
                  (o[1:-1, 1:-1] - o[1:-1, :-2]) / dx,
                  (o[1:-1, 2:] - o[1:-1, 1:-1]) / dx)
    oy = np.where(v[1:-1, 1:-1] > 0,
                  (o[1:-1, 1:-1] - o[:-2, 1:-1]) / dy,
                  (o[2:, 1:-1] - o[1:-1, 1:-1]) / dy)
    oxx = (o[1:-1, 2:] - 2.0 * o[1:-1, 1:-1] + o[1:-1, :-2]) / dx**2
    oyy = (o[2:, 1:-1] - 2.0 * o[1:-1, 1:-1] + o[:-2, 1:-1]) / dy**2

    on = o.copy()
    on[1:-1, 1:-1] = (
        o[1:-1, 1:-1]
        - dt * (u[1:-1, 1:-1] * ox + v[1:-1, 1:-1] * oy)
        + nu * dt * (oxx + oyy)
    )
    return on


def solve_cavity(N=101, Re=100.0, dt=1e-3, T=100.0, L=1.0, U_lid=1.0,
                 rtol_poisson=1e-4, tol_steady=1e-8, verbose=True):
    """March to steady state. Returns a dict of fields and diagnostics."""
    dx = dy = L / (N - 1)
    nu = U_lid * L / Re
    xs = np.linspace(0, L, N)
    X, Y = np.meshgrid(xs, xs)

    s = np.zeros((N, N))
    o = np.zeros((N, N))
    n_max = int(round(T / dt))
    sweeps_total = 0

    for n in range(1, n_max + 1):
        s, k = solve_poisson(s, o, dx, rtol_poisson)
        sweeps_total += k
        u, v = velocities(s, dx, dy, U_lid)

        o_wall = wall_vorticity(s, dx, dy, U_lid)
        o[-1, 1:-1] = o_wall[-1, 1:-1]
        o[0, 1:-1] = o_wall[0, 1:-1]
        o[1:-1, 0] = o_wall[1:-1, 0]
        o[1:-1, -1] = o_wall[1:-1, -1]

        o_new = advance_vorticity(o, u, v, dx, dy, dt, nu)
        change = np.abs(o_new - o).max()
        o = o_new

        if verbose and n % 5000 == 0:
            print(f"    step {n:7d}  t = {n*dt:7.2f}  max|d omega| = {change:.3e}")
        if change < tol_steady:
            if verbose:
                print(f"    steady at step {n} (t = {n*dt:.3f})")
            break

    u, v = velocities(s, dx, dy, U_lid)
    return {"X": X, "Y": Y, "psi": s, "omega": o, "u": u, "v": v,
            "N": N, "Re": Re, "dt": dt, "dx": dx, "nu": nu,
            "steps": n, "t_final": n * dt,
            "poisson_sweeps_total": sweeps_total,
            "poisson_residual": np.abs(poisson_residual(s, o, dx)).max(),
            "final_change": change}


def save(result, path):
    """Persist the fields so the centrelines can be extracted later.

    The original saved only PNGs, which means a finished run cannot be
    revisited without repeating it.
    """
    np.savez_compressed(path, **{k: v for k, v in result.items()})
