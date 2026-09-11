"""How much viscosity does first-order upwind ADD to the cavity solver?

This is the module's central question, and it decides how the Ghia comparison
should be read.

THE ARGUMENT

Module 01 derived the modified equation for FTBS: the scheme does not solve

    u_t + c u_x = 0

it solves

    u_t + c u_x = nu_num u_xx,    nu_num = |c| dx (1 - C) / 2,    C = |c| dt/dx

That nu_num is not a metaphor. It is a diffusion coefficient with the same
units and the same effect as physical viscosity, and in the vorticity transport
equation it sits directly alongside the physical nu.

So the solver is not running at

    Re = U L / nu

it is running at something closer to

    Re_eff = U L / (nu + nu_num)

CONSEQUENCE

nu_num depends on dx, not on nu. So as the target Re rises (nu falls) at fixed
grid, nu_num stays put and eventually dominates. At Re = 1000 on a 101x101 grid
with |u| ~ 1 near the lid, nu_num is about 4.5x the physical viscosity.

This is why agreement with Ghia should be good at Re = 100 and degrade sharply
with Re -- and it is predictable BEFORE running anything, which turns a
disappointing result into a quantitative one.

The cell Reynolds number Re_cell = |u| dx / nu is the standard way to express
the same thing: once Re_cell exceeds ~2, first-order upwind's added diffusion is
comparable to the physical diffusion it is supposed to be resolving.
"""

import numpy as np


def numerical_viscosity(u_mag, dx, dt):
    """nu_num for first-order upwind, from the module 01 modified equation.

    Parameters
    ----------
    u_mag : local speed |u| (the advecting velocity, not the advected quantity)
    dx, dt: grid spacing and timestep

    The (1 - C) factor means nu_num is LARGEST as dt -> 0. Shrinking the
    timestep to be safe makes the numerical diffusion worse, not better.
    Refining dx is what reduces it.
    """
    C = u_mag * dt / dx
    return u_mag * dx * (1.0 - C) / 2.0


def cell_reynolds(u_mag, dx, nu):
    """Re_cell = |u| dx / nu -- the grid Peclet number for this problem.

    Rule of thumb: below ~2 the physical diffusion is resolved on the grid;
    above it, the upwind scheme's own diffusion takes over.
    """
    return u_mag * dx / nu


def effective_reynolds(U, L, nu, u_mag, dx, dt):
    """Re based on total (physical + numerical) diffusion."""
    return U * L / (nu + numerical_viscosity(u_mag, dx, dt))


def report(cases, U=1.0, L=1.0, dt=1e-3, u_mag=None):
    """Print the table for a list of (label, nu, N) cases.

    u_mag defaults to 1.0 (lid speed -- the worst case, near the moving wall).
    Pass u_mag=0.2 for a representative interior speed.
    """
    um = 1.0 if u_mag is None else u_mag
    print(f"{'case':10s}{'grid':>8s}{'dx':>9s}{'nu':>10s}"
          f"{'nu_num':>10s}{'ratio':>8s}{'Re_eff':>9s}{'Re_cell':>9s}")
    print("-" * 73)
    rows = []
    for label, nu, N in cases:
        dx = L / (N - 1)
        nn = numerical_viscosity(um, dx, dt)
        rows.append((label, nu, N, dx, nn, nn / nu,
                     effective_reynolds(U, L, nu, um, dx, dt),
                     cell_reynolds(um, dx, nu)))
        print(f"{label:10s}{N:6d}^2{dx:9.4f}{nu:10.4f}"
              f"{nn:10.5f}{nn/nu:8.2f}{rows[-1][6]:9.0f}{rows[-1][7]:9.2f}")
    return rows
