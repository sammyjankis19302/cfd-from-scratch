# 1D linear advection: du/dt + c*du/dx = 0, first-order upwind (FTBS)
# Same code as my first script, just split into functions and with the bugs fixed.
# See BUGS.md for what was wrong and why.

import numpy as np
import matplotlib.pyplot as plt


def grid(L, nx):
    """nx points covering [0, L). Endpoint left out: on a ring x=L is x=0."""
    dx = L / nx
    x = np.zeros(nx)
    for i in range(nx):
        x[i] = i * dx
    return x, dx


def initial(x):
    """u(x,0) = sin(pi*x)"""
    return np.sin(np.pi * x)


def exact(x, c, t):
    """u(x,t) = sin(pi*(x - c*t)). Derived in theory/01-analytical-solution.md"""
    return np.sin(np.pi * (x - c * t))


def step(u_old, CFL):
    """One timestep. Every point is updated; point 0 takes its left
    neighbour from the end of the grid (periodic)."""
    nx = len(u_old)
    u_new = np.zeros(nx)
    for i in range(nx):
        if i == 0:
            left = nx - 1
        else:
            left = i - 1
        u_new[i] = u_old[i] - CFL * (u_old[i] - u_old[left])
    return u_new


def solve(u, CFL, nt):
    """March nt steps."""
    for n in range(nt):
        u = step(u, CFL)
    return u


def main():
    L, nx, nt, dt, c = 10.0, 100, 1000, 0.001, 1.0

    x, dx = grid(L, nx)
    CFL = c * dt / dx
    T = nt * dt
    print("c =", c, " dx =", dx, " dt =", dt, " CFL =", CFL, " T =", T)

    u = solve(initial(x), CFL, nt)
    s = exact(x, c, T)
    e = np.abs(s - u)          # absolute error, not signed

    print("numerical  max =", round(u.max(), 6), " min =", round(u.min(), 6))
    print("analytical max =", round(s.max(), 6), " min =", round(s.min(), 6))
    print("amplitude lost =", round(100 * (1 - u.max() / s.max()), 2), "%")
    print("max error      =", round(e.max(), 6))

    plt.plot(x, s, label="analytical", linewidth=8, color="lightblue")
    plt.plot(x, u, label="numerical", color="red")
    plt.title("1D advection   du/dt + c*du/dx = 0   (CFL = %g)" % CFL)
    plt.xlabel("x")
    plt.ylabel("u")
    plt.legend()
    plt.savefig("advection_1d.png", dpi=150)
    plt.show()


main()
