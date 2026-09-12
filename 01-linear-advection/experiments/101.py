import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# 1D LINEAR ADVECTION EQUATION
#
#       u_t + a u_x = 0
#
# BTCS:
#
# (u_j^(n+1) - u_j^n)/dt
#     + a (u_(j+1)^(n+1) - u_(j-1)^(n+1))/(2 dx) = 0
#
# ============================================================


# ------------------------------------------------------------
# PARAMETERS
# ------------------------------------------------------------

a = 1.0                  # Advection velocity
L = 1.0                  # Domain length

Nx = 200                 # Number of grid points
dx = L / Nx

CFL = 2.0                # Can be > 1: BTCS is unconditionally stable
dt = CFL * dx / abs(a)

T_final = 0.5            # Final simulation time
Nt = int(np.ceil(T_final / dt))

# Adjust dt so that we finish exactly at T_final
dt = T_final / Nt

# Actual CFL number
C = a * dt / dx

print("dx =", dx)
print("dt =", dt)
print("Number of time steps =", Nt)
print("CFL =", C)


# ------------------------------------------------------------
# GRID
# ------------------------------------------------------------

x = np.linspace(0.0, L, Nx, endpoint=False)


# ------------------------------------------------------------
# INITIAL CONDITION: JUMP / SQUARE PULSE
#
# u = 1 for 0.25 <= x < 0.50
# u = 0 elsewhere
#
# ------------------------------------------------------------

u = np.zeros(Nx)

u[(x >= 0.25) & (x < 0.50)] = 1.0

u_initial = u.copy()


# ------------------------------------------------------------
# BUILD BTCS MATRIX
#
# From:
#
# u_j^(n+1)
# + (C/2) u_(j+1)^(n+1)
# - (C/2) u_(j-1)^(n+1)
# = u_j^n
#
# So:
#
#              A u^(n+1) = u^n
#
# ------------------------------------------------------------

A = np.zeros((Nx, Nx))

for j in range(Nx):

    # Main diagonal
    A[j, j] = 1.0

    # Right neighbor: j+1
    A[j, (j + 1) % Nx] = C / 2.0

    # Left neighbor: j-1
    A[j, (j - 1) % Nx] = -C / 2.0


# ------------------------------------------------------------
# TIME INTEGRATION
# ------------------------------------------------------------

for n in range(Nt):

    # Solve:
    #
    # A u^(n+1) = u^n
    #
    u = np.linalg.solve(A, u)


# ------------------------------------------------------------
# EXACT SOLUTION
#
# For periodic advection:
#
# u(x,t) = u0(x - a*t)
#
# ------------------------------------------------------------

x_shifted = (x - a * T_final) % L

u_exact = np.zeros(Nx)

u_exact[(x_shifted >= 0.25) & (x_shifted < 0.50)] = 1.0


# ------------------------------------------------------------
# PLOT
# ------------------------------------------------------------

plt.figure(figsize=(9, 5))

plt.plot(
    x,
    u_initial,
    "--",
    linewidth=2,
    label="Initial condition"
)

plt.plot(
    x,
    u_exact,
    linewidth=2,
    label="Exact solution"
)

plt.plot(
    x,
    u,
    linewidth=2,
    label="BTCS + Central Difference"
)

plt.xlabel("x")
plt.ylabel("u")
plt.title(
    f"1D Advection: BTCS + Central Difference\n"
    f"Nx = {Nx}, CFL = {C:.2f}, T = {T_final}"
)

plt.legend()
plt.grid(True)
plt.tight_layout()

plt.show()
