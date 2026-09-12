import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

# ============================================================
# 1D LINEAR ADVECTION EQUATION
#
#               u_t + a u_x = 0
#
# Crank-Nicolson in time
# Central Difference in space
#
# Periodic boundary conditions
# ============================================================


# ------------------------------------------------------------
# PARAMETERS
# ------------------------------------------------------------

a = 1.0              # Advection speed
L = 1.0              # Length of domain

Nx = 400             # Number of grid points
dx = L / Nx

CFL = 0.8
dt = CFL * dx / abs(a)

T_final = 0.4

# Number of time steps
Nt = int(np.ceil(T_final / dt))

# Adjust dt so final time is exact
dt = T_final / Nt

# Actual CFL number
C = a * dt / dx

print("dx =", dx)
print("dt =", dt)
print("CFL =", C)
print("Nt =", Nt)


# ------------------------------------------------------------
# GRID
# ------------------------------------------------------------

x = np.linspace(0.0, L, Nx, endpoint=False)


# ------------------------------------------------------------
# INITIAL CONDITION
#
# Square pulse
# ------------------------------------------------------------

u = np.zeros(Nx)

u[(x >= 0.20) & (x < 0.40)] = 1.0

u_initial = u.copy()


# ============================================================
# CONSTRUCT MATRICES
#
# CN scheme:
#
# (u^(n+1) - u^n)/dt
#       + a/2 [D0 u^(n+1) + D0 u^n] = 0
#
# Therefore:
#
#       A u^(n+1) = B u^n
#
# where
#
# A = I + (C/2) * (central difference operator * dx)
# B = I - (C/2) * (central difference operator * dx)
#
# Explicitly:
#
# A[j,j+1] =  C/4
# A[j,j-1] = -C/4
#
# B[j,j+1] = -C/4
# B[j,j-1] =  C/4
# ============================================================


# Main diagonals
main = np.ones(Nx)

# Periodic upper/lower diagonals
upper_A = (C / 4.0) * np.ones(Nx - 1)
lower_A = (-C / 4.0) * np.ones(Nx - 1)

upper_B = (-C / 4.0) * np.ones(Nx - 1)
lower_B = (C / 4.0) * np.ones(Nx - 1)


# Construct sparse matrices
A = diags(
    diagonals=[main, upper_A, lower_A],
    offsets=[0, 1, -1],
    shape=(Nx, Nx),
    format='lil'
)

B = diags(
    diagonals=[main, upper_B, lower_B],
    offsets=[0, 1, -1],
    shape=(Nx, Nx),
    format='lil'
)


# ------------------------------------------------------------
# PERIODIC BOUNDARY TERMS
# ------------------------------------------------------------

# Matrix A
A[0, Nx - 1] = -C / 4.0
A[Nx - 1, 0] = C / 4.0

# Matrix B
B[0, Nx - 1] = C / 4.0
B[Nx - 1, 0] = -C / 4.0


# Convert to efficient sparse format
A = A.tocsc()
B = B.tocsc()


# ------------------------------------------------------------
# TIME INTEGRATION
# ------------------------------------------------------------

for n in range(Nt):

    # Right-hand side
    rhs = B @ u

    # Solve A u^(n+1) = rhs
    u = spsolve(A, rhs)


# ------------------------------------------------------------
# EXACT SOLUTION
#
# u(x,t) = u0(x - a*t)
#
# Periodic domain
# ------------------------------------------------------------

x_shift = (x - a * T_final) % L

u_exact = np.zeros(Nx)

u_exact[(x_shift >= 0.20) & (x_shift < 0.40)] = 1.0


# ------------------------------------------------------------
# PLOT
# ------------------------------------------------------------

plt.figure(figsize=(10, 5))

plt.plot(
    x,
    u_initial,
    '--',
    linewidth=2,
    label='Initial condition'
)

plt.plot(
    x,
    u_exact,
    linewidth=2,
    label='Exact solution'
)

plt.plot(
    x,
    u,
    linewidth=1.5,
    label='Crank-Nicolson + Central Difference'
)

plt.xlabel('x')
plt.ylabel('u')

plt.title(
    '1D Advection: Purely Dispersive Numerical Error\n'
    f'CN + Central Difference, CFL = {C:.2f}'
)

plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
