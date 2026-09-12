import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

# ============================================================
# 1D ADVECTION EQUATION
#
#       u_t + a u_x = 0
#
# BTCS + Central Difference:
#
# (u_j^(n+1) - u_j^n)/dt
#   + a [u_(j+1)^(n+1) - u_(j-1)^(n+1)]/(2dx) = 0
#
# Periodic boundary conditions
# ============================================================


# -------------------------
# PARAMETERS
# -------------------------
a = 1.0
L = 1.0

Nx = 400
dx = L / Nx

CFL = 2.0
dt = CFL * dx / abs(a)

T_final = 0.4
Nt = int(np.ceil(T_final / dt))

# Adjust dt to reach exactly T_final
dt = T_final / Nt
C = a * dt / dx

print(f"dx  = {dx}")
print(f"dt  = {dt}")
print(f"CFL = {C}")
print(f"Nt  = {Nt}")


# -------------------------
# SPATIAL GRID
# -------------------------
x = np.linspace(0, L, Nx, endpoint=False)


# -------------------------
# INITIAL CONDITION
# Square pulse / jump
# -------------------------
u = np.zeros(Nx)

u[(x >= 0.2) & (x <= 0.4)] = 1.0

u_initial = u.copy()


# -------------------------
# CONSTRUCT MATRIX A
#
# A u^(n+1) = u^n
#
# Main diagonal:       1
# Upper diagonal:     +C/2
# Lower diagonal:     -C/2
# -------------------------
A = lil_matrix((Nx, Nx))

for j in range(Nx):

    # u_j^(n+1)
    A[j, j] = 1.0

    # +(C/2) u_(j+1)^(n+1)
    A[j, (j + 1) % Nx] = C / 2.0

    # -(C/2) u_(j-1)^(n+1)
    A[j, (j - 1) % Nx] = -C / 2.0


# Convert to CSC format for efficient solving
A = A.tocsc()


# -------------------------
# TIME MARCHING
# -------------------------
for n in range(Nt):

    # Solve A u^(n+1) = u^n
    u = spsolve(A, u)


# -------------------------
# EXACT SOLUTION
#
# u(x,t) = u0(x - a*t)
# with periodic wrapping
# -------------------------
x0 = (x - a * T_final) % L

u_exact = np.zeros(Nx)

u_exact[(x0 >= 0.2) & (x0 <= 0.4)] = 1.0


# -------------------------
# PLOT
# -------------------------
plt.figure(figsize=(10, 5))

plt.plot(
    x,
    u_initial,
    '--',
    linewidth=2,
    label='Initial'
)

plt.plot(
    x,
    u_exact,
    linewidth=2,
    label='Exact'
)

plt.plot(
    x,
    u,
    linewidth=2,
    label='BTCS + CD'
)

plt.xlabel('x')
plt.ylabel('u')

plt.title(
    f'BTCS + Central Difference for 1D Advection\n'
    f'N = {Nx}, CFL = {C:.2f}'
)

plt.legend()
plt.grid()
plt.show()
