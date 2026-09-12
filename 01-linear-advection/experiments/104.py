import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# 1D ADVECTION: LEAPFROG + CENTRAL DIFFERENCE
#
# u_t + a u_x = 0
#
# u_j^(n+1) = u_j^(n-1)
#             - C [u_(j+1)^n - u_(j-1)^n]
#
# Periodic boundary conditions
# ============================================================

# Parameters
a = 1.0
L = 1.0

Nx = 400
dx = L / Nx

CFL = 0.8
dt = CFL * dx / abs(a)

T_final = 0.4
Nt = int(np.ceil(T_final / dt))

dt = T_final / Nt
C = a * dt / dx

# Grid
x = np.linspace(0.0, L, Nx, endpoint=False)

# Initial square pulse
u0 = np.zeros(Nx)
u0[(x >= 0.2) & (x < 0.4)] = 1.0

# u at time n = 0
u_old = u0.copy()

# ============================================================
# STARTUP STEP
# Lax-Friedrichs:
#
# u_j^1 = 0.5(u_(j+1)^0 + u_(j-1)^0)
#         - C/2 (u_(j+1)^0 - u_(j-1)^0)
# ============================================================

u = (
    0.5 * (np.roll(u_old, -1) + np.roll(u_old, 1))
    - 0.5 * C * (np.roll(u_old, -1) - np.roll(u_old, 1))
)

# Leapfrog time integration
for n in range(1, Nt):

    u_new = (
        u_old
        - C * (np.roll(u, -1) - np.roll(u, 1))
    )

    u_old = u
    u = u_new


# Exact solution
x_shift = (x - a * T_final) % L

u_exact = np.zeros(Nx)
u_exact[(x_shift >= 0.2) & (x_shift < 0.4)] = 1.0


# Plot
plt.figure(figsize=(10, 5))

plt.plot(x, u0, "--", linewidth=2, label="Initial")
plt.plot(x, u_exact, linewidth=2, label="Exact")
plt.plot(x, u, linewidth=1.5,
         label="Leapfrog + Central Difference")

plt.xlabel("x")
plt.ylabel("u")
plt.title(
    f"1D Advection: Leapfrog Scheme\n"
    f"CFL = {C:.2f}"
)

plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
