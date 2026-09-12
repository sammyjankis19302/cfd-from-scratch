import numpy as np
import matplotlib.pyplot as plt

# Parameters
a = 1.0
L = 1.0

Nx = 200
dx = L / Nx

CFL = 0.2
dt = CFL * dx / a

T_final = 0.3
Nt = int(T_final / dt)

x = np.linspace(0, L, Nx, endpoint=False)

# Square pulse
u = np.zeros(Nx)
u[(x >= 0.2) & (x < 0.4)] = 1.0

u0 = u.copy()

# FTCS
for n in range(Nt):

    u = (
        u
        - 0.5 * CFL
        * (np.roll(u, -1) - np.roll(u, 1))
    )

# Plot
plt.figure(figsize=(10, 5))

plt.plot(x, u0, "--", label="Initial")
plt.plot(x, u, label="FTCS")

plt.xlabel("x")
plt.ylabel("u")
plt.title("FTCS for 1D Advection — Unstable")
plt.legend()
plt.grid()
plt.show()
