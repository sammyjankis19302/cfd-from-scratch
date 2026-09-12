import numpy as np
import matplotlib.pyplot as plt

# --------------------------------
# Parameters
# --------------------------------
L = 10.0
nx = 1001

x = np.linspace(0, L, nx)
dx = x[1] - x[0]

dt = 0.0002
nt = 2500

# --------------------------------
# Initial condition
# --------------------------------
u = np.sin(np.pi * x)

# Times at which we want to save the solution
save_times = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]

solutions = {}
solutions[0.0] = u.copy()

# --------------------------------
# Time integration
# --------------------------------
for n in range(1, nt + 1):

    u_old = u.copy()

    # Lax-Friedrichs scheme
    for j in range(1, nx - 1):

        flux_right = 0.5 * u_old[j + 1]**2
        flux_left  = 0.5 * u_old[j - 1]**2

        u[j] = (
            0.5 * (u_old[j + 1] + u_old[j - 1])
            - (dt / (2 * dx)) * (flux_right - flux_left)
        )

    # Boundary conditions
    u[0] = 0.0
    u[-1] = 0.0

    # Current time
    t = n * dt

    # Save requested times
    for save_time in save_times:
        if save_time not in solutions:
            if abs(t - save_time) < dt / 2:
                solutions[save_time] = u.copy()

# --------------------------------
# Plot
# --------------------------------
plt.figure(figsize=(10, 6))

for t in save_times:
    plt.plot(x, solutions[t], label=f"t = {t}")

plt.xlabel("x")
plt.ylabel("u")
plt.title("Inviscid Burgers' Equation")
plt.legend()
plt.grid()
plt.show()
