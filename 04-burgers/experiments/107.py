import numpy as np
import matplotlib.pyplot as plt

L = 10
nx = 101
nt = 1000

dt = 0.001
dx = L/(nx-1)

nu = .00001

x = np.linspace(0, L, nx)

u = np.sin(np.pi*x)

for n in range(nt):

    u_old = u.copy()

    for i in range(1, nx-1):

        # nonlinear advection
        if u_old[i] > 0:

            advection = (
                u_old[i] * dt/dx
                * (u_old[i] - u_old[i-1])
            )

        else:

            advection = (
                u_old[i] * dt/dx
                * (u_old[i+1] - u_old[i])
            )

        # diffusion
        diffusion = (
            nu * dt/dx**2
            * (
                u_old[i+1]
                - 2*u_old[i]
                + u_old[i-1]
            )
        )

        u[i] = u_old[i] - advection + diffusion

plt.plot(x, u)
plt.xlabel("x")
plt.ylabel("u")
plt.title("1D Burgers equation")
plt.show()
