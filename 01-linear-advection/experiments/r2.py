import numpy as np
import matplotlib.pyplot as plt
L = 10
nx = 1001

dx = L/(nx-1)

c = -1
dt = 0.001

T = 1
nt = int(T/dt)

x = np.linspace(0,L,nx)

frequencies = [2, 10, 50]

plt.figure(figsize=(10,6))

for f in frequencies:

    # Initial condition
    u = np.sin(f*np.pi*x)

    for i in range(nt):

        u_new = u.copy()

        for j in range(nx-1):
            u[j] = u_new[j] - c*(dt/dx)*(u_new[j+1]-u_new[j])

        # Boundary condition
        u[-1] = np.sin(f*np.pi*(L-c*(i+1)*dt))

    # Analytical solution
    s = np.sin(f*np.pi*(x-c*T))

    plt.plot(x,s,'--',label=f'Analytical f={f}')
    plt.plot(x,u,label=f'Numerical f={f}')

plt.xlabel("x")
plt.ylabel("u")
plt.legend()
plt.grid()
plt.show()
