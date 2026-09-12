import numpy as np
import matplotlib.pyplot as plt
L=10
nt=1000
dt=0.001
nx=1001
dx=L/(nx-1)
x=np.linspace(0,L,nx)
# print(x,L)
u=np.sin(np.pi*x)
print(u)
CFL=np.max(np.abs(u))*(dt/dx)
print(CFL)
for i in range(nt):
    u_old=u.copy()
    for j in range(1,nx-1):
        if u_old[j]>0:
            u[j]=u_old[j]-(u_old[j]*(dt/dx)*(u_old[j]-u_old[j-1]))
        else:
            u[j]=u_old[j]-(u_old[j]*(dt/dx)*(u_old[j+1]-u_old[j]))
plt.plot(x,u)
plt.show()
