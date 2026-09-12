import numpy as np
import matplotlib.pyplot as plt
L=10
nx=101
nt=1000
dt=0.001
dx=L/(nx-1)
D=0.1      #diffusion constant
x=np.linspace(0,L,nx)
u=np.sin(np.pi*x)
for i in range(nt):
    u_old=u.copy()
    for j in range(1,nx-1):
        u[j]=u_old[j]-(u_old[j]*(dt/dx)*(u_old[j]-u_old[j-1]))+D*(dt/dx**2)*(u_old[j+1]-2*u_old[j]+u_old[j-1])
plt.plot(x,u)
plt.show()
Re=(L*u)/D #Showing us the difference between advection and diffusion term ------ Reynold's number gives us this (advection/diffusion),u * u_x ~ U**2 / L,nu * u_xx ~ nu * U / L**2,Re = U * L / nu
print(Re)      #(U**2 / L) * u_t + (U**2 / L) * u * u_x = (nu * U / L**2) * u_xx    ,u_t + u * u_x = (1 / Re) * u_xx-----Non dinmentionized burgers eqaution 

