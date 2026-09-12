# lets do the 2d nonliner advection eq (omggggg)
import numpy as np
import matplotlib.pyplot as plt
L=10
nt=2000
dt=0.0001
nx=101
ny=101
dx=L/(nx-1)
dy=L/(ny-1)
x=np.linspace(0,L,nx)
y=np.linspace(0,L,ny)
X,Y = np.meshgrid(x,y)
u=np.sin(np.pi*X)*np.sin(np.pi*Y)
v=np.sin(np.pi*X)*np.sin(np.pi*Y)
print(u)
for b in range(1,nt):
    u_old=u.copy()
    v_old=v.copy()
    for j in range(1,ny):
        for i in range(1,nx):
            u[j,i]=u_old[j,i]-u_old[j,i]*(dt/dx)*(u_old[j,i]-u_old[j,i-1])-v_old[j,i]*(dt/dy)*(u_old[j,i]-u_old[j-1,i])
            v[j,i]=v_old[j,i]-u_old[j,i]*(dt/dx)*(v_old[j,i]-v_old[j,i-1])-v_old[j,i]*(dt/dy)*(v_old[j,i]-v_old[j-1,i])
print("u shape =", u.shape)
plt.figure(figsize=(8,6))
plt.contourf(X,Y,u,levels=50)
plt.title("numerical_nonlinear")
plt.show()


