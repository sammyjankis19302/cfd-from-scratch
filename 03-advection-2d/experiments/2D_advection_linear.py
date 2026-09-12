# 2-D advection eq (Hopefully this also workssssssss.)     :)))))
import numpy as np
import matplotlib.pyplot as plt
L=10
nx=101
ny=101
dx=L/(nx-1)
dy=L/(ny-1)
cx=10
cy=10
dt=0.0004
nt=500
CFLx=cx*dt/dx
CFLy=cy*dt/dy
print("Total CFL =", CFLx + CFLy)
print("CFLx=",CFLx,"CFLy=",CFLy)
x=np.linspace(0,L,nx)
y=np.linspace(0,L,ny)
X,Y=np.meshgrid(x,y)
# print(X,Y)
u=np.sin(np.pi*X)*np.sin(np.pi*Y)
print(u)
for f in range(nt):
    u_old=u.copy()
    for j in range(1,ny):
        for i in range(1,nx):
            u[j,i]=u_old[j,i]-CFLx*(u_old[j,i]-u_old[j,i-1])-CFLy*(u_old[j,i]-u_old[j-1,i])
T=dt*nt
u_a=np.sin(np.pi*(X-cx*T))*np.sin(np.pi*(Y-cy*T))
print("u shape =", u.shape)
print("u_a shape =", u_a.shape)
plt.figure(figsize=(8,6))
plt.contourf(X,Y,u,levels=50)
plt.title("numerical")
plt.show()

plt.figure(figsize=(8,6))
plt.contourf(X,Y,u_a,levels=50)
plt.title("analytical")
plt.show()


