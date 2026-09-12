#2d non-linear advection diffusion with upwinding on botht the velocities (2D Burger's eqaution )(This is a bit tricky boyyy)(Hopefully it works).
import numpy as np
import matplotlib.pyplot as plt
L=10
nx=101
ny=101
nt=100
dt=0.001
dx=L/(nx-1)
dy=L/(ny-1)
nu=0.001
x=np.linspace(0,L,nx)
y=np.linspace(0,L,ny)
X,Y=np.meshgrid(x,y)
u=np.sin(np.pi*X)*np.sin(np.pi*Y)
v=np.sin(np.pi*X)*np.sin(np.pi*Y)
for n in range(nt) :
    u_old=u.copy()
    v_old=v.copy()
    for j in range(1,ny-1):
        for i in range(1,nx-1):
            if u_old[j,i]>0:
                ux=(u_old[j,i]-u_old[j,i-1])/dx  
            else:
                ux=(u_old[j,i+1]-u_old[j,i])/dx
            if v_old[j,i]>0:
                uy=(u_old[j,i]-u_old[j-1,i])/dy 
            else :
                uy=(u_old[j+1,i]-u_old[j,i])/dy 
            uxx=(u_old[j,i+1]-2*u_old[j,i]+u_old[j,i-1])/dx**2
            uyy=(u_old[j+1,i]-2*u_old[j,i]+u_old[j-1,i])/dy**2
            u[j,i]=u_old[j,i]-dt*(ux*u_old[j,i]+uy*v_old[j,i])+nu*(uxx+uyy)*dt
            if u_old[j,i]>0:
                vx=(v_old[j,i]-v_old[j,i-1])/dx  
            else:
                vx=(v_old[j,i+1]-v_old[j,i])/dx
            if v_old[j,i]>0:
                vy=(v_old[j,i]-v_old[j-1,i])/dy 
            else :
                vy=(v_old[j+1,i]-v_old[j,i])/dy 
            vxx=(v_old[j,i+1]-2*v_old[j,i]+v_old[j,i-1])/dx**2
            vyy=(v_old[j+1,i]-2*v_old[j,i]+v_old[j-1,i])/dy**2
            v[j,i]=v_old[j,i]-dt*(vx*u_old[j,i]+vy*v_old[j,i])+nu*(vxx+vyy)*dt
plt.figure()
plt.contourf(X,Y,u)
plt.colorbar()
plt.show()
plt.figure()
plt.contourf(X,Y,v)
plt.colorbar()
plt.show()

                




