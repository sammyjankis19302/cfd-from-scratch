# NOTE: line 1 below reads 'mport numpy as np' -- a missing 'i'. This file
# does not parse. Preserved exactly as written; the fix is one character.
mport numpy as np
import matplotlib.pyplot as plt
L=1
nx=101
ny=101
dx=L/(nx-1)
dy=L/(ny-1)
x=np.linspace(0,L,nx)
y=np.linspace(0,L,ny)
X,Y=np.meshgrid(x,y)
s=np.zeros((ny,nx))
o=np.zeros((ny,nx))
u=np.zeros((ny,nx))
v=np.zeros((ny,nx))
s[0,:]=0
s[:,0]=0
s[-1,:]=0
s[:,-1]=0
U=1
j=ny-1
for i in range(1,nx-1):
    s_j=s[j,i]
    s_j1=s[j-1,i]
    s_j2=s[j-2,i]

    s_ghost=(3*U*dy-1.5*s_j+3*s_j1-0.5*s_j2)                #above the top wall using Ghia paper

    o[j,i]=-(s_ghost-2*s_j+s_j1)/(dy**2)                    # Wall vorticity from the moving-lid boundary condition U = 1
U=0
j=0
for i in range(1,nx-1):
    s_j=s[j,i]
    s_j1=s[j+1,i]
    s_j2=s[j+2,i]

    s_ghost=(-1.5*s_j+3*s_j1-0.5*s_j2)                      #below the bottom wall using Ghia paper

    o[j,i]=-(s_ghost-2*s_j+s_j1)/(dy**2)                    # Wall vorticity from the moving-lid boundary condition U = 0
V=0
i=0
for j in range(1,ny-1):
    s_i=s[j,i]
    s_i1=s[j,i+1]
    s_i2=s[j,i+2]

    s_ghost=(-1.5*s_i+3*s_i1-0.5*s_i2)                      #Left side of the left wall using Ghia paper

    o[j,i]=-(s_ghost-2*s_i+s_i1)/(dx**2)                    # Wall vorticity from the moving-lid boundary condition V = 0
V=0
i=nx-1
for j in range(1,ny-1):
    s_i=s[j,i]
    s_i1=s[j,i-1]
    s_i2=s[j,i-2]   

    s_ghost=(-1.5*s_i+3*s_i1-0.5*s_i2)                      #Right side of the right wall using Ghia paper

    o[j,i]=-(s_ghost-2*s_i+s_i1)/(dx**2)                    # Wall vorticity from the moving-lid boundary condition V=  0
for j in range(1,ny-1):
    for i in range(1,nx-1):
        u[j,i]=(s[j+1,i]-s[j-1,i])/(2*dy)
        v[j,i]=-(s[j,i+1]-s[j,i-1])/(2*dx)                  #The minus sign broooooooooo extremely important dont even try to do some sht
o_new=o.copy()
o_old=o.copy()
dt = 0.001
nu = 0.01
for j in range(1,ny-1):
    for i in range(1,nx-1):
        if u[j,i]>0:
            o_x=(o_old[j,i]-o_old[j,i-1])/dx
        else:
            o_x=(o_old[j,i+1]-o_old[j,i])/dx
        if v[j,i]>0:
            o_y=(o_old[j,i]-o_old[j-1,i])/dy
        else:
            o_y=(o_old[j+1,i]-o_old[j,i])/dy
        o_xx=(o_old[j,i+1]-2*o_old[j,i]+o_old[j,i-1])/(dx**2)
        o_yy=(o_old[j+1,i]-2*o_old[j,i]+o_old[j-1,i])/(dy**2)
        o_new[j,i]=o_old[j,i]-dt*u[j,i]*o_x-dt*v[j,i]*o_y+nu*dt*(o_xx+o_yy)
o=o_new.copy()
s_new=s.copy()
for itt in range(1000):
    for j in range(1,ny-1):
        for i in range(1,nx-1):
            s_new[j,i]=0.25*(s[j+1,i]+s[j-1,i]+s[j,i+1]+s[j,i-1]+o[j,i]*dx**2)
    s=s_new.copy()
    






        






