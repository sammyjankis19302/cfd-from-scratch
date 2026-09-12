import numpy as np
import matplotlib.pyplot as plt
import time
start_time = time.time()
L=1
nx=201
ny=201
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
nt=100000
tol=1e-8
for nn in range(nt):
    s_new=s.copy()
    for itt in range(10000):
        for j in range(1,ny-1):
            for i in range(1,nx-1):
                s_new[j,i]=0.25*(s[j+1,i]+s[j-1,i]+s[j,i+1]+s[j,i-1]+o[j,i]*dx**2)
        errrr=  np.max(np.abs(s_new-s))
        s=s_new.copy()
        if errrr<1e-6:
            break
    for j in range(1,ny-1):
        for i in range(1,nx-1):
            u[j,i]=(s[j+1,i]-s[j-1,i])/(2*dy)
            v[j,i]=-(s[j,i+1]-s[j,i-1])/(2*dx)                  #The minus sign broooooooooo extremely important dont even try to do some sht
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

        o[j,i]=-(s_ghost-2*s_i+s_i1)/(dx**2)                    # Wall vorticity from the moving-lid boundary condition V=0
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
    errrr=np.max(np.abs(o_new-o_old))
    if nn % 100==0:
        print(nn,errrr)
    o=o_new.copy()
    elapsed_time = time.time() - start_time
    if errrr < tol:
        print("Converged!")
        print("Iteration:", nn)
        print("Error:", errrr)
        print("Time taken:", elapsed_time, "seconds")
        print("Time taken:", elapsed_time / 60, "minutes")
        print("Time taken:", elapsed_time / 3600, "hours")
        break
plt.figure()
plt.contourf(X, Y, s,)
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Stream__function")
plt.savefig("pictures/Re_100_t100_201x/stream_function")
plt.show()

plt.figure()
plt.contourf(X, Y, o, )
plt.colorbar(label="Vorticity")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Vorticity")
plt.axis("equal")
plt.savefig("pictures/Re_100_t100_201x/vorticity_t100.png")
plt.show()

plt.figure()
plt.streamplot(X, Y, u, v)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Lid-Driven Cavity Flow")
plt.savefig("pictures/Re_100_t100_201x/v_adn_u_function_t100.png")
plt.show()






