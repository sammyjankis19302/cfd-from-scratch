# NOTE: itt=5 below. Five Jacobi sweeps cannot converge a 101x101 Poisson
# problem -- it needs a few thousand (see ../README.md: Jacobi is O(N^2)).
# So this always prints 'diverged heheheh'. Looks like a debug leftover;
# 5000 was probably intended. Preserved exactly as written.
#
# Also note: this is the first script here with a NON-ZERO source term,
# o[20:31,20:31]=1 -- the step from Laplace to Poisson, and what the
# cavity's grad^2(psi) = -omega actually needs.
import numpy as np
import matplotlib.pyplot as plt
L=10
nx=101
ny=101
dx=L/(nx-1)
dy=L/(ny-1)
x=np.linspace(0,L,nx)
y=np.linspace(0,L,ny)
X,Y=np.meshgrid(x,y)
s=np.zeros((nx,ny))
s[0,:]=0
s[:,0]=0
s[-1,:]=0
s[:,-1]=0
o=np.zeros((nx,ny))
o[20:31,20:31]=1
tol=1e-5
itt=5
for tt in range(itt):
    s_old=s.copy()
    s_new=s.copy()
    for j in range(1,ny-1):
        for i in range(1,nx-1):
            s_new[j,i]=0.25*(s_old[j-1,i]+s_old[j+1,i]+s_old[j,i+1]+s_old[j,i-1]+o[j,i]*dx**2)
    err=np.max(np.abs(s_new-s_old))
    s=s_new.copy()
    if err<tol:
        print("converged","\nerror=",err,"\nMax itterations=",itt+1)
        break
else:
    print("diverged heheheh","\nerror=",err)
plt.contourf(X, Y, s)
plt.colorbar(label="psi value")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Streamfunction from Poisson Equation")
plt.show()



