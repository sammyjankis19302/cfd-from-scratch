## 1d 1st order Advection eqaution du/dt+cdu/dx=0 (Hopefully it works)
import numpy as np
import matplotlib.pyplot as plt
#def f(x):
 #   return sin(x)
L=10
nx=1001
nt=500             #deltaT=nt*dt
dt=0.001            #mooving =c*deltaT
dx=(L/(nx-1))
c=10.05
CFL=c*dt/dx
print('CFL=',CFL,"c=",c,)
x=np.linspace(0,L,nx)
u=np.sin(np.pi*x)
r=np.ones(L)
o=u.copy()
print("Maximum Amplitide of initial solution (Numeraical)",max(o),"Minimum is ",min(o))
# plt.plot(x,o)
# plt.show()
print(u)
for i in range(nt):
    u_new=u.copy()
    # u[0] = np.sin(-np.pi*c*nt*dt)
    # u[0]=np.sin(-np.pi*(i+1)*c*dt)
    for j in range(1,nx):
        u[j]=u_new[j]-CFL*(u_new[j]-u_new[j-1])
#print(u)

a=max(u)
print("Max Numeraical Amplitude of final solution is ",a,"Minimum is",min(u))
# plt.plot(x,u)
# plt.show()
s=np.sin(np.pi*(x-c*nt*dt)) 
#print(s)
print("Maximum analytical Amplitide is",max(s),"Minimum is",min(s))
# plt.plot(x,s)
# plt.show()
#print(len(s),len(o),len(u))
e=s-u
print(e)
print("Error max is ",max(e),"Min is ",min(e))
plt.plot(x,s,label="analytical solution",linewidth="9",color="Lightblue")  
# plt.plot(x,o,label="Initial solution",color="Green")
# plt.plot(x,e,label="Error",color="Orange")
plt.plot(x,u,label="Numeraical solution",color="Red")
plt.title("1D 1st Order advection equation ----du/dt+c*du/dx=0")
plt.legend()
plt.show()
