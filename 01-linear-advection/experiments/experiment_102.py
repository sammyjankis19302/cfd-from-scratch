## 1d 1st order Advection eqaution du/dt+cdu/dx=0 (Hopefully it works)
import numpy as np
import matplotlib.pyplot as plt
#def f(x):
 #   return sin(x)
L=10
nx=10001
nx1=1000001
nt=100             #deltaT=nt*dt
dt=0.0001            #mooving =c*deltaT
dx=(L/(nx-1))
dx1=(L/(nx1-1))
c=-1
c1=-.1
CFL=c*dt/dx

CFL1=c1*dt/dx1
print('CFL=',CFL,"c=",c,)

print('CFL2=',CFL1,"c=",c1,)
x=np.linspace(0,L,nx)
x1 = np.linspace(0,L,nx1)
print(x)
u=np.sin(50*np.pi*x)*np.ex;(-5*(x-0.5*L)**2)
u1=np.sin(50*np.pi*x1)*np.exp(-5*(x1-0.5*L)**2)
r=np.ones(L)
o=u.copy()
print("Maximum Amplitide of initial solution (Numeraical)",max(o),"Minimum is ",min(o))
# plt.plot(x,o)
# plt.show()
#print(u)

for i in range(nt):
    u_new=u.copy()
    # u[nx-1] = np.sin(-np.pi*c*nt*dt)
    u[nx-1]=np.sin(150*np.pi*(L-((i+1)*c*dt)))*np.exp(-5*(L-0.5*L)**2)
    for j in range(0,nx-1):
        u[j]=u_new[j]-c*((dt/dx)*(u_new[j+1]-u_new[j]))
#print(u)

for i in range(nt):
    u_new1=u1.copy()
    # u[nx-1] = np.sin(-np.pi*c*nt*dt)
    u1[nx1-1]=np.sin(150*np.pi*(L-((i+1)*c1*dt)))*np.exp(-5*(L-0.5*L)**2)
    for j in range(0,nx1-1):
        u1[j]=u_new1[j]-c1*((dt/dx1)*(u_new1[j+1]-u_new1[j]))
#print(u)
print(u[L])
a=max(u)
print("Max Numeraical Amplitude of final solution is ",a,"Minimum is",min(u))
# plt.plot(x,u)
# plt.show()
s=np.sin(150*np.pi*(x-c*nt*dt))*np.exp(-5*(x-c*nt*dt-0.5*L)**2) 
#print(s)
print("Maximum analytical Amplitide is",max(s),"Minimum is",min(s))
# plt.plot(x,s)
# plt.show()
#print(len(s),len(o),len(u))
e=s-u
print(e)
print("Error max is ",max(e),"Min is ",min(e))
plt.plot(x,s,label="analytical solution",linewidth="1",color="Blue")  

# plt.plot(x,o,label="Initial solution",color="Green")
# plt.plot(x,e,label="Error",color="Orange")
plt.plot(x,u,label="Numeraical solution",color="Red")
plt.plot(x1,u1,label="Sir's Numeraical solution",color="Green")
plt.title("1D 1st Order advection equation ----du/dt+c*du/dx=0")
plt.legend()
plt.show()


