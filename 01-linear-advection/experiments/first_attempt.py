# 1D_advection linear equation with sin wave as the initial condition
import numpy as np
import matplotlib.pyplot as plt

L = 10
nx = 101
nt = 1000
dt = 0.001
dx = L / (nx - 1)
x = np.linspace(0, L, nx)
u = np.sin(
    np.pi * x
)  # Initial solution or my sin wave , I can keep it anything I need to advect
c = 1
print("Speed=", c)
CFL = c * dt / dx
print("CFl=", CFL)
for j in range(nt):
    u_old = u.copy()
    for i in range(1, nx - 1):
        u[i] = u_old[i] - (CFL * (u_old[i] - u_old[i - 1]))
# T = nt * dt
# u_analytical = np.sin(np.pi * (x - c * T))
# print("max of numerical is ", max(u))
# print("max of analytical is ", max(u_analytical))
# plt.plot(x, u, label="numerical")
# plt.plot(x, u_analytical, label="analytical")
# plt.savefig("1dadvection")
a = max(u)
print("Max Numeraical Amplitude of final solution is ", a, "Minimum is", min(u))
# plt.plot(x,u)
# plt.show()
s = np.sin(np.pi * (x - c * nt * dt))
# print(s)
print("Maximum analytical Amplitide is", max(s), "Minimum is", min(s))
# plt.plot(x,s)
# plt.show()
# print(len(s),len(o),len(u))
e = s - u
print(e)
print("Error max is ", max(e), "Min is ", min(e))
plt.plot(x, s, label="analytical solution", linewidth="9", color="Lightblue")
# plt.plot(x,o,label="Initial solution",color="Green")
# plt.plot(x,e,label="Error",color="Orange")
plt.plot(x, u, label="Numeraical solution", color="Red")
plt.title("1D 1st Order advection equation ----du/dt+c*du/dx=0")
plt.legend()
plt.savefig("hehe")
