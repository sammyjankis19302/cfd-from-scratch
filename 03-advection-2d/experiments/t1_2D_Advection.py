# NOTE: this file does not parse. Line 55 reads
#     for j in range(1, ny):i
# The trailing 'i' makes the loop a one-liner, so the nested 'for k' block
# below it is an IndentationError. Preserved exactly as written.
# 2D_advection_linear.py in this folder is the working version.
import numpy as np
import matplotlib.pyplot as plt

# 2D 1st order Advection equation
# du/dt + cx*du/dx + cy*du/dy = 0

L = 10

nx = 101
ny = 101

nt = 500

dt = 0.001

dx = L/(nx-1)
dy = L/(ny-1)

cx = 10.0
cy = 10.0

CFLx = cx*dt/dx
CFLy = cy*dt/dy

print("CFLx =", CFLx)
print("CFLy =", CFLy)

# x and y coordinates
x = np.linspace(0, L, nx)
y = np.linspace(0, L, ny)

# Create 2D grid
X, Y = np.meshgrid(x, y)

# Initial condition
u = np.sin(np.pi*X) * np.sin(np.pi*Y)

# Save initial solution
o = u.copy()

print("Maximum Amplitude of initial solution =", np.max(o))
print("Minimum Amplitude of initial solution =", np.min(o))


# ------------------------------------------------
# Numerical solution
# ------------------------------------------------

for i in range(nt):

    u_new = u.copy()

    for j in range(1, ny):i

        for k in range(1, nx):

            u[j, k] = u_new[j, k] \
                     - CFLx*(u_new[j, k] - u_new[j, k-1]) \
                     - CFLy*(u_new[j, k] - u_new[j-1, k])


# ------------------------------------------------
# Final numerical solution
# ------------------------------------------------

print("Maximum Numerical Amplitude =", np.max(u))
print("Minimum Numerical Amplitude =", np.min(u))


# ------------------------------------------------
# Analytical solution
# ------------------------------------------------

t = nt*dt

s = np.sin(np.pi*(X-cx*t)) * np.sin(np.pi*(Y-cy*t))

print("Maximum Analytical Amplitude =", np.max(s))
print("Minimum Analytical Amplitude =", np.min(s))


# ------------------------------------------------
# Error
# ------------------------------------------------

e = s-u

print("Error max is =", np.max(e))
print("Error min is =", np.min(e))


# ------------------------------------------------
# Plot Numerical Solution
# ------------------------------------------------

plt.figure()
plt.contourf(X, Y, u, 50)
plt.colorbar()
plt.title("2D Numerical Solution")
plt.xlabel("x")
plt.ylabel("y")
plt.show()


# ------------------------------------------------
# Plot Analytical Solution
# ------------------------------------------------

plt.figure()
plt.contourf(X, Y, s, 50)
plt.colorbar()
plt.title("2D Analytical Solution")
plt.xlabel("x")
plt.ylabel("y")
plt.show()


# ------------------------------------------------
# Plot Error
# ------------------------------------------------

plt.figure()
plt.contourf(X, Y, e, 50)
plt.colorbar()
plt.title("Error")
plt.xlabel("x")
plt.ylabel("y")
plt.show()
