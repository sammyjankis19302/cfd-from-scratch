# 2D Laplace equation solved with Jacobi iteration.
#
# My original hand-written version, kept unchanged as the starting point.
# The Jacobi update, the Dirichlet boundaries and the two-array pattern are
# all correct. Two things about it are documented in ../README.md:
#
#   1. `print("Ittt =", max_iterations + 1)` always prints 501 regardless of
#      what happened. It actually converged at iteration 152. Should be ti+1.
#
#   2. The stopping test uses max|s_new - s_old| -- the CHANGE per sweep --
#      rather than the residual, which is what "solved" actually means. For
#      Jacobi the two differ by exactly 4/dx^2, so this overestimates
#      convergence by 400x here and by 25,583x at N=81.
#
# Both are fixed in ../src/poisson.py.

import numpy as np
import matplotlib.pyplot as plt

L = 10
nx = 11
ny = 11
dx = L / (nx - 1)
dy = L / (ny - 1)
x = np.linspace(0, L, nx)
y = np.linspace(0, L, ny)
X, Y = np.meshgrid(x, y)
s = np.zeros((ny, nx))
print(s)
s[:, 0] = 0
s[:, -1] = 0
s[0, :] = 0
s[-1, :] = 1
tolerance = 1e-5
max_iterations = 500
for ti in range(max_iterations):
    s_old = s.copy()
    s_new = s.copy()
    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            s_new[j, i] = 0.25 * (
                s_old[j, i + 1] + s_old[j, i - 1] + s_old[j - 1, i] + s_old[j + 1, i]
            )
    error = np.max(np.abs(s_new - s_old))
    # Update solution
    s = s_new.copy()
    # Check convergence
    if error < tolerance:
        print("Converged  letsss partyyyyy!!!")
        print("Ittt =", max_iterations + 1)
        print("Final error =", error)
        break
else:
    print("Did not converge bruh.")
    print("Final error ===", error)
plt.contourf(X, Y, s)
plt.colorbar()
plt.show()
