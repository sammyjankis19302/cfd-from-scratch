#Vorticity is 0 which means that del**2 phi=0, so 
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
tolerance = 1e-10
max_iterations = 500
for ti in range(max_iterations):
    s_old = s.copy()
    s_new = s.copy()
    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            s_new[j, i] = 0.25 * (
                s_old[j, i + 1] + s_old[j, i - 1] + s_old[j - 1, i] + s_old[j + 1, i]
            )
    err = np.max(np.abs(s_new - s_old))
    s = s_new.copy()
    if err < tolerance:
        print("Convergeddddddd  ")
        print("Ittt =", max_iterations + 1)
        print("Final errorrrr =", err)
        break
else:
    print("Did not convergeeee")
    print("Final errorrr ===", err)
plt.contourf(X, Y, s)
plt.colorbar()
plt.show()
