import matplotlib.pyplot as plt
import numpy as np
K1=np.linspace(0.000000000000001,np.pi,100)
# print(K)
exact=np.ones(len(K1))
print(len(K1)==len(exact))
C2D2_real=(np.sin(K1/2)/(K1/2))**2
C2D4_real=(30-32*np.cos(K1)+2*np.cos(2*K1))/(12*K1**2)
F2D_real = (np.sin(K1/2)/(K1/2))**2 * np.cos(K1)
B2D_real = (np.sin(K1/2)/(K1/2))**2 * np.cos(K1)
F2D_img = (np.sin(K1/2)/(K1/2))**2 * np.sin(K1)
B2D_img = -(np.sin(K1/2)/(K1/2))**2 * np.sin(K1)

plt.plot(K1, np.zeros_like(K1), label="Exact_img,CD-2_img,CD-2_img")
plt.plot(K1,exact,label="Exact")
plt.plot(K1,F2D_img,label="FD_img")
plt.plot(K1,B2D_img,label="BD_img")
plt.plot(K1,C2D2_real,label="CD2_real",color="cyan",linewidth=8)
plt.plot(K1,C2D4_real,label="CD4_real",color="green",linewidth=8)
plt.plot(K1,B2D_real,label="BD_real",color="red",linewidth=5)
plt.plot(K1,F2D_real,label="FD_real",color="blue",linewidth=2)
plt.legend()
plt.title("2nd Order")
plt.show()

