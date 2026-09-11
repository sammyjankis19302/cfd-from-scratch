import matplotlib.pyplot as plt
import numpy as np
K=np.linspace(0.000000000000001,np.pi,100)
# print(K)
exact=np.ones(len(K))
print(len(K)==len(exact))
CD2_real=np.sin(K)/K
FD_real=np.sin(K)/K
BD_real=np.sin(K)/K
CD4_real=(8*np.sin(K)-np.sin(2*K)) / (6*K)
FD_img=(1-np.cos(K))/K
BD_img=-(1-np.cos(K))/K
plt.plot(K, np.zeros_like(K), label="Exact_img,CD-2_img,CD-4_img")
plt.plot(K,exact,label="Exact")
plt.plot(K,FD_img,label="FD_img")
plt.plot(K,BD_img,label="BD_img")
plt.plot(K,CD2_real,label="CD2_real",color="cyan",linewidth=8)
plt.plot(K,CD4_real,label="CD4_real",color="green",linewidth=8)
plt.plot(K,BD_real,label="BD_real",color="red",linewidth=5)
plt.plot(K,FD_real,label="FD_real",color="blue",linewidth=2)
plt.legend()
plt.title("1st Order")
plt.show()






