import numpy as np
import matplotlib.pyplot as plt

laser_data = np.loadtxt("misc/laser.txt")
t = laser_data[:,0]
Ax = laser_data[:,1]
Ay = laser_data[:,2]
Az = laser_data[:,3]



plt.plot(t,Ax,color = "k",label = "Ax",linewidth = 0.8)
plt.plot(t,Ay,color = "brown",label = "Ay",linewidth = 0.8)
plt.plot(t,Az,color = "blue",label = "Az",linewidth = 0.8)
plt.legend()
plt.savefig("images/laser.png",dpi=200)