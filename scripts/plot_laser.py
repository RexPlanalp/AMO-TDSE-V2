import numpy as np
import matplotlib.pyplot as plt

laser_data = np.loadtxt("misc/laser.txt")
t = laser_data[:,0]
Ax = laser_data[:,1]
Ay = laser_data[:,2]
Az = laser_data[:,3]

plt.plot(t,Ax,color = "k",label = "Ax")
plt.plot(t,Ay,color = "brown",label = "Ay")
plt.plot(t,Az,color = "blue",label = "Az")
plt.legend()
plt.savefig("images/laser.png")