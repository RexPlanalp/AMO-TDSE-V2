import numpy as np
import matplotlib.pyplot as plt

laserData = np.loadtxt("misc/laser.txt")

t = laserData[:,0]
Ax = laserData[:,1]
Ay = laserData[:,2]
Az = laserData[:,3]

fig,ax = plt.subplots()


ax.plot(t,Ax,color = "k",label = "Ax",linewidth = 0.8)
ax.plot(t,Ay,color = "brown",label = "Ay",linewidth = 0.8)
ax.plot(t,Az,color = "blue",label = "Az",linewidth = 0.8)
ax.set_ylabel("Vector Potential")
ax.set_xlabel("Time (au)")
ax.legend()
fig.savefig("images/laser.png",dpi=200)