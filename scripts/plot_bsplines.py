import numpy as np
import matplotlib.pyplot as plt
import json


with open("input.json") as f:
    data = json.load(f)

n_bspline = data["BSpline"]["nbasis"]


bspline_data = np.loadtxt("misc/bsplines.txt")
dbspline_data = np.loadtxt("misc/dbsplines.txt")
meta_data = np.loadtxt("misc/bspline_metadata.txt")

Nr,dr = int(meta_data[0]),meta_data[1]

real = bspline_data[:,0]
imag = bspline_data[:,1]

dreal = dbspline_data[:,0]
dimag = dbspline_data[:,1]

r = np.array([i * dr for i in range(Nr)])



fig,ax = plt.subplots()
for i in range(n_bspline):
    ax.plot(r, real[i*Nr:(i+1)*Nr], color='black',linewidth = 0.8)  
    ax.plot(r, imag[i*Nr:(i+1)*Nr], color='brown',linewidth = 0.8)  
ax.set_xlabel('r (au)')
ax.set_ylabel("BSpline")
ax.tight_layout()
fig.savefig("images/bsplines.png",dpi = 200)


fig,ax = plt.subplots()
for i in range(n_bspline):
    ax.plot(r, dreal[i*Nr:(i+1)*Nr], color='black',linewidth = 0.8)  
    ax.plot(r, dimag[i*Nr:(i+1)*Nr], color='brown',linewidth = 0.8)  
ax.set_xlabel('r (au)')
ax.set_ylabel("d-BSpline")
ax.tight_layout()
ax.savefig("images/dbsplines.png",dpi = 200)
   



