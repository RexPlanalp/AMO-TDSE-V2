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



plt.figure(figsize=(8,8))
for i in range(n_bspline):
    plt.plot(r, real[i*Nr:(i+1)*Nr], color='black',linewidth = 0.8)  # real part
    plt.plot(r, imag[i*Nr:(i+1)*Nr], color='brown',linewidth = 0.8)  # imag part
plt.xlabel('Point index')
plt.ylabel('Value')
plt.title('B-spline basis functions (real:black, imag:brown)')
plt.tight_layout()
plt.show()
plt.savefig("images/bsplines.png",dpi = 200)


plt.figure(figsize=(8,8))
for i in range(n_bspline):
    plt.plot(r, dreal[i*Nr:(i+1)*Nr], color='black',linewidth = 0.8)  # real part
    plt.plot(r, dimag[i*Nr:(i+1)*Nr], color='brown',linewidth = 0.8)  # imag part
plt.xlabel('Point index')
plt.ylabel('Value')
plt.title('B-spline basis functions (real:black, imag:brown)')
plt.tight_layout()
plt.show()
plt.savefig("images/dbsplines.png",dpi = 200)
   



