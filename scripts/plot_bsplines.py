import re
import numpy as np
import matplotlib.pyplot as plt

bdata = []
with open('misc/bsplines.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        toks = re.findall(r'\(([^)]+)\)', line)
        row = [list(map(float, t.split(','))) for t in toks]
        bdata.append(row)

bdata = np.array(bdata) 

dbdata = []
with open('misc/dbsplines.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        toks = re.findall(r'\(([^)]+)\)', line)
        row = [list(map(float, t.split(','))) for t in toks]
        dbdata.append(row)

dbdata = np.array(dbdata) 


n_basis, n_pts, _ = bdata.shape
x = np.arange(n_pts)          


plt.figure(figsize=(8,8))
for i in range(n_basis):
    plt.plot(x, bdata[i,:,0], color='black',linewidth = 0.8)  # real part
    plt.plot(x, bdata[i,:,1], color='brown',linewidth = 0.8)  # imag part
plt.xlabel('Point index')
plt.ylabel('Value')
plt.title('B-spline basis functions (real:black, imag:brown)')
plt.tight_layout()
plt.show()
plt.savefig("images/bsplines.png",dpi = 200)


plt.figure(figsize=(8,8))
for i in range(n_basis):
    plt.plot(x, dbdata[i,:,0], color='black',linewidth = 0.8)  # real part
    plt.plot(x, dbdata[i,:,1], color='brown',linewidth = 0.8)  # imag part
plt.xlabel('Point index')
plt.ylabel('Value')
plt.title('B-spline basis functions (real:black, imag:brown)')
plt.tight_layout()
plt.show()
plt.savefig("images/dbsplines.png",dpi = 200)


