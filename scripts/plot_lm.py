import json
import matplotlib.pyplot as plt
import numpy as np

with open("input.json") as f:
    data = json.load(f)

lmax = data["Angular"]["lmax"]


lmMap = np.loadtxt("misc/lm_map.txt",dtype = np.int32)
lvalues = lmMap[:,0]
mvalues = lmMap[:,1]
blockIndices = lmMap[:,2]

space_size =lmax + 1
space = np.zeros((space_size, 2 * lmax + 1))

for i in range(len(lvalues)):
    space[lmax - int(lvalues[i]), int(mvalues[i]) + lmax] = 1

fig,ax = plt.subplots()
ax.imshow(np.flipud(space), cmap='gray', interpolation='none', origin='lower')
ax.set_xlabel('m')
ax.set_ylabel('l')
ax.set_xticks([i for i in range(0, 2 * lmax + 1, 10)])  
ax.set_xticklabels([str(i - lmax) for i in range(0, 2 * lmax + 1, 10)])  
ax.set_title('Reachable (white)  Points in l-m Space')
fig.savefig("images/lm_space.png")