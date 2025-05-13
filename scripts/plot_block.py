
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import json


# Set norm to either Log or Linear

with open("input.json") as f:
    data = json.load(f)

lmax = data["Angular"]["lmax"]

lmMap = np.loadtxt("misc/lm_map.txt")
blockNorms = np.loadtxt("misc/block_norms.txt")

fig,ax = plt.subplots(figsize=(10, 8))
space_size = lmax + 1
space = np.zeros((space_size, 2 * lmax + 1))

lvalues = lmMap[:,0]
mvalues = lmMap[:,1]

blockIndices = blockNorms[:,0]
probVals = blockNorms[:,1]
for i in range(len(lvalues)):
    space[lmax - int(lvalues[i]), int(mvalues[i]) + lmax] = probVals[i]

space[space==0] = np.min(space[space!=0])
cax = ax.imshow(space, cmap='inferno', interpolation='nearest', norm=Normalize())  
ax.set_xlabel('m')
ax.set_ylabel('l')
ax.set_xticks([0, lmax, 2 * lmax])  
ax.set_xticklabels([-lmax, 0, lmax])
ax.set_yticks([0, lmax])  
ax.set_yticklabels([lmax, 0])  
ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False) 
fig.colorbar(cax, ax=ax, shrink=0.5)
ax.set_title('Heatmap of Probabilities for l and m Values')
fig.savefig("images/block_norms.png",dpi = 200)
fig.clf()