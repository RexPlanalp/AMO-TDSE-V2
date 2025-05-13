
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

import json


with open("input.json") as f:
    data = json.load(f)


lmax = data["Angular"]["lmax"]
lm_to_block_txt = np.loadtxt("misc/lm_map.txt")
probabilities_txt = np.loadtxt("misc/block_norms.txt")

fig,ax = plt.subplots(figsize=(10, 8))
space_size =lmax + 1
space = np.zeros((space_size, 2 * lmax + 1))

column1 = lm_to_block_txt[:,0]
column2 = lm_to_block_txt[:,1]
column3 = lm_to_block_txt[:,2]


column4 = probabilities_txt[:,0]
column5 = probabilities_txt[:,1]
for i in range(len(column1)):
    space[lmax - int(column1[i]), int(column2[i]) + lmax] = column5[i]

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
fig.savefig("images/pyramid.png",dpi = 200)
fig.clf()