import json
import matplotlib.pyplot as plt
import numpy as np

with open("input.json") as f:
    data = json.load(f)


lmax = data["Angular"]["l_max"]
lm_to_block_txt = np.loadtxt("misc/lm_map.txt")
fig,ax = plt.subplots()
space_size =lmax + 1
space = np.zeros((space_size, 2 * lmax + 1))

column1 = lm_to_block_txt[:,0]
column2 = lm_to_block_txt[:,1]
column3 = lm_to_block_txt[:,2]
for i in range(len(column1)):
    space[lmax - int(column1[i]), int(column2[i]) + lmax] = 1

ax.imshow(np.flipud(space), cmap='gray', interpolation='none', origin='lower')
ax.set_xlabel('m')
ax.set_ylabel('l')
ax.set_xticks([i for i in range(0, 2 * lmax + 1, 10)])  # Positions for ticks
ax.set_xticklabels([str(i - lmax) for i in range(0, 2 * lmax + 1, 10)])  # Labels from -lmax to lmax
ax.set_title('Reachable (white) and Unreachable (black) Points in l-m Space')
fig.savefig("images/lm_space.png")