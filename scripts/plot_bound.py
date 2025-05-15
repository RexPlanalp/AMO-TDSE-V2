import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm, Normalize
import json

with open("input.json") as f:
    input = json.load(f)

nmax = input["TISE"]["nmax"]
initialState = input["TDSE"]["initialNLM"]

bound_pops = np.loadtxt("misc/bound_pops.txt")

rows,cols = bound_pops.shape

pop_matrix = np.full((nmax+1, nmax+1), np.nan)

for i in range(rows):
    n = int(bound_pops[i,0])
    l = int(bound_pops[i,1])
    pop = bound_pops[i,2]

    if (n == initialState[0] and l == initialState[1]):
        continue
   
    pop_matrix[n,l] = pop
for n in range(nmax+1):
    for l in range(n):  # Only allow l < n
        if np.isnan(pop_matrix[n, l]):  # Ensure only invalid positions are masked
            pop_matrix[n, l] = None  
fig, ax = plt.subplots(figsize=(8, 6))
cmap = cm.viridis
cmap.set_bad(color='white')  # Make masked values white

im = ax.imshow(pop_matrix, cmap=cmap, origin="lower", norm=Normalize(), aspect='auto')

# Add colorbar
cbar = plt.colorbar(im)
cbar.set_label("Population")

# Axis labels and formatting
ax.set_xlabel("l (Orbital Quantum Number)")
ax.set_ylabel("n (Principle Angular Momentum)")
ax.set_title("Population Heatmap (Bound States)")
ax.set_xticks(range(nmax+1))
ax.set_yticks(range(nmax+1))

plt.savefig("images/bound_pops.png")