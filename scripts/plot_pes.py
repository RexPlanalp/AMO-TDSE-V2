import numpy as np 

import matplotlib.pyplot as plt 
from matplotlib.colors import LogNorm, Normalize 
import matplotlib.colors as mcolors 

import json
import sys

with open("input.json") as f:
    input = json.load(f)


slice = input["Observables"]["PES"]["slice"]
norm = input["Observables"]["PES"]["norm"]
threshold = input["Observables"]["PES"]["threshold"]

# Plot Angle Integrated Spectrum
angleIntegratedData = np.loadtxt("misc/pes.txt")

energyData = np.array(angleIntegratedData[:,0])
pesData = np.array(angleIntegratedData[:,1])

fig,ax = plt.subplots()
ax.semilogy(energyData,pesData,label = "Angle Integrated Spectrum",color = "k",linewidth = 0.8)
ax.set_ylabel("Probability")
ax.set_xlabel("Energy (au)")
ax.legend()
fig.savefig("images/pes.png",dpi = 200)


# Plot Angle Resolved Spectrum
angleResolvedData = np.loadtxt("misc/pad.txt")

energyData = np.array(angleResolvedData[:,0])
momentumData = np.sqrt(2*energyData)

thetaData = np.array(angleResolvedData[:,1])
phiData = np.array(angleResolvedData[:,2])
padData = np.array(angleResolvedData[:,3])

kx = momentumData*np.sin(thetaData)*np.cos(phiData)
ky = momentumData*np.sin(thetaData)*np.sin(phiData)
kz = momentumData*np.cos(thetaData)

maxVal = np.max(padData)
minVal= maxVal*threshold

cmap = "hot_r"

fig,ax = plt.subplots()

if (slice == "XZ"):
    xLabel = "kx"
    yLabel = "kz"

    xData = kx
    yData = kz

elif (slice == "XY"):
    xLabel = "kx"
    yLabel = "ky"

    xData = kx
    yData = ky

if (norm == "log"):
    plotNorm = mcolors.LogNorm(minVal,maxVal)
elif (norm == "linear"):
    plotNorm = mcolors.Normalize(minVal,maxVal)

scatterPlot = ax.scatter(xData,yData,c=padData,norm=plotNorm,cmap=cmap)
ax.set_aspect("equal",adjustable = "box")
ax.set_xlabel(xLabel)
ax.set_ylabel(yLabel)
fig.colorbar(scatterPlot,ax=ax)
fig.savefig("images/pad.png",dpi = 200)
