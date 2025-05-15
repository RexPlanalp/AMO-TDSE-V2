import numpy as np
import matplotlib.pyplot as plt
import json


with open("input.json") as f:
    input = json.load(f)



polarization = input["Laser"]["polarization"]

I_max = input["Laser"]["I"]/3.51E16
w = input["Laser"]["w"]
N = input["Laser"]["N"]

if input["Atom"]["potential"] == "H":
    Ip = -0.5
elif input["Atom"]["potential"] == "Ar":
    Ip = -0.5791546178
elif input["Atom"]["potential"] == "He":
    Ip = -0.9443

Up = I_max/(4*w**2)

cut_off = 3.17*Up - Ip


hhgOutput = np.loadtxt("misc/hhg_data.txt")



# Load data
hhgData = hhgOutput[:, [1, 3, 5]].T
laserData = hhgOutput[:, [2,4,6]].T
totalTime = hhgOutput[:, 0]

# Initialize dipole_acceleration as a 2D array to store each component
dipole_acceleration = np.zeros((3, len(totalTime)))

# Calculate dipole acceleration for each component
for i in range(3):
    dipole_acceleration[i, :] = -hhgData[i, :] + np.gradient(laserData[i, :], totalTime)

# Plot the dipole acceleration components (optional)
plt.figure(figsize=(8, 6))
for i in range(3):
    plt.plot(totalTime, dipole_acceleration[i, :], label=f'Component {i+1}')
plt.xlabel('Time (a.u.)')
plt.ylabel('Dipole Acceleration (a.u.)')
plt.title('Dipole Acceleration Components')
plt.legend()
plt.savefig("images/dipole_accel.png")
plt.clf()

# Apply a window function (Blackman window) to each component
window = np.blackman(len(totalTime))
windowed_dipole = dipole_acceleration * window  # Broadcasting applies the window to each row

# Perform FFT on each component
dipole_fft = np.fft.fft(windowed_dipole, axis=1)
frequencies = np.fft.fftfreq(len(totalTime), d=totalTime[1] - totalTime[0])

# Compute the power spectrum (magnitude squared of FFT) for each component
power_spectrum = np.abs(dipole_fft)**2

# Sum the power spectra of all components to get the total power spectrum
total_power_spectrum = np.sum(power_spectrum, axis=0)

# Only keep the positive frequencies
positive_freq_idx = frequencies > 0
frequencies = frequencies[positive_freq_idx]
total_power_spectrum = total_power_spectrum[positive_freq_idx]

# Plot the total harmonic spectrum
plt.figure(figsize=(8, 6))
plt.semilogy(frequencies * 2*np.pi / w, total_power_spectrum, color='b')
plt.axvline(cut_off/w, color='r', linestyle='--', label='Cut-off Energy')
plt.xlim([0, 60])        # Adjust based on your data's frequency range
plt.ylim([1e-4, 1e4])     # Adjust based on the power spectrum's range
plt.xlabel('Frequency (atomic units)')
plt.ylabel('Intensity (arb. units)')
plt.title('Harmonic Spectrum')
plt.grid(True, which='both', ls='--', lw=0.5)
plt.legend()
plt.savefig("images/harmonic_spectrum.png")
plt.clf()
