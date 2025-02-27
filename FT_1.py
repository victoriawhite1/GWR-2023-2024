
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

#Load h_plus and h_cross from the saved files

h_plus = np.loadtxt('h_plus.txt')
h_cross = np.loadtxt('h_cross.txt')
f,advanced_LIGO = np.loadtxt('advanced_LIGO.txt',unpack = True) 
d = 3.09e22 #10 kpc in cm
#compute transform as needed:
#put h_plus and h_cross into a single array

h_combined = (h_plus + h_cross) /(np.sqrt(2)*d) 

#Compute the Fourier transform

fourier_transform = np.fft.rfft(h_combined)

#The Fourier transform will be a complex array
#compute the magnitude or phase if needed

magnitude = np.abs(fourier_transform) 
phase = np.angle(fourier_transform)

#you can also compute the frequencies associated with the Fourier transform
# The frequency bins will range from 0 to the Nyquist frequency
num_samples = len(h_combined)

sampling_rate = 16384 #Hertz
frequency_bins = np.fft.rfftfreq(num_samples, d=1/sampling_rate)
low_idx = np.squeeze(np.argwhere(frequency_bins >= 5) )[0]
high_idx = np.squeeze(np.argwhere(frequency_bins >= 5000))[0]
freq = frequency_bins[low_idx:high_idx]
magnitude2 = magnitude[low_idx:high_idx]
#resample the advance_LIGO
f_sf = interp1d(f,advanced_LIGO)
new_advanced_LIGO = f_sf(freq)

#calculate SNR from sigma_squared

T = 1 / (freq[2] - freq[1])#observation time 
sigma_squared = (T/2) * new_advanced_LIGO**2
SNR_squared = (2* np.sum(magnitude2**2 / sigma_squared))
print(np.sqrt(SNR_squared))


#create plots
myfig,ax = plt.subplots()
ax.loglog(frequency_bins,magnitude)
#ax.loglog(freq, np.sqrt(sigma_squared), 'g')
ax.loglog(freq, new_advanced_LIGO, 'm')
ax.set_xlabel("Frequency [Hz]") 
ax.set_ylabel("Magnitude of Fourier Transform")
plt.show() 
