import numpy as np
from scipy.interpolate import interp1d

import os
import matplotlib.pyplot as plt


zip_path = '/Users/victoriawhite/Downloads/quad_2023'
#load quadrupole data
def load_quadrupole_data(zip_path, model_name):
    """
    Load quadrupole data from the specified model name in the zip file.
    """
    model_folder = os.path.join(zip_path, model_name)  #using mass 15 as example
    data_file = os.path.join(model_folder, 'quadrupole.dat')

    quadrupole_data = np.loadtxt(data_file)
    return quadrupole_data

#calculate GW strain by finding h+ and hx
def calculate_strain(q_theta_theta, q_phi_phi, q_theta_phi, dt, D):
    """
    Calculate the strain components h_+ and h_x.
    """
    G = 6.67430e-8  # Gravitational constant
    c = 2.99792e10    # Speed of light
    U_MSUN =  1.9891e33
    # Calculate h_+
    h_plus = (G*U_MSUN / (c**4 * D)) * (np.diff(q_theta_theta) / dt - np.diff(q_phi_phi) / dt)

    # Calculate h_cross
    h_cross = 2 * (G*U_MSUN / (c**4 * D)) * np.diff(q_theta_phi) / dt

    return h_plus, h_cross

#convert q, theta and phi to spherical coordinates
def convert_to_spherical(q, theta, phi):
    """
    Convert quadrupole moment from Cartesian to spherical coordinates.
    """
    q_theta_theta = (q[:, 6] * np.cos(phi) ** 2 +
                    q[:, 9] * np.sin(phi) ** 2 +
                    2 * q[:, 7] * np.sin(phi) * np.cos(phi)) * np.cos(theta) ** 2 + \
                    q[:, 11] * np.sin(theta) ** 2 - \
                    2 * (q[:, 8] * np.cos(phi) + q[:, 10] * np.sin(phi)) * np.sin(theta) * np.cos(theta)

    q_phi_phi = (q[:, 6] * np.sin(phi) ** 2 +
                 q[:, 9] * np.cos(phi) ** 2 - \
                 2 * q[:, 7] * np.sin(phi) * np.cos(phi))

    q_theta_phi = (q[:, 9] - q[:, 6]) * np.cos(theta) * np.sin(phi) * np.cos(phi) + \
                  q[:, 7] * np.cos(theta)*(np.cos(phi) ** 2 - np.sin(phi) ** 2)  + \
                  q[:, 8] * np.sin(theta) * np.sin(phi) - q[:, 10] * np.sin(theta) * np.cos(phi)

    return q_theta_theta, q_phi_phi, q_theta_phi

def generate_SNR(q, time,D, theta, phi):
    	

#Convert quadrupole moments to spherical coordinates
    q_theta_theta, q_phi_phi, q_theta_phi = convert_to_spherical(q, theta, phi)

    #Calculate time step
    dt = np.diff(time)

    #Calculate strain components
    h_plus, h_cross = calculate_strain(q_theta_theta, q_phi_phi, q_theta_phi, dt, D)

    f, advanced_LIGO = np.loadtxt('advanced_LIGO.txt', unpack=True)
    h_combined = (h_plus + h_cross) / (np.sqrt(2))

    #Compute the Fourier transform
    fourier_transform = np.fft.rfft(h_combined)

    #The Fourier transform will be a complex array
    #Compute the magnitude and phase
    magnitude = np.abs(fourier_transform)
    h_f = magnitude * (1 / 16384)
    phase = np.angle(fourier_transform)

    #Compute the frequencies associated with the Fourier transform
    #The frequency bins will range from 0 to the Nyquist frequency
    num_samples = len(h_combined)
    sampling_rate = 16384  # Hertz
    frequency_bins = np.fft.rfftfreq(num_samples, d=1 / sampling_rate)
    low_idx = np.squeeze(np.argwhere(frequency_bins >= 5))[0]
    high_idx = np.squeeze(np.argwhere(frequency_bins >= 5000))[0]
    freq = frequency_bins[low_idx:high_idx]
    magnitude2 = magnitude[low_idx:high_idx]

    # Resample the advanced_LIGO
    f_sf = interp1d(f, advanced_LIGO)
    new_advanced_LIGO = f_sf(freq)

    # Calculate SNR from sigma_squared
    T = 1 / (freq[2] - freq[1])  # Observation time
    sigma_squared = (T / 2) * new_advanced_LIGO ** 2
    SNR_sq = 4 * np.sum((h_f[low_idx:high_idx] ** 2) / new_advanced_LIGO ** 2) / T
    SNR_squared = (2 * np.sum(magnitude2 ** 2 / sigma_squared))
    SNR1 = np.sqrt(SNR_sq)
    return SNR1

def generate_source_orientations(nphi=1000, ntheta=1000):
    #Generate random distribution of sources across the sky.
    phi_array = np.random.uniform(0, 2 * np.pi, nphi)
    costheta_array = np.random.uniform(-1, 1, ntheta)
    theta_array = np.arccos(costheta_array)
    return theta_array, phi_array

#define angle_avegrage_SNR so you can call it within the main function
def angle_avg_snr(model_name, D):
    """#add function description"""
     # Load quadrupole data for the specified model
    quadrupole_data = load_quadrupole_data(zip_path, model_name)
    nt = len(quadrupole_data) - 1 #how long quadrupole data array is
# Extract time-changing quadrupole moments and other parameters
    time = quadrupole_data[:, 0]
    q = quadrupole_data[:, 1:]
   

    
    # Calculate angle average SNR
    theta_array, phi_array = generate_source_orientations(nphi=250, ntheta=250)


    # Store SNR values
    snr_values = []


    # Loop through the orientations to compute SNR averaged over all orientations
    for theta, phi in zip(theta_array, phi_array):
        SNR = generate_SNR(q, time, D, theta, phi)


        print("theta=", theta, "phi=", phi, "SNR:", SNR)
        snr_values.append(SNR)

    avg_snr = np.mean(snr_values)



    return avg_snr



def main():
   
    model_names = ['9.5', '11', '15.01', '23']


    D = 3.09e22  # source distance in centimeters of 10kpc

    for model in model_names: 
        avg_snr = angle_avg_snr(model, D)



        print("Avg SNR over source orientation for model %s: %1.4f"%(model, avg_snr))
    plotflag = False
    if plotflag:
        # create figure
        # myfig,ax = plt.subplots(projection='3d')
        myfig = plt.figure()
        ax1 = myfig.add_subplot(2, 2, 1,projection='3d')
        plot1=ax1.scatter(np.sin(theta_array) * np.cos(phi_array), np.sin(theta_array) * np.sin(phi_array), np.cos(theta_array),c=snr_values1)
        cbar = plt.colorbar(plot1)
        cbar.set_label("SNR",labelpad=10)
 
        ax2 = myfig.add_subplot(2, 2, 2,projection='3d')
        plot2=ax2.scatter(np.sin(theta_array) * np.cos(phi_array), np.sin(theta_array) * np.sin(phi_array), np.cos(theta_array),c=snr_values2)
        cbar2 = plt.colorbar(plot2)
        cbar2.set_label("SNR",labelpad=10)
    
        ax3 = myfig.add_subplot(2, 2, 3,projection='3d')
        plot3=ax3.scatter(np.sin(theta_array) * np.cos(phi_array), np.sin(theta_array) * np.sin(phi_array), np.cos(theta_array),c=snr_values3)
        cbar3 = plt.colorbar(plot3)
        cbar3.set_label("SNR",labelpad=10)
      
        ax4 = myfig.add_subplot(2, 2, 4,projection='3d')
        plot4=ax4.scatter(np.sin(theta_array) * np.cos(phi_array), np.sin(theta_array) * np.sin(phi_array), np.cos(theta_array),c=snr_values4)
        cbar4 = plt.colorbar(plot4)
        cbar4.set_label("SNR",labelpad=10)
     
    
    plt.show()
    # Plot
 #  myfig, ax = plt.subplots()
   #ax.loglog(frequency_bins, 2 * np.sqrt(frequency_bins) * magnitude * (1 / 16384))
   #ax.loglog(freq, new_advanced_LIGO, 'm')
   #ax.set_ylabel("Magnitude of Fourier Transform")
    #plt.show()

    # Print or return the strain components
    # create array with 3 columns
    # outarr = np.zeros((nt,3))
    # outarr[:,0] = time[1:]
    # outarr[:,1] = h_plus
    # outarr[:,2] = h_cross
    # outfile = "gw_strain_model23_theta%1.3f_phi%1.3f_D%1.3e.txt"%(theta,phi,D)
    # np.savetxt(outfile,outarr)

# source_dist.py

# use theta and phi to get random distribution of sources across the sky
# theta ranges 0 to pi
# phi ranges 0 to 2pi

# created theta and phi arrays earlier
# phi_array = np.random.uniform (0, 2*np.pi, 1000)
# print(phi_array)

# print(np.mean(phi_array))


#costheta_array = np.random.uniform (-1, 1, 1000)
# print(costheta_array)
# print(np.mean(costheta_array))

# theta_array = np.arccos(costheta_array)
# print(theta_array)
# print(np.mean(theta_array))

# Generate source orientations before running main function
theta_array, phi_array = generate_source_orientations()

# create figure
# myfig,ax = plt.subplots(projection='3d')
#myfig = plt.figure()
#ax = myfig.add_subplot(projection='3d')

#ax.scatter(np.sin(theta_array) * np.cos(phi_array), np.sin(theta_array) * np.sin(phi_array), np.cos(theta_array))
#plt.show()

# Generate source orientations before running main function
#theta_array, phi_array = generate_source_orientations()

if __name__ == "__main__":
    main()



