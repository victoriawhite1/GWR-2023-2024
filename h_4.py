import numpy as np
import os
import matplotlib.pyplot as plt

def load_quadrupole_data(zip_path, model_name):
    """
    Load quadrupole data from the specified model name in the zip file.
    """
    model_folder = os.path.join(zip_path, model_name, '11')  #using mass 23 as example
    data_file = os.path.join(model_folder, 'quadrupole.dat')
    
    quadrupole_data = np.loadtxt(data_file)
    return quadrupole_data

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

#

def main():
    # Set the path to the quad_2023 directory
    zip_path = '/Users/victoriawhite/Downloads'
    model_name = 'quad_2023'  #Change this to the desired model folder

    # Load quadrupole data for the specified model
    quadrupole_data = load_quadrupole_data(zip_path, model_name)
    nt = len(quadrupole_data)-1  #how long quadrupole data array  is
 
    # Extract time-changing quadrupole moments and other parameters
    time = quadrupole_data[:, 0]
    q = quadrupole_data[:, 1:]
    theta = np.pi/2  #(in radians) inclination angle
    phi = 0.0    #(radians) azimuthal angle
    D = 3.09e22    #source distance in centimeters of 10kpc
    
    # Convert quadrupole moments to spherical coordinates
    q_theta_theta, q_phi_phi, q_theta_phi = convert_to_spherical(q, theta, phi)
    
    # Calculate time step
    dt = np.diff(time)
    
    # Calculate strain components
    h_plus, h_cross = calculate_strain(q_theta_theta, q_phi_phi, q_theta_phi, dt, D)
    
    # Print or return the strain components
    print("h_plus:", h_plus)
    print("h_cross:", h_cross)

    #create array with 3 columns
    print(len(quadrupole_data))
    outarr =  np.zeros((nt,3))
    outarr[:,0] = time[1:]
    outarr[:,1] = h_plus
    outarr[:,2] = h_cross
    outfile = "gw_strain_model23_theta%1.3f_phi%1.3f_D%1.3e.txt"%(theta,phi,D)  
    np.savetxt(outfile,outarr) 

if __name__ == "__main__":
    main()


def compare_angles():
    #compare the outfiles for theta as zero and pi/2
    time,h_plus,h_cross = np.loadtxt("gw_strain_model23_theta1.571_phi0.000_D1.000e+00.txt", unpack=True)  #fill 3 arrays with information from the file
    return time, h_plus, h_cross

time, h_plus, h_cross = compare_angles()
# Save h_plus and h_cross to a file
np.savetxt('h_plus.txt', h_plus)
np.savetxt('h_cross.txt', h_cross)
 #create a figure-object with axes
myfig,ax = plt.subplots() 
ax.plot(time,h_plus)
ax.plot(time,h_cross) 
plt.show()



