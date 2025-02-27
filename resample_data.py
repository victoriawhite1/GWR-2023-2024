import numpy as np
from scipy.interpolate import interp1d
import os

def resample_and_save_data_for_all_masses():
    mass_names = ['9a', '9b', '9.25', '9.5', '11', '15.01', '23']
    base_path = "/Users/victoriawhite/Downloads/quad_2023"  #base path to the quad_2023 directory
    for mass_name in mass_names:
        resample_and_save_data_for_mass(base_path, mass_name)

def resample_and_save_data_for_mass(base_path, mass_name):
    mass_folder_path = os.path.join(base_path, mass_name)
    quadrupole_data_path = os.path.join(mass_folder_path, "quadrupole.dat")
    t, qxx, qxy, qxz, qyy, qyz, qzz, _, _, _, _, _, _ = np.loadtxt(quadrupole_data_path, unpack=True)

    # Resample data to 16384 Hz
    new_t = np.linspace(t.min(), t.max(), int((t.max() - t.min()) * 16384))
    f_qxx = interp1d(t, qxx)
    f_qxy = interp1d(t, qxy)
    f_qxz = interp1d(t, qxz)
    f_qyy = interp1d(t, qyy)
    f_qyz = interp1d(t, qyz)
    f_qzz = interp1d(t, qzz)
    new_qxx = f_qxx(new_t)
    new_qxy = f_qxy(new_t)
    new_qxz = f_qxz(new_t)
    new_qyy = f_qyy(new_t)
    new_qyz = f_qyz(new_t)
    new_qzz = f_qzz(new_t)

    # Save resampled data
    resampled_data = np.column_stack((new_t, new_qxx, new_qxy, new_qxz, new_qyy, new_qyz, new_qzz))
    output_file_path = os.path.join(base_path, mass_name, "quadrupole_resampled.dat")
    np.savetxt(output_file_path, resampled_data, fmt='%1.16e', header="t qxx qxy qxz qyy qyz qzz")

if __name__ == "__main__":
    resample_and_save_data_for_all_masses()






