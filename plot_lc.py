'''
Purpose: Plot light curve from the corresponding file in the Magnitudes directory
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np

# Request details of what is to be plotted
galaxy = input("Galaxy? ")
star_name = input("Name of star: ")
wavelength = input("3p6um or 4p5um? ")
period = float(input("Period of star: "))

mag_file = '/home/ac833/Magnitudes/' + str(galaxy.upper()) + '/' + str(star_name.upper()) + '_' + wavelength + '.txt'

df = pd.read_csv(mag_file, header=0, delim_whitespace=True)

# List of all HMJDs for the observations
hmjd = []

# Get HMJDs
for i in range(1,25):
    
    if i < 10:
        filename = "/home/ac833/Data/LMC/BCD/HV00872/ch1/e0" + str(i) + "/HV00872_3p6um_e0" + str(i) + "_d1_cbcd_dn.fits"
    else: filename = "/home/ac833/Data/LMC/BCD/HV00872/ch1/e" + str(i) + "/HV00872_3p6um_e" + str(i) + "_d1_cbcd_dn.fits"
        
    image = fits.open(filename)
    header = image[0].header
    
    hmjd.append(header['HMJD_OBS'])


# Convert to numpy arrays
hmjd_array = np.array(hmjd)

# Get phase
phase = hmjd_array/period - np.floor(hmjd_array/period)

# Get all relevant info from df
mag = pd.concat([df['Mag'],df['Mag'], df['Mag']], axis=0)
err = pd.concat([df['Error'], df['Error'], df['Error']], axis=0)
std_dev = pd.concat([df['Std_dev'], df['Std_dev'], df['Std_dev']], axis=0)

# Negative phase
neg_phase = phase - 1

# Positive phase
pos_phase = phase + 1

# Concat phases
full_phase = np.concatenate((neg_phase, phase, pos_phase), axis=0)

# Plot
fig = plt.figure(figsize=[13,4])

axes = fig.add_axes([0.1,0.1,0.9,0.9])
axes.scatter(full_phase, mag, alpha=0.5, c='purple', linestyle='None')
axes.errorbar(full_phase, mag, yerr=std_dev, alpha=0.2, c='purple', linestyle='None')
axes.set_xlabel('Phase')

if wavelength == '3p6um':
	y_axis = '[3.6]'
else: y_axis = '[4.5]'
axes.set_ylabel(y_axis)
axes.set_title(star_name.upper())
axes.invert_yaxis()

# Save file
new_file = '/home/ac833/Light_Curves/' + galaxy.upper() + '/' + star_name.upper() + '_' + wavelength + '.png'
print new_file
plt.savefig(new_file)

plt.show()