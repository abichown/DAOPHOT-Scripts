'''
Purpose: Correct the instrumental magnitudes via zero point, standard aperture, aperture and location corrections
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

import numpy as np
import pandas as pd 
from astropy.io import fits
from math import log10 
import os
import sys

# Zero point correction - applied to both ap and alf mags
def zero_point(row_of_df, start_dither):

	# Image name - only need to do it once as they are all using same telescope etc
	image = stem + '_d1_cbcd_dn.fits'

	# Calculate zmag
	fits_image = fits.open(image)
	hdr = fits_image[0].header
	fluxconv = hdr['fluxconv'] * (10 ** 6) 
	px_ste = 3.3845 * (10 ** -11) # px size in steradians using online converters

	if channel == 1:
		F0 = 280.9
	elif channel == 2:
		F0 = 179.7
	else: F0 = 'Invalid'

	zmag = round(2.5 * log10(F0/(fluxconv*px_ste)), 2)

	# Load alf files, apply zmag and write out to file
	for alf in alf_files:

		# Load data
		id, x, y, m, e = np.loadtxt(alf, skiprows=4, usecols=(0,1,2,3,4), unpack=True)

		# Apply zmag to m
		m = m - 25 + zmag

		# Write to new file
		filename = alf.replace('.alf', '.alf_zp')
		f = open(filename, 'w')
		f.write("ID  X  Y  Mag  Error \n")

		for z in range(0, len(id)):
			f.writelines("%d %.3f %.3f %.4f %.4f \n" % (id[z], x[z], y[z], m[z], e[z]))

		f.close()

	# Load ap files, apply zmag and write out to file
	for ap in ap_files:

		# Load data - this is in df's because of the rubbish way .ap files are laid out
		num_lines = sum(1 for line in open(ap)) # num lines in the ap file

		first_rows = filter(lambda x: x%3 == 1, range(4, num_lines))
		second_rows = filter(lambda x: x%3 == 2, range(4,num_lines))
		third_rows = filter(lambda x: x%3 == 0, range(4, num_lines))

		ap_firstlines = pd.read_csv(ap, skiprows=[0,1,2,3]+second_rows+third_rows, names=['ID', 'X', 'Y', 'Mag'], delim_whitespace=True, header=None)
		ap_secondlines = pd.read_csv(ap, skiprows=[0,1,2,3]+first_rows+third_rows, names=['Sky', 'Std sky', 'Skew sky', 'Error'], delim_whitespace=True, header=None)

		# Concat and drop unneccesary columns
		ap_df = pd.concat((ap_firstlines, ap_secondlines), axis=1)
		ap_df.drop(['Sky', 'Std sky', 'Skew sky'], inplace=True, axis=1)

		# Apply zmag to 'Mag' column
		for z in range(0, len(ap_df)):
			if ap_df['Mag'][z]!= 99.999:
				ap_df.ix[z,'Mag'] = ap_df['Mag'][z] - 25 + zmag

				# Convert to fl

		# Write to new file
		filename = ap.replace('.ap', '.ap_zp')
		f = open(filename, 'w')
		f.write("ID  X  Y  Mag  Error \n")

		for z in range(0, len(ap_df)):
			f.writelines("%d %.3f %.3f %.4f %.4f \n" % (ap_df['ID'][z], ap_df['X'][z], ap_df['Y'][z], ap_df['Mag'][z], ap_df['Error'][z]))

		f.close()

	return(0)

# Get onto standard aperture system - applied to ap mags only
def std_aper():
	return(0)

# Aperture correction - uses ap and alf mags to calculate correction value then applies to psf mags only
def ap_corr():
	return(0)

# Location correction - only need to apply to alf mags as won't be using ap mags after this
def loc_corr():
	return(0)

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
				MAIN PART
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Read in list of stars that need correcting
df = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star', 'Channel', 'Epoch'])

for i in range(0, len(df)):
	for j in [1,6]: 

		# Initial setup
		galaxy = df['Galaxy'][i]
		target_name = df['Star'][i]

		if df['Channel'][i] == 1:
			channel = 1
			wavelength = '3p6um'
		elif df['Channel'][i] == 2:
			channel = 2
			wavelength = '4p5um'
		else: channel = 'Invalid channel'

		if df['Epoch'][i] < 10:
			epoch_number = '0' + str(df['Epoch'][i])
		else: epoch_number = str(df['Epoch'][i])

		if j == 1:
			field = 1
		elif j == 6:
			field = 2
		else: field = "Invalid field"

		# Find absolute path of where images are
		home = '/home/ac833/Data/'
		cwd = home + str(galaxy) + '/BCD/' + target_name + '/ch' + str(df['Channel'][i]) + '/e' + str(epoch_number) + '/'

		stem = target_name + '_' + wavelength + '_e' + epoch_number

		# Change working directory to where the data is
		os.chdir(cwd)

		# Set up names of ap and alf files - have 5 alf and 5 ap files for each epoch
		# Put them in two lists to make it easier to loop over later
		ap_files = []
		alf_files = []

		for dither in range(j, j+5):
			ap_files.append(stem + '_d' + str(dither) + '_cbcd_dn.ap')
			alf_files.append(stem + '_d' + str(dither) + '_cbcd_dn.alf')

		# Zero point
		zero_point(i,j)

		# Standard aperture
		#std_aper(i,j)

		# Aperture correction
		#ap_corr(i,j)

		# Location correction
		#loc_corr(i,j)
