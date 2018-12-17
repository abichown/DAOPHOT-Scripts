'''
Purpose: Calibrate the PSF instrumental magnitudes from ALLFRAME to put onto the standard IRAC Vega system (Reach 2005).
1. Location correction - corrects for effects of IRAC
2. Aperture correction - converts PSF magnitudes onto (3,12,20) aperture system
3. Standard aperture calibration - converts (3,12,20) onto (10,12,20) standard system
4. Zero magnitude correction - corrects for zero magnitude point
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

import numpy as np
import pandas as pd 
from astropy.io import fits
from math import log10 
import os
import sys
import pexpect



def ap_corr(row_of_df, start_dither):

	# Applies to '.alf' magnitudes only
	for alf in alf_files:

		# Change filenames - need to add the _loc because the files have already been corrected for location
		#alf = alf.replace('.alf', '.alf_loc')
		als = alf.replace('.alf', '.als')
		ap = alf.replace('.alf', '.ap')

		# Image name
		image = alf.replace('.alf', '.fits')

		# lst file name
		lst = alf.replace('.alf', '.lst')


		# Open files into dataframes
		# The lst dataframe contains the aperture magnitudes and ALL stars will be used from this file
		# The als dataframe contains the psf magnitudes and only a SUBSET of these stars with IDs matching those of the lst df will be used

		df_lst = pd.read_csv(lst, skiprows=3, header=None, names=['ID_lst', 'X_lst', 'Y_lst', 'ap_mag', 'ap_error'], usecols=(0,1,2,3,4), delim_whitespace=True)
		df_als = pd.read_csv(als, skiprows=3, header=None, names=['ID_als', 'X_als', 'Y_als', 'psf_mag', 'psf_error'], usecols=(0,1,2,3,4), delim_whitespace=True)

		pd.set_option('display.max_rows', 500)

		# Sort rows by ID number
		df_lst.sort_values('ID_lst', axis=0, inplace=True)
		df_als.sort_values('ID_als', axis=0, inplace=True)

		# Reset index to start at 0
		df_lst.reset_index(drop=True, inplace=True)
		df_als.reset_index(drop=True, inplace=True)


		# Cross-check the two dataframes to only keep stars that appear in both dataframes
		for i in range(0, len(df_als)):

			match = 0 

			for j in range(0, len(df_lst)):

				# Check whether IDs match
				if df_als['ID_als'][i] == df_lst['ID_lst'][j]:
					match += 1

			# If they match, do nothing as want to keep this row. Otherwise, delete the row in the als dataframe
			if match == 0:
				df_als.drop(i, axis=0, inplace=True)

		# Reset index to start at 0
		df_lst.reset_index(drop=True, inplace=True)
		df_als.reset_index(drop=True, inplace=True)

		# Now do the same process but removing any stars in the lst file not in the als file as the stars need to have magnitudes in both dataframes
		for i in range(0, len(df_lst)):

			match = 0

			for j in range(0, len(df_als)):

				# Check whether IDs match
				if df_lst['ID_lst'][i] == df_als['ID_als'][j]:
					match += 1

			# If they match, do nothing as want to keep this row. Otherwise, delete the row in the als dataframe
			if match == 0:
				df_lst.drop(i, axis=0, inplace=True)

		# Again, reset index to start at 0
		df_lst.reset_index(drop=True, inplace=True)
		df_als.reset_index(drop=True, inplace=True)	

		# Now df_als and df_lst have the same stars so can concatenate
		df_combined = pd.concat((df_als, df_lst), axis=1)

		# Remove columns not needed
		df_combined.drop(['X_lst', 'Y_lst', 'ID_lst'], axis=1, inplace=True)

		# Now create new column in df_combined to calculate the difference between aperture and psf magnitudes
		df_combined['mag_difference'] = df_combined['ap_mag'] - df_combined['psf_mag']	

		# Calculate aperture correction value 
		apc = round(df_combined['mag_difference'].mean(), 3)

		# # Open ALLFRAME magnitude file (.alf_loc) and apply this apc value to these magnitudes
		# # This converts the PSF magnitudes from ALLFRAME onto the (3,12,20) aperture magnitude system
		# # The standard aperture calibration step is then required to convert these (3,12,20) magnitudes onto the standard IRAC Vega (10,12,20) system

		df_alf = pd.read_csv(alf, header=None, skiprows=3, usecols=[0,1,2,3,4], names=['ID', 'x', 'y', 'mag', 'error'], delim_whitespace=True)

		df_alf['mag'] += apc # ADD apc value to magnitudes

		# # Write out new corrected magnitudes to file
		new_filename = alf.replace('.alf', '.alf_apc')
		df_alf.to_csv(new_filename, sep=' ', index=False)

	return(0)

def zp_and_std_ape(row_of_df, start_dither):

	for alf in alf_files:

		# Working on the .alf_apc magnitudes which have been location and aperture corrected
		alf = alf.replace('.alf', '.alf_apc')

		# Calculate zmag value
		image = alf.replace('.alf_apc', '.fits')
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

		# To go from (3,12,20) to (10,12,20) we have a correction value from IRAC handbook
		# For this aperture and annulus size, it is the same value for [3.6] and [4.5]
		std_corr = 1.112

		# Load alf_apc magnitudes into a dataframe
		df_alf = pd.read_csv(alf, header=0, delim_whitespace=True)

		for z in range(0, len(df_alf)):

			# Convert magnitude to flux
			flux =  F0 * (10 ** (-df_alf['mag'][z]/2.5))

			# Apply std aperture correction value
			flux *= std_corr

			# Convert back to magnitude
			df_alf.ix[z, 'mag'] = -2.5 * log10(flux/F0)

		# Apply zmag to magnitudes
		df_alf['mag'] = df_alf['mag'] + zmag - 25

		# Round off magnitudes to 4 dp
		for z in range(0, len(df_alf)):

			df_alf['mag'][z] = round(df_alf['mag'][z], 4)

		# Write out to file - 'cal' for calibrated
		new_filename = alf.replace('.alf_apc', '.alf_zp')
		df_alf.to_csv(new_filename, index=False, sep=' ')

	return(0)

def loc_corr(row_of_df, start_dither):

	# Open correction file corresponding to correct channel
	if channel == 1:
		corr_file = '/home/ac833/Spitzer-corr-images/ch1_al_s192.fits'
		corr_image = fits.open(corr_file)
		corr_data = corr_image[0].data
	elif channel == 2:
		corr_file = '/home/ac833/Spitzer-corr-images/ch2_al_s192.fits'
		corr_image = fits.open(corr_file)
		corr_data = corr_image[0].data
	else: corr_file = 'Corr file not found'

	# Change of file names
	# Not required here as this is the first step of the calibration procedure

	# Loop over all '.alf' files and correct
	for alf in alf_files:	

		alf = alf.replace('.alf', '.alf_zp')

		# File that needs correcting is the alf file as this is the first correction
		# Load data into df
		df = pd.read_csv(alf, header=0, delim_whitespace=True)	

		# Convert to flux, correct and convert back to mag for each star
		for index in range(0, len(df)):

			# Find the x and y coordinates corresponding to that star
			x_coord = int(np.floor(float(df['x'][index])))
			y_coord = int(np.floor(float(df['y'][index])))

			# Checks to reset edge of frame stars
			if x_coord >= 256:
				x_coord = 256
			if x_coord <= 1:
				x_coord = 1
			if y_coord >= 256:
				y_coord = 256
			if y_coord <= 1:
				y_coord = 1

			# Find the corr value at that (x_coord, y_coord) pixel
			corr_val = corr_data[x_coord-1, y_coord-1]

			# Convert to flux, apply corr val, and convert back to mag
			if df['mag'][index] != 99.999:
				flux = 10 ** (df['mag'][index]/-2.5) # convert to flux
				flux = flux * corr_val # apply correction
				df.loc[index, 'mag'] = -2.5 * log10(flux) # convert back to mag and update df

		# Make new filename by replacing .alf with .alf_loc
		alf = alf.replace('.alf_zp', '.alf_cal')

		# Write csv to file
		df.to_csv(alf, sep=' ', index=False)

	# # Loop over all '.als' files and correct
	# for als in als_files:	

	# 	# File that needs correcting is the alf file as this is the first correction
	# 	# Load data into df
	# 	df = pd.read_csv(als, header=None, skiprows=3, delim_whitespace=True, usecols=(0,1,2,3,4), names=['ID', 'X', 'Y', 'mag', 'error'])	

	# 	# Convert to flux, correct and convert back to mag for each star
	# 	for index in range(0, len(df)):

	# 		# Find the x and y coordinates corresponding to that star
	# 		x_coord = int(np.floor(float(df['X'][index])))
	# 		y_coord = int(np.floor(float(df['Y'][index])))

	# 		# Find the corr value at that (x_coord, y_coord) pixel
	# 		corr_val = corr_data[x_coord-1, y_coord-1]

	# 		# Convert to flux, apply corr val, and convert back to mag
	# 		if df['mag'][index] != 99.999:
	# 			flux = 10 ** (df['mag'][index]/-2.5) # convert to flux
	# 			flux = flux * corr_val # apply correction
	# 			df.loc[index, 'mag'] = -2.5 * log10(flux) # convert back to mag and update df

	# 	# Make new filename by replacing .alf with .alf_loc
	# 	als = als.replace('.als', '.als_loc')

	# 	# Write csv to file
	# 	df.to_csv(als, sep=' ', index=False)		

	# # Loop over all aperture files as these also need to be corrected
	# for ap in ap_files:	

	# 	# Load data - this is in multiple df's because of the rubbish way .ap files are formatted
	# 	num_lines = sum(1 for line in open(ap)) # num lines in the ap file

	# 	first_rows = filter(lambda x: x%3 == 1, range(4, num_lines))
	# 	second_rows = filter(lambda x: x%3 == 2, range(4,num_lines))
	# 	third_rows = filter(lambda x: x%3 == 0, range(4, num_lines))

	# 	ap_firstlines = pd.read_csv(ap, skiprows=[0,1,2,3]+second_rows+third_rows, names=['ID', 'X', 'Y', 'mag'], delim_whitespace=True, header=None)
	# 	ap_secondlines = pd.read_csv(ap, skiprows=[0,1,2,3]+first_rows+third_rows, names=['Sky', 'Std sky', 'Skew sky', 'error'], delim_whitespace=True, header=None)

	# 	# Concat and drop unneccesary columns
	# 	ap_df = pd.concat((ap_firstlines, ap_secondlines), axis=1)
	# 	ap_df.drop(['Sky', 'Std sky', 'Skew sky'], inplace=True, axis=1)

	# 	# Convert to flux, correct and convert back to mag for each star
	# 	for index in range(0, len(ap_df)):

	# 		# Find the x and y coordinates corresponding to that star
	# 		x_coord = int(np.floor(float(ap_df['X'][index])))
	# 		y_coord = int(np.floor(float(ap_df['Y'][index])))

	# 		corr_val = corr_data[x_coord-1, y_coord-1]

	# 		if ap_df['mag'][index] != 99.999:
	# 			flux = 10 ** (ap_df['mag'][index]/-2.5) 
	# 			flux = flux * corr_val
	# 			ap_df.loc[index, 'mag'] = -2.5 * log10(flux)
		
	# 	# Make new filename by replacing _apc with _all
	# 	ap = ap.replace('.ap', '.ap_loc')

	# 	# Write csv to file
	# 	ap_df.to_csv(ap, sep=' ', index=False)	

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
				MAIN PART
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Read in list of stars that need correcting
df = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star', 'Period', 'RA', 'Dec', 'Channel'])

for i in range(0, len(df)):
	for j in [1,6]: #[1,6] when got both fields working 

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

		if j == 1:
			field = 1
		elif j == 6:
			field = 2
		else: field = "Invalid field"

	   	# Get number of epochs
		if galaxy == 'LMC':
	   		num_epochs = 24
		elif galaxy == 'SMC':
	   		num_epochs = 12
		else: num_epochs = 0

	   	# Run over each epoch
		for epoch in range(1, num_epochs+1):

			if epoch < 10:
				epoch_number = '0' + str(epoch)
			else: epoch_number = str(epoch)

			# Find absolute path of where images are
			home = '/home/ac833/Data/'
			cwd = home + str(galaxy) + '/BCD/' + target_name + '/ch' + str(df['Channel'][i]) + '/e' + str(epoch_number) + '/'

			stem = target_name + '_' + wavelength + '_e' + epoch_number

			# Change working directory to where the data is
			os.chdir(cwd)

			print "Working on: " + str(target_name) + "    Epoch: " + str(epoch_number) + "     Field: " + str(field)

			# Set up lists of ap, alf and als files to loop over later
			ap_files = []
			alf_files = []
			als_files = []

			for dither in range(j, j+5):
				ap_files.append(stem + '_d' + str(dither) + '_cbcd_dn.ap')
				alf_files.append(stem + '_d' + str(dither) + '_cbcd_dn.alf')
				als_files.append(stem + '_d' + str(dither) + '_cbcd_dn.als') 

			# Aperture correction
			ap_corr(i,j)

			# Standard aperture calibration and zero magnitude point correction
			zp_and_std_ape(i,j)

			# Location correction
			loc_corr(i,j)