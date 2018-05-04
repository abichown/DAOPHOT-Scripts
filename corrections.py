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
# Also standard aperture - applied to only ap mags
def zp_and_std_ape(row_of_df, start_dither):

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
	for alf in field_alf_files:

		# Load data
		df = pd.read_csv(alf, header=None, names=['ID', 'X', 'Y', 'M1', 'E1', 'M2', 'E2', 'M3', 'E3', 'M4', 'E4', 'M5', 'E5', 'chi', 'sharp'], skiprows=3, delim_whitespace=True)

		# Apply zmag to m (only to non 99 values)
		for index in range(0, len(df)):
			if df['M1'][index] != 99.9999:
				df.loc[index, 'M1'] = df['M1'][index] - 25 + zmag
			if df['M2'][index] != 99.9999:
				df.loc[index, 'M2'] = df['M2'][index] - 25 + zmag				
			if df['M3'][index] != 99.9999:
				df.loc[index, 'M3'] = df['M3'][index] - 25 + zmag
			if df['M4'][index] != 99.9999:
				df.loc[index, 'M4'] = df['M4'][index] - 25 + zmag
			if df['M5'][index] != 99.9999:
				df.loc[index, 'M5'] = df['M5'][index] - 25 + zmag

		# Write to new file
		filename = alf.replace('.alf', '.alf_zp')
		df.to_csv(filename, sep=' ')

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

				# Convert mag to flux
				flux = F0 * (10 ** (-ap_df['Mag'][z]/2.5))

				# Apply standard aperture calibration factor of 1.112
				flux = 1.112 * flux

				# Convert back to mag
				ap_df.ix[z, 'Mag'] = -2.5 * log10(flux/F0)

		# Write to new file
		filename = ap.replace('.ap', '.ap_zp')
		f = open(filename, 'w')
		f.write("ID  X  Y  Mag  Error \n")

		for z in range(0, len(ap_df)):
			f.writelines("%d %.3f %.3f %.4f %.4f \n" % (ap_df['ID'][z], ap_df['X'][z], ap_df['Y'][z], ap_df['Mag'][z], ap_df['Error'][z]))

		f.close()

	return(0)

# Aperture correction - uses ap and als mags for dither 1 to calculate correction value then applies to alf mags only
def ap_corr(row_of_df, start_dither):

	# Get names for dither 1 als and lst files - don't need ap file because lst has ap mags
	als = stem + '_d1_cbcd_dn.als'
	lst = stem + '_d1_cbcd_dn.lst'

	# Open files into df's
	df_lst = pd.read_csv(lst, skiprows=3, header=None, names=['ID_lst', 'X_lst', 'Y_lst', 'ap_mag', 'ap_error'], usecols=(0,1,2,3,4), delim_whitespace=True)
	df_als = pd.read_csv(als, skiprows=3, header=None, names=['ID_als', 'X_als', 'Y_als', 'psf_mag', 'psf_error'], usecols=(0,1,2,3,4), delim_whitespace=True)

	# Sort values by ID
	df_lst.sort_values('ID_lst', axis=0, inplace=True)
	df_als.sort_values('ID_als', axis=0, inplace=True)

	# Then cross-check the lst and als dataframes to only keep stars that appear in both
	for i in range(0, len(df_als)):
		match = 0
		for j in range(0, len(df_lst)):
			# Check whether IDs match
			if df_als['ID_als'][i] == df_lst['ID_lst'][j]:
				match += 1

		# If they match do nothing, else remove row
		if match == 0:
			df_als.drop(i, axis=0, inplace=True)

	# Reset index to start at 0
	df_als.reset_index(drop=True, inplace=True)
	df_lst.reset_index(drop=True, inplace=True)

	# Now get rid of any lst stars not in als
	for i in range(0, len(df_lst)):
		match = 0
		for j in range(0, len(df_als)):
			# Check whether IDs match
			if df_lst['ID_lst'][i] == df_als['ID_als'][j]:
				match += 1

		# If they match, do nothing, else remove row
		if match == 0:
			df_lst.drop(i, axis=0, inplace=True)

	# Reset index again
	df_als.reset_index(drop=True, inplace=True)
	df_lst.reset_index(drop=True, inplace=True)

	# Now df_als and df_lst have same stars so can concat 
	df = pd.concat((df_als, df_lst), axis=1)

	# Remove columns not needed
	df.drop(['X_lst', 'Y_lst', 'ID_lst'], axis=1, inplace=True) # keep als x and y as better

	# Calculate difference in mag
	df['Difference'] = df['ap_mag'] - df['psf_mag']

	# Calculate average difference
	corr_val = round(df['Difference'].mean(), 3)

	# Apply this correction value to all alf files
	for alf in field_alf_files:

		# Filename
		filename = alf.replace('.alf', '.alf_zp')

		# Open file into df
		data = pd.read_csv(filename, delim_whitespace=True, header=0)

		# Add corr_val to mag
		for index in range(0, len(data)):
			if data['M1'][index] != 99.9999:
				data.loc[index, 'M1'] = data['M1'][index] + corr_val
			if data['M2'][index] != 99.9999:
				data.loc[index, 'M2'] = data['M2'][index] + corr_val				
			if data['M3'][index] != 99.9999:
				data.loc[index, 'M3'] = data['M3'][index] + corr_val
			if data['M4'][index] != 99.9999:
				data.loc[index, 'M4'] = data['M4'][index] + corr_val
			if data['M5'][index] != 99.9999:
				data.loc[index, 'M5'] = data['M5'][index] + corr_val

		# Output to a new file suffiz .alf_apc
		filename = filename.replace('_zp', '_apc')
		data.to_csv(filename, sep=' ')

	return(0)

# Location correction - only need to apply to alf mags as won't be using ap mags after this
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

	# Loop over all 5 dithers and correct
	for alf in field_alf_files:

		# File that needs correcting is the alf_apc file
		filename = alf.replace('.alf','.alf_apc')

		# Load data into df
		df = pd.read_csv(filename, header=0, delim_whitespace=True)

		# Convert to flux, correct and convert back to mag for each star
		for index in range(0, len(df)):

			# Find the x and y coordinates corresponding to that star
			x_coord = int(np.floor(float(df['X'][index])))
			y_coord = int(np.floor(float(df['Y'][index])))

			# Find the corr value at that (x_coord, y_coord) pixel
			if (x_coord > 0 and x_coord < 256):
				x_coord = x_coord
			elif (x_coord <= 0):
				x_coord = 0
			else: x_coord = 255 # y_coord > 256

			if (y_coord > 0 and y_coord < 256):
				y_coord = y_coord
			elif (y_coord <= 0):
				y_coord = 0
			else: y_coord = 255 # y_coord > 256

			corr_val = corr_data[x_coord, y_coord]

			# For each mag column, convert to flux, apply correction and convert back to flux
			for index in range(0, len(df)):
				if df['M1'][index] != 99.9999:
					flux = 10 ** (df['M1'][index]/-2.5)
					flux = flux * corr_val
					df.loc[index, 'M1'] = -2.5 * log10(flux)
				if df['M2'][index] != 99.9999:
					flux = 10 ** (df['M2'][index]/-2.5)
					flux = flux * corr_val
					df.loc[index, 'M2'] = -2.5 * log10(flux)				
				if df['M3'][index] != 99.9999:
					flux = 10 ** (df['M3'][index]/-2.5)
					flux = flux * corr_val
					df.loc[index, 'M3'] = -2.5 * log10(flux)
				if df['M4'][index] != 99.9999:
					flux = 10 ** (df['M4'][index]/-2.5)
					flux = flux * corr_val
					df.loc[index, 'M4'] = -2.5 * log10(flux)
				if df['M5'][index] != 99.9999:
					flux = 10 ** (df['M5'][index]/-2.5)
					flux = flux * corr_val
					df.loc[index, 'M5'] = -2.5 * log10(flux)			

		# Make new filename by replacing _apc with _all
		filename = filename.replace('_apc', '_all')

		# Write csv to file
		df.to_csv(filename, sep=' ')

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

		print "Working on: " + str(target_name) + "    Epoch: " + str(epoch_number) + "     Field: " + str(field)

		# Set up names of ap, alf and field alf files - have 5 alf and 5 ap files for each epoch and 2 field alf files
		# Put them in three lists to make it easier to loop over later
		ap_files = []
		alf_files = []
		field_alf_files = []

		for dither in range(j, j+5):
			ap_files.append(stem + '_d' + str(dither) + '_cbcd_dn.ap')
			alf_files.append(stem + '_d' + str(dither) + '_cbcd_dn.alf')

		for f in [1,2]:
			field_alf_files.append(stem + '_f' + str(f) + '.alf')


		# Zero point
		zp_and_std_ape(i,j)

		# Aperture correction
		ap_corr(i,j)

		# Location correction
		loc_corr(i,j)
