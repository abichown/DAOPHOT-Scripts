'''
Purpose: Select candidate PSF stars from star list. 
Then do series of tests to remove stars that wouldn't make good psf stars
The following tests are carried out:
1. Is the star too close to the edge of the frame?
2. Is the star bright enough?
3. Does the star look star-like?
4. Is it just a hot pixel?
The stars that pass all the test are outputted to the lst file to be used in the psf creation process
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

# Import modules
import pandas as pd
import numpy as np
import sys
import os
import shutil

# find files 
stars = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star','Channel','Epoch'])

# Iterate over every line in text file
for i in range(0, len(stars)):

	# Grab galaxy, star name, channel and epoch
	galaxy = stars['Galaxy'][i]
	target_name = stars['Star'][i]

	if stars['Channel'][i] == 1:
		wavelength = '3p6um'
	elif stars['Channel'][i] == 2:
		wavelength = '4p5um'
	else: wavelength = 'channel not defined'

	if stars['Epoch'][i] < 10:
		epoch_number = '0' + str(stars['Epoch'][i])
	else: epoch_number = str(stars['Epoch'][i])

    # Find absolute path of image	
	home = '/home/ac833/Data/'
    
	cwd = home + str(galaxy) + '/BCD/' + target_name + '/ch' + str(stars['Channel'][i]) + '/e' + str(epoch_number) + '/'
	stem = target_name + '_' + wavelength + '_e' + epoch_number

	# Change directory to where image is
	os.chdir(cwd)

	print "Working on: " + stem

	# Loop over all 10 dithers
	for i in range(1,11):

		filename = stem + '_d' + str(i) + '_cbcd_dn' # get filename to use for this dither

		# Read in FITS image
		hdulist = fits.open(filename+'.fits')

		# Access the primary header-data unit (HDU)
		hdu = hdulist[0]
		data = hdu.data

		# Read in candidate PSF stars
		df = pd.read_csv(filename+'.lst', delim_whitespace=True, skiprows=3, header=None, names=['ID', 'X', 'Y', 'Mag', 'Error'], index_col=0)

		# Carry out all the tests on each star in the df
		for index, row in df.iterrows():

			execute = 1

			# TEST 1 : TOO CLOSE TO EDGE OF FRAME

			# If X < 15 or X >241, drop row
			if row['X'] < 15 or row['X'] > 241:
				df.drop(index, inplace=True)
				execute = 0 # don't need to carry out rest of tests

			if execute == 1:

				# If Y < 15 or Y > 241, drop row
				if row['Y'] < 15 or row['Y'] > 241:
					df.drop(index, inplace=True)
					execute = 0 # don't need to carry out rest of tests

			# TEST 2 : NOT BRIGHT ENOUGH

			# Get x and y coords of the star in question
			x_coord = int(round(row['X'] - 1)) # zero-indexed in data and must be rounded to nearest integer
			y_coord = int(round(row['Y'] - 1)) # zero-indexed in data and must be rounded to nearest integer

			if execute == 1:
				if data[y_coord, x_coord] < 150:
					df.drop(index, inplace=True)
					execute = 0 # don't need to carry out rest of tests

			# TEST 3 : NOT STAR-LIKE

			# Value of pixel at centre
			centre = data[y_coord, x_coord]

			# Set up of rings for the nearest neighbours
			ring1 = [data[(y_coord, x_coord - 1)], data[(y_coord, x_coord + 1)], data[(y_coord - 1, x_coord)], data[(y_coord + 1, x_coord)]]

			ring2 = [data[(y_coord, x_coord - 2)], data[(y_coord + 1, x_coord - 1)], data[(y_coord + 2, x_coord)], data[(y_coord + 1, x_coord + 1)], data[(y_coord, x_coord + 2)], data[(y_coord - 1, x_coord + 1)], data[(y_coord - 2, x_coord)], data[(y_coord - 1, x_coord - 1)]]

			ring3 = [data[(y_coord, x_coord - 3)], data[(y_coord + 1, x_coord - 2)], data[(y_coord + 2, x_coord - 1)], data[(y_coord + 3, x_coord)], data[(y_coord + 2, x_coord + 1)], data[(y_coord + 1, x_coord + 2)], data[(y_coord, x_coord + 3)], data[(y_coord - 1, x_coord + 2)], data[(y_coord - 2, x_coord + 1)], data[(y_coord - 3, x_coord)], data[(y_coord - 2, x_coord - 1)], data[(y_coord - 1, x_coord - 2)]]

			# Get average values of the rings
			ave_1 = sum(ring1)/len(ring1)
			ave_2 = sum(ring2)/len(ring2)
			ave_3 = sum(ring3)/len(ring3)

			if execute == 1:

				# Is centre > ave1 > ave2 > ave3?
				if ave_1 > centre or ave_2 > centre or ave_3 > centre or ave_2 > ave_1 or ave_3 > ave_1 or ave_3 > ave_2: 
					df.drop(index, inplace=True)
					execute = 0 # don't need to do rest of tests

			# TEST 4 : ELIMINATE HOT PIXELS

			if execute == 1:

				# Check ratio between centre and ring1 is > 10%
				if ave1 / centre < 0.1:
					df.drop(index, inplace=True)
					execute = 0

			# Write out final list of stars to the lst file in the correct format

			# Get header of lst file
			f = open(filename+'.lst', 'r')
			header = f.read().splitlines()[0:3]
			f.close()

			# Now overwrite this file
			f = open(filename+'.lst', 'w')
			f.writelines(header[0] + '\n' + header[1] + '\n' + header[2] + '\n')

			df.to_csv(f, sep=' ', mode='a', header=None)

			f.close()






