'''
UPDATE: am testing whether after all the tests choosing the 6 brightest stars works any better
Purpose: Select candidate PSF stars from star list. 
Then do series of tests to remove stars that wouldn't make good psf stars
The following tests are carried out:
1. Is the star too close to the edge of the frame?
2. Is the star bright enough?
3. Does the star look star-like?
4. Is it just a hot pixel?
The brightest 6 stars that pass all the test are outputted to the lst file to be used in the psf creation process
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

# Import modules
import pandas as pd
import numpy as np
import sys
import os
import shutil
import astropy.io.fits as fits

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

		print "Dither " + str(i)

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

			# If X < 6 or X >250, drop row
			if row['X'] < 6 or row['X'] > 250:
				df.drop(index, inplace=True)
				print "Deleting star %d because it is too close to edge of frame" % index
				execute = 0 # don't need to carry out rest of tests

			if execute == 1:

				# If Y < 6 or Y > 250, drop row
				if row['Y'] < 6 or row['Y'] > 250:
					df.drop(index, inplace=True)
					print "Deleting star %d because it is too close to edge of frame" % index
					execute = 0 # don't need to carry out rest of tests

			# TEST 2 : NOT BRIGHT ENOUGH

			# Get x and y coords of the star in question
			x_coord = int(round(row['X'] - 1)) # zero-indexed in data and must be rounded to nearest integer
			y_coord = int(round(row['Y'] - 1)) # zero-indexed in data and must be rounded to nearest integer

			if execute == 1:
				if data[y_coord, x_coord] < 100:
					df.drop(index, inplace=True)
					print "Deleting star %d because it is not bright enough" % index
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
					print "Deleting star %d because it isn't star-like" % index
					execute = 0 # don't need to do rest of tests

			# TEST 4 : ELIMINATE HOT PIXELS

			if execute == 1:

				# Check ratio between centre and ring1 is > 10%
				if ave_1 / centre < 0.1:
					df.drop(index, inplace=True)
					print "Deleting star %d because it is just a hot pixel" % index
					execute = 0


			# TEST 5 : CLOSE NEIGHBOURS

			if execute == 1:

				# Initially set neighbour = 0 i.e. assume no neighbours
				neighbour = 0

				# Write out coordinates of star
				coords = []

				for i in range(-6,7):
					for j in range(-6,7):
						coords.append((y_coord+j, x_coord+i))

				# Check 20 x 20 grid for stars nearby
				for i in range(-10,11):
					for j in range(-10,11):

						value = 100

						# Check the value isn't in the star's area
						if (y_coord+j, x_coord+i) not in coords:
							if y_coord+j < 256 and y_coord >= 0:
								if x_coord+i < 256 and x_coord+i >= 0:

									# Check for neighbours
									if data[(y_coord+j, x_coord+i)] > value:
										detection = data[(y_coord+j, x_coord+i)]
										neighbour = 1

				if neighbour == 1:
					print "Deleting star %d because there is a pixel value of %f nearby" % (index, detection)
					df.drop(index, inplace=True)  

			# Write out final list of stars to the lst file in the correct format

			# Get header of lst file
			f = open(filename+'.lst', 'r')
			header = f.read().splitlines()[0:3]
			f.close()

			# Now overwrite this file
			f = open(filename+'.lst', 'w')
			f.writelines(header[0] + '\n' + header[1] + '\n' + header[2] + '\n')

			# Send stars to lst file
			df.to_csv(f, sep=' ', mode='a', header=None) #.iloc[0:6]

			f.close()






