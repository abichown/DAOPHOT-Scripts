'''
Purpose: Extract magnitude and error of the V*. At the moment, V* always in centre of frame
so program looks for star near (133,122). Also, V* only present in [3.6] field 2 and 
[4.5] field 1 so only looks there.
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

# Import modules
import os
import sys
import pandas as pd
import numpy as np
from astropy import wcs
from astropy.io import fits
from math import sqrt


# Read star list
df = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star', 'Period', 'RA', 'Dec', 'Channel'])

for i in range(0, len(df)):

	# Then get galxy name
	galaxy = df['Galaxy'][i]
	#galaxy = df['Galaxy'][index]

	# Get star name
	star = df['Star'][i]

    # Get channel and convert to wavelength
	if df['Channel'][i] == 1:
		channel = '1'
		wavelength = '3p6um'
		field = '2' # we want the field with the V* in it
		start_dither = '6'
	elif df['Channel'][i] == 2:
		channel = '2'
		wavelength = '4p5um'
		field = '1' # we want the field with the V* in it
		start_dither = '1'
	else: wavelength = 'channel not defined'

	# Number of epochs
	if galaxy == 'LMC':
		num_epochs = 24
	elif galaxy == 'SMC':
		num_epochs = 12
	else: num_epochs = 0

	# Get RA and Dec of Cepheid
	ra = df['RA'][i]
	dec = df['Dec'][i]

	# Open file to write mag and error to 
	filename = '/home/ac833/Magnitudes/' + galaxy + '/' + star + '_' + wavelength +'.txt'
	f = open(filename, 'w')
	f.write("Epoch Mag Error Std_err \n")

	count = 0

	print star

	for epoch in range(1,num_epochs+1):

		# Get epoch in correct form
		if epoch < 10:
			epoch = '0' + str(epoch)
		else: epoch = str(epoch)

		# Change cwd to folder with data in
		cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star + '/ch' + channel + '/e' + epoch + '/'
		os.chdir(cwd)

		# Convert RA and Dec to pixel coordinates using astropy.wcs
		hdulist = fits.open(star + '_' + wavelength + '_e' + epoch + '_d' + start_dither + '_cbcd_dn.fits')		
		w = wcs.WCS(hdulist[0].header) # Parse the WCS keywords in the primary HDU
		world = np.array([[ra, dec]], np.float_) # the world coordinates of the target Cepheid
		pix = w.wcs_world2pix(world,1) # convert world coordinates to pixel coordinates in an array. For just list of coords want pix[0]


		# Coordinates we want to compare .ave file with
		pix_coord = pix[0] # just as a list. x = pix_coord[0] and y = pix_coord[1]

		print "Pixel coordinates from NED for epoch " + str(epoch) + " : x = " + str(pix_coord[0]) + "y = " + str(pix_coord[1])

		# Open file containing average magnitudes
		ave_file = star + '_' + wavelength + '_f' + field + '.ave'

		# # THIS METHOD OF GETTING THE STAR USES ASTROPY.WCS BUT CURRENTLY ISN'T WORKING ON HV00914
		# if (os.path.isfile(ave_file)):

		# 	ave = pd.read_csv(ave_file, header=0, delim_whitespace=True)

		# 	potential_id = []
		# 	potential_x = []
		# 	potential_y = []
		# 	potential_mags = []

		# 	for i in range(0, len(ave)):

		# 		# Find stars that could potentially be the target star because they are in the right area of the frame 
		# 		if abs(ave['X'][i]-pix_coord[0]) <= 8 and abs(ave['Y'][i]-pix_coord[1]) <= 8:

		# 			potential_id.append(ave['ID'][i])
		# 			potential_x.append(ave['X'][i])
		# 			potential_y.append(ave['Y'][i])
		# 			potential_mags.append(ave['m_ave'][i])

		# 	# Now want to work out which star is closest to the world coordinates

		# 	# Define function to determine minimum 
		# 	def sq_diff(x, y):
		# 		return sqrt(x**2 + y**2)

		# 	potential_sq_diff = []

		# 	# Compute this difference value for all the stars that could potentially be the Cepheid
		# 	for i in range(0,len(potential_id)):
		# 		potential_sq_diff.append(sq_diff(potential_x[i]-pix_coord[0], potential_y[i]-pix_coord[1]))

		# 	sq_val = 100
		# 	true_index = 99

		# 	# Obtain the index which has the smallest difference value - this is the star we want
		# 	for i in range(0,len(potential_id)):
		# 		if potential_sq_diff[i] < sq_val:
		# 			sq_val = potential_sq_diff[i]
		# 			true_index = i

		# 	cepheid_mag = potential_mags[true_index]
				
		# 	#cepheid_mag = min(potential_mags) # this is the brightest of the stars in that region and is the Cepheid

		# 	# Find row which has this magnitude in it
		# 	stop = False
		# 	i = 0

		# 	while stop == False:

		# 		if cepheid_mag == ave['m_ave'][i]:
					
		# 			# Write this to row of the dataframe to file
		# 			f.writelines("%s %f %f %f \n" % (epoch, float(ave['m_ave'][i]), float(ave['e_ave'][i]), float(ave['std_err'][i])))
		# 			count += 1
		# 			stop = True

		# 		else: i += 1

		# THIS METHOD OF GETTING THE STAR USES ASTROPY.WCS BUT CURRENTLY ISN'T WORKING ON HV00914
		if (os.path.isfile(ave_file)):

			ave = pd.read_csv(ave_file, header=0, delim_whitespace=True)

			# potential_id = []
			# potential_x = []
			# potential_y = []
			# potential_mags = []

			for i in range(0, len(ave)):

				# Find stars that could potentially be the target star because they are in the right area of the frame 
				if abs(ave['X'][i]-pix_coord[0]) <= 3 and abs(ave['Y'][i]-pix_coord[1]) <= 3:

					f.writelines("%s %f %f %f %f %f %f \n" % (epoch, float(ave['ID'][i]), float(ave['X'][i]), float(ave['Y'][i]), float(ave['m_ave'][i]), float(ave['e_ave'][i]), float(ave['std_err'][i])))
					count += 1

		# # THIS METHOD OF GETTING STAR ASSUMES IT IS NEAR THE CENTRE OF THE IMAGE
		# if (os.path.isfile(ave_file)):
		# 	ave = pd.read_csv(ave_file, header=0, delim_whitespace=True)

		# 	# Search for star near (133,122)
		# 	# I've played with these ranges to ensure HV00872 ch1 and ch2 data get included
		# 	# Might need to be changed
		# 	for i in range(0, len(ave)):
		# 		if (ave['X'][i] < 135.00 and ave['X'][i] > 128.00):
		# 			if (ave['Y'][i] < 125.00 and ave['Y'][i] > 119.00):
		# 				# write to file
		# 				f.writelines("%s %f %f %f \n" % (epoch, float(ave['m_ave'][i]), float(ave['e_ave'][i]), float(ave['std_err'][i])))
		# 				count += 1

	print count

	f.close()
