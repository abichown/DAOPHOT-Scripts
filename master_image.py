'''
Purpose: Create a medianed image from all the 120 (24 x 5) or 60 (12 x 5) frames for a field for a channel.
Then obtain the master star list which gets put through ALLFRAME in the next script.
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

# Import modules
import pandas as pd
import os
import pexpect
import shutil
import numpy as np
import astropy.io.fits as fits
import sys

import matplotlib.pyplot as plt
from kneed import DataGenerator, KneeLocator

# Get list of fields to make master images for 
df = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star', 'Period', 'RA', 'Dec', 'Channel'])

# No longer need these as won't have duplicates
#df.drop_duplicates(inplace=True)
#df.reset_index(drop=True, inplace=True)

for i in range(0,len(df)):

	galaxy = str(df['Galaxy'][i]) 
	star_name = str(df['Star'][i]) 
	channel = str(df['Channel'][i]) 

	if channel == '1':
		wavelength = '3p6um'
	else: wavelength = '4p5um'

	if galaxy == 'LMC':
		num_epochs = 24
	elif galaxy == 'SMC':
		num_epochs = 12
	else: num_epochs = 0 # will error later on in script

	num_images = num_epochs * 5 # for a 5 dither pattern like we have for CHP data

	print galaxy, star_name, wavelength, num_epochs

	# Create list of all files to be used to make the image
	field1_files = []
	field2_files = []

	# Populate lists

	# Field 1
	for epoch in range(1,num_epochs+1):
		for dither in range(1,6):

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			dither = str(dither)

			filename = star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.ap'
			field1_files.append(filename)


	# Field 2
	for epoch in range(1,num_epochs+1):
		for dither in range(6,11):

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			dither = str(dither)

			filename = star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.ap'
			field2_files.append(filename)


	temp = '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/temp/'

	# Delete temp folder if it already exists - NEED TO ADD IF IT EXISTS TO DELETE
	if (os.path.isdir(temp)):
		shutil.rmtree(temp)

	# Make temp folder
	os.mkdir(temp)

	for start_dither in [6]: # [1,6] just checking field 2 works right now 

		if start_dither == 1:
			field = '1'
		else: field = '2'

		# Copy all FITS images, psfs and phot files to temp folder
		for epoch in range(1,num_epochs+1):

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			cwd = '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'

			for dither in range(start_dither,start_dither+5):

				dither = str(dither)

				shutil.copyfile(cwd + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.ap', temp + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.ap') #.ap
				shutil.copyfile(cwd + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.fits', temp + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.fits')
				#shutil.copyfile(cwd + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.psf', temp + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.psf')

		# Change to temp folder 
		os.chdir(temp)

		# Use DAOMATCH to get initial transformations
		daomatch = pexpect.spawn('daomatch', timeout=30)

		fout = file('log.txt','w')
		daomatch.logfile = fout

		# Input epoch 1 first dither of field as master input file
		daomatch.expect("Master input file:")

		if field == '1':
			daomatch.sendline(field1_files[0]) # this is epoch 1 dither 1 file
			daomatch.expect("Output file name")
			daomatch.sendline(star_name + '_' + wavelength + '_f1.mch')
		else: 
			daomatch.sendline(field2_files[0]) # this is epoch 1 dither 1 file
			daomatch.expect("Output file name")
			daomatch.sendline(star_name + '_' + wavelength + '_f2.mch')		


		for j in range(1,num_images):

			daomatch.expect("Next input file")

			if field == '1':
				daomatch.sendline(field1_files[j]+'/')
			else:
				daomatch.sendline(field2_files[j]+'/')


		daomatch.expect("Next input file")
		daomatch.sendline("") # exit


		# Use DAOMASTER to refine the transformations
		daomaster = pexpect.spawn('daomaster')

		fout = file('daomaster_log.txt','w')
		daomaster.logfile = fout

		print "Running DAOMASTER"

		daomaster.expect("File with list of input files:")
		daomaster.sendline(star_name + '_' + wavelength + '_f' + field + '.mch')
		daomaster.expect("Minimum number, minimum fraction, enough frames:")
		daomaster.sendline("1, 0.5, " + str(num_images)) # play around with these values
		daomaster.expect("Maximum sigma:")
		daomaster.sendline("99") # play around with this value
		daomaster.expect("Your choice:")
		daomaster.sendline("6") # solve for 6 degrees of freedom
		daomaster.expect("Critical match-up radius:")
		daomaster.sendline("7") 

		for j in range(1,num_images):

			if field == '1':
				daomaster.expect(field1_files[j])
				daomaster.sendline("")
			else:
				daomaster.expect(field2_files[j])
				daomaster.sendline("")

		# Repeat with decreasing match up size
		for match_up in range(7,-1,-1):
			daomaster.expect("New match-up radius")
			daomaster.sendline(str(match_up))


		# Options for different output files - only want transformations according to cookbook
		daomaster.expect("Assign new star IDs?")
		daomaster.sendline("y") # assign new ids so all frames have same ids
		daomaster.expect("A file with mean magnitudes and scatter?")
		daomaster.sendline("n")
		daomaster.expect("A file with corrected magnitudes and errors?")
		daomaster.sendline("n")
		daomaster.expect("A file with raw magnitudes and errors?")
		daomaster.sendline("n")
		daomaster.expect("A file with the new transformations?")
		daomaster.sendline("y")
		daomaster.expect("Output file name")
		daomaster.sendline(star_name + '_' + wavelength + '_f' + field + '.mch')
		daomaster.expect("New output file name")
		daomaster.sendline("")
		daomaster.expect("A file with the transfer table?")
		daomaster.sendline("e") # exits rest of options

		print "Running MONTAGE2"

		# Use MONTAGE2 to actually make image
		montage2 = pexpect.spawn('montage2')

		# Set up log file
		fout = file('montage_log.txt','w')
		montage2.logfile = fout

		montage2.expect("File with transformations:")
		montage2.sendline(star_name + '_' + wavelength + '_f' + field + '.mch')
		montage2.expect("Image-name suffix:")
		montage2.sendline("")
		montage2.expect("Minimum number of frames, percentile:")
		montage2.sendline("1,0.5") # play around with minimum number of frames
		montage2.expect("X limits of output image:")
		montage2.sendline("e")
		montage2.expect("Y limits of output image:")
		montage2.sendline("e")
		montage2.expect("Expansion factor:")
		montage2.sendline("1") # creates image with same scale as bcd images
		montage2.expect("Determine sky from overlap region?")
		montage2.sendline("y")
		montage2.expect("Name for output image")
		montage2.sendline(star_name + '_' + wavelength + '_f' + field + '.fits')
		montage2.expect("Good bye")
		montage2.close(force=True)

		mch_file = star_name + '_' + wavelength + '_f' + field +'.mch'
		shutil.copy(mch_file, star_name + '_' + wavelength + '_f' + field +'_full.mch') # want to keep list of all frames

		# Check whether weightings are good
		check = False

		while check == False:

			f = open('montage_log.txt', 'r')

			weights = []

			for line in f:
				y = line.split()

				for i in range(1,len(y)):
					if y[i] == 'weight':
						weights.append(y[i+2]) # this is the weights value

			weights = np.array(weights) # convert to array
			weights = weights.astype(float) # convert from strings to floats

			bad_frame_list = []

			# Check for bad frames
			for i in weights:
				if i < 200:
					bad_frame_list.append(i)

			if len(bad_frame_list) > 0:
				bad = True 
			else: 
				bad = False
				check = True # no bad frames so don't need to do anything

			# If there was a bad weighting, then want to remove this image from the match file
			if bad == True:

				bad_frame_value = max(weights) # value of the bad frame

				g = open('montage_log.txt', 'r')

				for line in g:
					y = line.split()

					for i in range(1,len(y)):
						if y[i] == 'weight':
							if float(y[i+2]) == bad_frame_value:
								bad_frame = y[i+4] # this gets the file name of the bad frame

				# Remove bad frame 

				with open(mch_file) as oldfile, open(mch_file+'_new', 'w') as newfile:
					for line in oldfile:
						if bad_frame not in line:
							newfile.write(line)

				shutil.copy(mch_file+'_new',  mch_file)

				# Re run MONTAGE2
				montage2 = pexpect.spawn('montage2')

				# Set up log file
				fout = file('montage_log.txt','w')
				montage2.logfile = fout

				montage2.expect("File with transformations:")
				montage2.sendline(mch_file)
				montage2.expect("Image-name suffix:")
				montage2.sendline("")
				montage2.expect("Minimum number of frames, percentile:")
				montage2.sendline("1,0.5") # play around with minimum number of frames
				montage2.expect("X limits of output image:")
				montage2.sendline("e")
				montage2.expect("Y limits of output image:")
				montage2.sendline("e")
				montage2.expect("Expansion factor:")
				montage2.sendline("1") # creates image with same scale as bcd images
				montage2.expect("Determine sky from overlap region?")
				montage2.sendline("y")
				montage2.expect("Name for output image")
				montage2.sendline(star_name + '_' + wavelength + '_f' + field + '.fits')
				montage2.expect("Good bye")
				montage2.close(force=True)				


		# Write down X and Y offsets
		log = open('montage_log.txt', 'r')
		lines = log.readlines()

		offsets = []
		
		for line in lines:
			if "Offsets" in line:

				offsets.append(line.split(' ')[-3])
				offsets.append(line.split(' ')[-2])


		# Use DAOPHOT to create star list
		shutil.copy('/home/ac833/daophot-options-files/daophot.opt', 'daophot.opt')
		shutil.copy('/home/ac833/daophot-options-files/photo.opt', 'photo.opt')
		shutil.copy('/home/ac833/daophot-options-files/allstar.opt', 'allstar.opt')

		daophot = pexpect.spawn('daophot') 

		print "Running DAOPHOT"

		# Set up logfile
		fout = file('daophot_log.txt','w')
		daophot.logfile = fout

		# Attach medianed image
		daophot.expect("Command:")
		daophot.sendline("at " + star_name + '_' + wavelength + '_f' + field +'.fits')
		daophot.expect("Command:")
		daophot.sendline("opt")
		daophot.expect("File with parameters")
		daophot.sendline("")
		daophot.expect("OPT>")
		#daophot.sendline("th=2")
		daophot.sendline("th=20") 
		daophot.expect("OPT>")
		daophot.sendline("")

		# daophot.expect("Command:")
		# daophot.sendline("fi")
		# daophot.expect("Number of frames averaged, summed:")
		# daophot.sendline(str(num_images)+",1")
		# daophot.expect("File for positions")
		# daophot.sendline("")

		# # Needs to store the number of detections for the threshold
		# for threshold in range(3,61):

		# 	daophot.expect("Are you happy with this?")
		# 	daophot.sendline("n")
		# 	daophot.expect("New threshold")
		# 	daophot.sendline(str(threshold))
		# 	daophot.expect("Output file name")
		# 	daophot.sendline("")
		# 	daophot.expect("New output file name")
		# 	daophot.sendline("")

		# # Last threshold = 60
		# daophot.expect("Are you happy with this?")
		# daophot.sendline("y")

		# # Open log file, extract number of detections at each threshold
		# log = open('daophot_log.txt', 'r')
		# split = log.read().split()

		# num_det = []

		# for i in range(0,58):
		# 	index = -8 -i*40
		# 	num_det.append(split[index])

		# # Reverse order so threshold increases and detections decrease
		# num_det = num_det[::-1]
		# thresholds = range(2,60)

		# print thresholds
		# print num_det

		# plt.clf()

		# plt.scatter(thresholds, num_det)
		# plt.show()

		# # Package that finds the elbow/knee of a curve - this will be the optimal detection threshold
		# kn = KneeLocator(thresholds, num_det, direction='decreasing')
		# print "Knee = " + str(kn.knee)


		# daophot.expect("Are you happy with this?")
		# daophot.sendline("y")

		# # Then run FIND with this desired threshold
		# daophot.expect("Command:")
		# daophot.sendline("opt")
		# daophot.expect("File with parameters")
		# daophot.sendline("")
		# daophot.expect("OPT>")
		# daophot.sendline("th="+str(kn.knee)) # now set to best threshold
		# daophot.expect("OPT>")
		# daophot.sendline("")

		daophot.expect("Command:")
		daophot.sendline("fi")
		daophot.expect("Number of frames averaged, summed:")
		daophot.sendline(str(num_images)+",1") 
		daophot.expect("File for positions")
		daophot.sendline("")
		# daophot.expect("New output file name")
		# daophot.sendline("")
		daophot.expect("Are you happy with this?")
		daophot.sendline("y")

		daophot.expect("Command:")
		daophot.sendline("ph")
		daophot.expect("File with aperture radii")
		daophot.sendline("")
		daophot.expect("PHO>")
		daophot.sendline("")
		daophot.expect("Input position file")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.coo')
		daophot.expect("Output file")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.ap')

		# Now choose brightest 20 stars in image to be candidate PSF stars
		daophot.expect("Command:")
		daophot.sendline("pi")
		daophot.expect("Input file name")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.ap')
		daophot.expect("Desired number of stars, faintest magnitude:")
		daophot.sendline("20,99")
		daophot.expect("Output file name")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.lst') 
		daophot.expect("Command:")
		daophot.sendline("ex")

		daophot.close(force=True)

		# Now run these candidate PSF stars through the series of tests to get rid of bad stars

		# Read in FITS image
		hdulist = fits.open(star_name + '_' + wavelength + '_f' + field +'.fits')

		# Access the primary header-data unit (HDU)
		hdu = hdulist[0]
		data = hdu.data

		# Obtain the length of the x and y axis of the image
		x_axis = hdulist[0].header['NAXIS1']
		y_axis = hdulist[0].header['NAXIS2']

		centre = [x_axis/2, y_axis/2] # centre of frame

		# Obtain lower and upper x and y limits for Test 1
		x_lo = centre[0] - (3*centre[0])/4
		x_up = centre[0] + (3*centre[0])/4
		y_lo = centre[1] - (3*centre[1])/4
		y_up = centre[1] + (3*centre[1])/4
		# x_lo = centre[0] - (9*centre[0])/10
		# x_up = centre[0] + (9*centre[0])/10
		# y_lo = centre[1] - (9*centre[1])/10
		# y_up = centre[1] + (9*centre[1])/10

		df2 = pd.read_csv(star_name + '_' + wavelength + '_f' + field + '.lst', delim_whitespace=True, skiprows=3, header=None, names=['ID', 'X', 'Y', 'Mag', 'Error'], index_col=0)

		deleted_stars = 0 

		# Carry out all the tests on each star in the df
		for index, row in df2.iterrows():

			execute = 1

			# TEST 1 : TOO CLOSE TO EDGE OF FRAME

			# If X < x_lo or X > x_up, drop row
			if row['X'] < x_lo or row['X'] > x_up:
				df2.drop(index, inplace=True)
				deleted_stars += 1
				#print "Deleting star %d because it is too close to edge of frame" % index
				execute = 0 # don't need to carry out rest of tests

			if execute == 1:

				# If Y < y_lo or Y > y_up, drop row
				if row['Y'] < y_lo or row['Y'] > y_up:
					df2.drop(index, inplace=True)
					deleted_stars += 1
					#print "Deleting star %d because it is too close to edge of frame" % index
					execute = 0 # don't need to carry out rest of tests

			# TEST 2 : NOT BRIGHT ENOUGH

			# Get x and y coords of the star in question
			x_coord = int(round(row['X'] - 1)) # zero-indexed in data and must be rounded to nearest integer
			y_coord = int(round(row['Y'] - 1)) # zero-indexed in data and must be rounded to nearest integer

			if execute == 1:
				if data[y_coord, x_coord] < 150:
					df2.drop(index, inplace=True)
					deleted_stars += 1
					#print "Deleting star %d because it is not bright enough" % index
					execute = 0 # don't need to carry out rest of tests

			# Write out final list of stars to the lst file in the correct format

			print "Stars remaining: " + str(20-deleted_stars)

			# Get header of lst file
			f = open(star_name + '_' + wavelength + '_f' + field + '.lst', 'r')
			header = f.read().splitlines()[0:3]
			f.close()

			# Now overwrite this file
			f = open(star_name + '_' + wavelength + '_f' + field + '.lst', 'w')
			f.writelines(header[0] + '\n' + header[1] + '\n' + header[2] + '\n')

			# Send stars to lst file
			df2.to_csv(f, sep=' ', mode='a', header=None) #.iloc[0:6]

			f.close()

		daophot = pexpect.spawn('daophot')

		daophot.expect("Command:")
		daophot.sendline("at " + star_name + '_' + wavelength + '_f' + field +'.fits')

		# Now create PSF model to be used in allstar in next step of this script
		daophot.expect("Command:")
		daophot.sendline("psf")
		daophot.expect("File with aperture results")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.ap')
		daophot.expect("File with PSF stars")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.lst')
		daophot.expect("File for the PSF")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.psf')

		daophot.expect("Command:")
		daophot.sendline("ex")

		daophot.close(force=True)

		# Open ALLSTAR
		allstar = pexpect.spawn('allstar')

		fout = file('allstar_log.txt', 'w')
		allstar.logfile = fout

		allstar.expect("OPT>")
		allstar.sendline("")
		allstar.expect("Input image name:")
		allstar.sendline(star_name + '_' + wavelength + '_f'+ field + '.fits')
		allstar.expect("File with the PSF")
		allstar.sendline(star_name + '_' + wavelength + '_f' + field + '.psf')

		allstar.expect("Input file")
		allstar.sendline(star_name + '_' + wavelength + '_f' + field + '.ap')
		allstar.expect("File for results")
		allstar.sendline(star_name + '_' + wavelength + '_f' + field + '.als')
		allstar.expect("Name for subtracted image")
		allstar.sendline(star_name + '_' + wavelength + '_f' + field + '_dns.fits')

		allstar.expect("Good bye")
		allstar.close(force=True)

		# Run DAOPHOT to add offsets back in
		daophot = pexpect.spawn('daophot')

		daophot.expect("Command:")
		daophot.sendline("off") # offsets to put x and y back in
		daophot.expect("Input file name:")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.als')
		daophot.expect("Additive offsets ID, DX, DY, DMAG:")
		daophot.sendline("0," + offsets[0] + "," + offsets[1] + ",0")
		daophot.expect("Output file name")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.mag')

		daophot.expect("Command:")
		daophot.sendline("ex")
		daophot.close(force=True)	

		# Change .ap file names in .mch file to .als
		f = open(star_name +'_' + wavelength + '_f'+field+'_full.mch', 'r')
		filedata = f.read()
		f.close()

		newdata = filedata.replace(".ap",".als")

		f = open(star_name +'_' + wavelength + '_f'+field+'_full.mch','w')
		f.write(newdata)
		f.close()	

		# Now copy the master image and the master star list to each of the epochs 
		for epoch in range(1,num_epochs+1):

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			shutil.copy(star_name + '_' + wavelength + '_f' + field + '.fits', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'+star_name + '_' + wavelength + '_f' + field + '_master.fits')
			shutil.copy(star_name + '_' + wavelength + '_f' + field + '_full.mch', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'+star_name + '_' + wavelength + '_f' + field + '_master.mch')
			shutil.copy(star_name + '_' + wavelength + '_f' + field + '.mag', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'+star_name + '_' + wavelength + '_f' + field + '_master.mag')
			shutil.copy(star_name + '_' + wavelength + '_f' + field + '.psf', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'+star_name + '_' + wavelength + '_f' + field + '_master.psf')

	# Delete temp folder
	shutil.rmtree(temp)