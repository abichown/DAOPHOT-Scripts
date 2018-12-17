'''
Purpose: FOR FIELD 1 I.E. OFF TARGET FIELDS ONLY!
Create a medianed image for each epoch - consists of 5 dithers per epoch.
Then obtain the master star list for that epoch which gets put through ALLFRAME in the next script.
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

	print galaxy, star_name, wavelength, num_epochs

	field = '1'
	start_dither = 1

	# Find absolute path of where images are
	home = '/home/ac833/Data/'

	# For each epoch, create median image of the 5 dithers
	for epoch in range(1, num_epochs+1):

		# Get epoch number into correct format
		if epoch < 10:
			epoch = '0' + str(epoch)
		else: epoch = str(epoch)

		# Change directory to where images are
		cwd = home + str(galaxy) + '/BCD/' + star_name + '/ch' + str(df['Channel'][i]) + '/e' + epoch + '/'
		os.chdir(cwd)

		# Copy options files to cwd
		shutil.copy('/home/ac833/daophot-options-files/daophot.opt', 'daophot.opt')
		shutil.copy('/home/ac833/daophot-options-files/photo.opt', 'photo.opt')
		shutil.copy('/home/ac833/daophot-options-files/allstar.opt', 'allstar.opt')

		# Run DAOMATCH
		daomatch = pexpect.spawn('daomatch')

		# Set up log file
		fout = file('daomatch_log.txt','w')
		daomatch.logfile = fout

		# Give DAOMATCH all the 5 dithers to be used in making the medianed image
		# The '/' after the file names lets DAOMATCH know that the scale is the same and it is only the rotation and shifts that might have changed
		daomatch.expect("Master input file:")
		daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch +'_d' + str(start_dither) + '_cbcd_dn.ap') # Give it the first BCD phot file
		daomatch.expect("Output file name")
		daomatch.sendline(star_name + '_' + wavelength +'_f'+field+'.mch')
		daomatch.expect("Next input file:")
		daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch +'_d' + str(start_dither+1) + '_cbcd_dn.ap/') # Give it the second BCD phot file
		daomatch.expect("Next input file")
		daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch +'_d' + str(start_dither+2) + '_cbcd_dn.ap/') # Give it the third BCD phot file
		daomatch.expect("Next input file")
		daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch +'_d' + str(start_dither+3) + '_cbcd_dn.ap/') # Give it the fourth BCD phot file
		daomatch.expect("Next input file")
		daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch +'_d' + str(start_dither+4) + '_cbcd_dn.ap/') # Give it the fifth BCD phot file
		daomatch.expect("Next input file")
		daomatch.sendline("") # exit

		# Run DAOMASTER
		daomaster = pexpect.spawn('daomaster')

		# Set up log file
		fout = file('daomaster_log.txt','w')
		daomaster.logfile = fout

		daomaster.expect("File with list of input files:")
		daomaster.sendline(star_name + '_' + wavelength +'_f'+field+'.mch')
		daomaster.expect("Minimum number, minimum fraction, enough frames:")
		daomaster.sendline("1, 0.5, 5") 
		daomaster.expect("Maximum sigma:")
		daomaster.sendline("99") 
		daomaster.expect("Your choice:")
		daomaster.sendline("6") # solve for 6 degrees of freedom
		daomaster.expect("Critical match-up radius:")
		daomaster.sendline("7") 
		
		for dither in range(start_dither+1,start_dither+5):
			daomaster.expect(star_name + '_' + wavelength + '_e' + epoch +'_d'+str(dither)+'_cbcd_dn.ap')
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
		daomaster.sendline("y")
		daomaster.expect("Output file name")
		daomaster.sendline(star_name + '_' + wavelength + '_f' + field + '.raw')
		daomaster.expect("A file with the new transformations?")
		daomaster.sendline("y")
		daomaster.expect("Output file name")
		daomaster.sendline(star_name + '_' + wavelength +'_f'+field+'.mch')
		daomaster.expect("New output file name")
		daomaster.sendline("")
		daomaster.expect("A file with the transfer table?")
		daomaster.sendline("e") # exits rest of options

		# Run MONTAGE2
		montage2 = pexpect.spawn('montage2')

		# Set up log file - need this for the offsets
		fout = file('montage_log.txt','w')
		montage2.logfile = fout

		montage2.expect("File with transformations:")
		montage2.sendline(star_name + '_' + wavelength +'_f'+field+'.mch')
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
		montage2.sendline(star_name + '_' + wavelength +'_f' + field + '.fits')

		# Write down X and Y offsets
		log = open('montage_log.txt', 'r')
		lines = log.readlines()

		offsets = []
		
		for line in lines:
			if "Offsets" in line:

				offsets.append(line.split(' ')[-3])
				offsets.append(line.split(' ')[-2])

		# Run DAOPHOT to get list of stars
		daophot = pexpect.spawn('daophot') 

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
		daophot.sendline("th=20") 
		daophot.expect("OPT>")
		daophot.sendline("")

		daophot.expect("Command:")
		daophot.sendline("fi")
		daophot.expect("Number of frames averaged, summed:")
		daophot.sendline("5,1") 
		daophot.expect("File for positions")
		daophot.sendline("")
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

		# Run candidate stars through tests to remove bad PSF stars

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
		# x_lo = centre[0] - (3*centre[0])/4
		# x_up = centre[0] + (3*centre[0])/4
		# y_lo = centre[1] - (3*centre[1])/4
		# y_up = centre[1] + (3*centre[1])/4
		x_lo = centre[0] - (9*centre[0])/10
		x_up = centre[0] + (9*centre[0])/10
		y_lo = centre[1] - (9*centre[1])/10
		y_up = centre[1] + (9*centre[1])/10

		df2 = pd.read_csv(star_name + '_' + wavelength + '_f' + field + '.lst', delim_whitespace=True, skiprows=3, header=None, names=['ID', 'X', 'Y', 'Mag', 'Error'], index_col=0)

		# Carry out all the tests on each star in the df
		for index, row in df2.iterrows():

			execute = 1

			# TEST 1 : TOO CLOSE TO EDGE OF FRAME

			# If X < x_lo or X > x_up, drop row
			if row['X'] < x_lo or row['X'] > x_up:
				df2.drop(index, inplace=True)
				print "Deleting star %d because it is too close to edge of frame" % index
				execute = 0 # don't need to carry out rest of tests

			if execute == 1:

				# If Y < y_lo or Y > y_up, drop row
				if row['Y'] < y_lo or row['Y'] > y_up:
					df2.drop(index, inplace=True)
					print "Deleting star %d because it is too close to edge of frame" % index
					execute = 0 # don't need to carry out rest of tests

			# TEST 2 : NOT BRIGHT ENOUGH

			# Get x and y coords of the star in question
			x_coord = int(round(row['X'] - 1)) # zero-indexed in data and must be rounded to nearest integer
			y_coord = int(round(row['Y'] - 1)) # zero-indexed in data and must be rounded to nearest integer

			if execute == 1:
				if data[y_coord, x_coord] < 150:
					df2.drop(index, inplace=True)
					print "Deleting star %d because it is not bright enough" % index
					execute = 0 # don't need to carry out rest of tests

			# Write out final list of stars to the lst file in the correct format

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

		# Make PSF model
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

		# Run ALLSTAR
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
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.mag') # this is the final file for star positions that we want 

		daophot.expect("Command:")
		daophot.sendline("ex")
		daophot.close(force=True)

		# Change .ap file names in .mch file to .als ready for next script
		f = open(star_name +'_' + wavelength + '_f'+field+'.mch', 'r')
		filedata = f.read()
		f.close()

		newdata = filedata.replace(".ap",".als")

		f = open(star_name +'_' + wavelength + '_f'+field+'.mch','w')
		f.write(newdata)
		f.close()	

		# Rename files we want for the next scripts to be consistent with the on-target version of this code
		os.rename(star_name + '_' + wavelength + '_f' + field + '.fits', star_name + '_' + wavelength + '_f' + field + '_master.fits')
		os.rename(star_name + '_' + wavelength + '_f' + field + '.mch', star_name + '_' + wavelength + '_f' + field + '_master.mch')
		os.rename(star_name + '_' + wavelength + '_f' + field + '.mag', star_name + '_' + wavelength + '_f' + field + '_master.mag')
		os.rename(star_name + '_' + wavelength + '_f' + field + '.psf', star_name + '_' + wavelength + '_f' + field + '_master.psf')