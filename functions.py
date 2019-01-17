'''
Purpose: Functions to be used in master_script.py
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

# Import modules
import sys
import pexpect
import shutil
import os
import fnmatch
import re
import pandas as pd
import time
import numpy as np
import astropy.io.fits as fits
from math import log10, sqrt, exp
from astropy import wcs

import matplotlib.pyplot as plt
from kneed import DataGenerator, KneeLocator

import matplotlib
import gloess_fits as gf
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
matplotlib.rc('text',usetex=True)
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Garamond']

# Aperture photometry
def aper_phot(star_name, galaxy, channel, wavelength, epoch_number):

	# Iterate over each of the 10 dithers for this epoch
	for j in range(1,11):

		dither = str(j) # current dither number

		# Filename to work on 
		file_stem = star_name + '_' + wavelength + '_e' + epoch_number + '_d' + dither + '_cbcd_dn'

		# Copy DAOPHOT options file to current working directory
		shutil.copy('/home/ac833/daophot-options-files/daophot.opt', 'daophot.opt')
		shutil.copy('/home/ac833/daophot-options-files/photo.opt', 'photo.opt')

		# Run DAOPHOT to perform the aperture photometry
		daophot = pexpect.spawn("daophot")

		fout = file('log.txt', 'w')
		daophot.logfile = fout		

		# Attach the image
		daophot.expect("Command:")
		daophot.sendline("at " + file_stem + '.fits')

		# Find the stars
		daophot.expect("Command:")
		daophot.sendline("fi")
		daophot.expect("Number of frames averaged, summed:")
		daophot.sendline("1,1")
		daophot.expect("File for positions")
		daophot.sendline("")
		daophot.expect("Are you happy with this?")
		daophot.sendline("y")

		# Perform aperture photometry
		daophot.expect("Command:")
		daophot.sendline("ph")
		daophot.expect("File with aperture radii")
		daophot.sendline("")
		daophot.expect("PHO>")
		daophot.sendline("")
		daophot.expect("Input position file")
		daophot.sendline(file_stem + ".coo")
		daophot.expect("Output file")
		daophot.sendline(file_stem + ".ap")

		# Choose candidate PSF stars
		daophot.expect("Command:")
		daophot.sendline("pi")
		daophot.expect("Input file name")
		daophot.sendline(file_stem + '.ap')
		daophot.expect("Desired number of stars, faintest magnitude:")
		daophot.sendline("15,99") # used for the aperture correction later
		daophot.expect("Output file name")
		daophot.sendline(file_stem + '.lst')   

		# Exit daophot
		daophot.expect("Command:")
		daophot.sendline("exit")
		daophot.close(force=True)

	return(0)

# Master image on target
def master_on_target(star_name, galaxy, channel, wavelength, epoch_number, num_epochs):

	num_images = num_epochs * 5 # 5 dithers per epoch 
	field = '2' # this image is referred to as field 2 while the off-target field is field 1
	start_dither = 6 # field 2 dithers run from d6 to d10

	# List of files to use
	files = []

	for epoch in range(1,num_epochs+1):
		for dither in range(6,11): # need dithers 6 - 10 as this is the on-target field

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			dither = str(dither)

			filename = star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.ap'
			files.append(filename)	

	# Create temporary folder to make image - will then be copied to each epoch folder
	temp = '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/temp/'

	# Delete temp folder if it already exists
	if (os.path.isdir(temp)):
		shutil.rmtree(temp)

	# Make temp folder
	os.mkdir(temp)

	# Copy all FITS images and aperture photometry files to temp folder
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

	# Change to temp folder 
	os.chdir(temp)

	################################################################################################
	# 							MAKE FILE OF COORDINATE TRANSFORMATIONS
	################################################################################################

	# Use DAOMATCH to get initial transformations
	daomatch = pexpect.spawn('daomatch')

	fout = file('log.txt','w')
	daomatch.logfile = fout

	# Input epoch 1 first dither of field as master input file
	daomatch.expect("Master input file:")
	daomatch.sendline(files[0]) # this is epoch 1 dither 1 file
	daomatch.expect("Output file name")
	daomatch.sendline(star_name + '_' + wavelength + '_f2.mch')	

	# Give it the rest of the images
	for j in range(1,num_images):

		daomatch.expect("Next input file")
		daomatch.sendline(files[j]+'/')

	daomatch.expect("Next input file")
	daomatch.sendline("") # exit

	# Use DAOMASTER to refine the transformations
	daomaster = pexpect.spawn('daomaster')

	fout = file('daomaster_log.txt','w')
	daomaster.logfile = fout

	daomaster.expect("File with list of input files:")
	daomaster.sendline(star_name + '_' + wavelength + '_f2.mch')
	daomaster.expect("Minimum number, minimum fraction, enough frames:")
	daomaster.sendline("1, 0.5, " + str(num_images)) # play around with these values
	daomaster.expect("Maximum sigma:")
	daomaster.sendline("99") 
	daomaster.expect("Your choice:")
	daomaster.sendline("6") # solve for 6 degrees of freedom
	daomaster.expect("Critical match-up radius:")
	daomaster.sendline("7") 

	for j in range(1,num_images):

		daomaster.expect(files[j])
		daomaster.sendline("")

	# Reduce the match up radius 
	for match_up in range(7,-1,-1):
		daomaster.expect("New match-up radius")
		daomaster.sendline(str(match_up))	

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
	daomaster.sendline(star_name + '_' + wavelength + '_f2.mch') # these are more refined transformations
	daomaster.expect("New output file name")
	daomaster.sendline("")
	daomaster.expect("A file with the transfer table?")
	daomaster.sendline("e") # exits rest of options

	################################################################################################
	# 							MAKE MEDIANED IMAGE
	################################################################################################

	montage2 = pexpect.spawn('montage2')

	# Set up log file
	fout = file('montage_log.txt','w')
	montage2.logfile = fout

	montage2.expect("File with transformations:")
	montage2.sendline(star_name + '_' + wavelength + '_f2.mch')
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
	montage2.sendline(star_name + '_' + wavelength + '_f2.fits')
	montage2.expect("Good bye")
	montage2.close(force=True)	

	# Now want to check whether a good medianed image was made
	# To do this, examine the weightings in the montage log file 
	# They should all be quite high. If there exists one which isn't then bad image made
	# In this case, the overly-weighted image gets discarded and a new image made
	# This repeats until all the weightings are acceptable

	# Want to keep a copy of the transformation file for all frames as this will be required later on
	# This step is just in case weightings are bad 
	mch_file = star_name + '_' + wavelength + '_f' + field +'.mch'
	shutil.copy(mch_file, star_name + '_' + wavelength + '_f' + field +'_full.mch') 

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

			bad_frame_value = max(weights) # value of the bad frame, usually 1000.00

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

			# Re-run MONTAGE2
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

	################################################################################################
	# 							X AND Y OFFSETS
	################################################################################################

	# Write down X and Y offsets - these need to be added back into master star list once it's been created
	log = open('montage_log.txt', 'r')
	lines = log.readlines()

	offsets = []
	
	for line in lines:
		if "Offsets" in line:

			offsets.append(line.split(' ')[-3])
			offsets.append(line.split(' ')[-2])


	# Copy relevant options files across to cwd
	shutil.copy('/home/ac833/daophot-options-files/daophot.opt', 'daophot.opt')
	shutil.copy('/home/ac833/daophot-options-files/photo.opt', 'photo.opt')
	shutil.copy('/home/ac833/daophot-options-files/allstar.opt', 'allstar.opt')


	################################################################################################
	# 							MAKE MASTER STAR LIST
	################################################################################################

	daophot = pexpect.spawn('daophot') 

	# Set up logfile
	fout = file('daophot_log.txt','w')
	daophot.logfile = fout

	# Attach medianed image and obtain star list
	daophot.expect("Command:")
	daophot.sendline("at " + star_name + '_' + wavelength + '_f' + field +'.fits')
	daophot.expect("Command:")
	daophot.sendline("opt")
	daophot.expect("File with parameters")
	daophot.sendline("")
	daophot.expect("OPT>")
	daophot.sendline("th=20") # set an appropriately high threshold for this highly S/N medianed image
	daophot.expect("OPT>")
	daophot.sendline("")	

	daophot.expect("Command:")
	daophot.sendline("fi")
	daophot.expect("Number of frames averaged, summed:")
	daophot.sendline(str(num_images)+",1") 
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

	################################################################################################
	# 							MAKE PSF MODEL
	################################################################################################

	# Choose 20 brightest stars as candidate PSF stars
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

	# Now run these candidate PSF stars through series of tests to get rid of bad stars

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

	psf_stars = pd.read_csv(star_name + '_' + wavelength + '_f' + field + '.lst', delim_whitespace=True, skiprows=3, header=None, names=['ID', 'X', 'Y', 'Mag', 'Error'], index_col=0)

	deleted_stars = 0 

	# Carry out all the tests on each star in the df 'psf_stars'
	for index, row in psf_stars.iterrows():

		execute = 1

		# TEST 1 : TOO CLOSE TO EDGE OF FRAME

		# If X < x_lo or X > x_up, drop row
		if row['X'] < x_lo or row['X'] > x_up:
			psf_stars.drop(index, inplace=True)
			deleted_stars += 1
			#print "Deleting star %d because it is too close to edge of frame" % index
			execute = 0 # don't need to carry out rest of tests

		if execute == 1:

			# If Y < y_lo or Y > y_up, drop row
			if row['Y'] < y_lo or row['Y'] > y_up:
				psf_stars.drop(index, inplace=True)
				deleted_stars += 1
				#print "Deleting star %d because it is too close to edge of frame" % index
				execute = 0 # don't need to carry out rest of tests

		# TEST 2 : NOT BRIGHT ENOUGH

		# Get x and y coords of the star in question
		x_coord = int(round(row['X'] - 1)) # zero-indexed in data and must be rounded to nearest integer
		y_coord = int(round(row['Y'] - 1)) # zero-indexed in data and must be rounded to nearest integer

		if execute == 1:
			if data[y_coord, x_coord] < 150:
				psf_stars.drop(index, inplace=True)
				deleted_stars += 1
				#print "Deleting star %d because it is not bright enough" % index
				execute = 0 # don't need to carry out rest of tests

	print "Stars remaining: " + str(20-deleted_stars)

	# Write out final list of stars to the lst file in the correct format

	# Get header of lst file
	f = open(star_name + '_' + wavelength + '_f' + field + '.lst', 'r')
	header = f.read().splitlines()[0:3]
	f.close()

	# Now overwrite this file
	f = open(star_name + '_' + wavelength + '_f' + field + '.lst', 'w')
	f.writelines(header[0] + '\n' + header[1] + '\n' + header[2] + '\n')

	# Send stars to lst file
	psf_stars.to_csv(f, sep=' ', mode='a', header=None)

	f.close()

	# Use these stars to create PSF model and run ALLSTAR

	daophot = pexpect.spawn('daophot')
	daophot.expect("Command:")
	daophot.sendline("at " + star_name + '_' + wavelength + '_f' + field +'.fits')

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

	allstar = pexpect.spawn('allstar')

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

	################################################################################################
	# 						ADD OFFSETS BACK IN TO MASTER STAR LIST
	################################################################################################

	daophot = pexpect.spawn('daophot')

	daophot.expect("Command:")
	daophot.sendline("off") # offsets to put x and y back in
	daophot.expect("Input file name:")
	daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.als')
	daophot.expect("Additive offsets ID, DX, DY, DMAG:")
	daophot.sendline("0," + offsets[0] + "," + offsets[1] + ",0")
	daophot.expect("Output file name")
	daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.mag') # this file is the master star list 
	daophot.expect("Command:")
	daophot.sendline("ex")
	daophot.close(force=True)

	################################################################################################
	# 					TIDYING UP AND COPYING FILES TO INDIVIDUAL EPOCHS
	################################################################################################

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

	return(0)

# Master image off target
def master_off_target(star_name, galaxy, channel, wavelength, epoch_number, num_epochs):

	num_images = 5 # for off-target field, medianed images are made for each epoch
	field = '1'
	start_dither = 1 # for off-target fields, dithers run from 1 - 5

	# Path to where all data is kept
	home = '/home/ac833/Data/'	

	# # For each epoch, create medianed image of the 5 dithers for field 1
	# for epoch in range(1, num_epochs+1):

	# 	# Get epoch number into correct format
	# 	if epoch < 10:
	# 		epoch = '0' + str(epoch)
	# 	else: epoch = str(epoch)

	# Change directory to where images are
	cwd = home + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number + '/'
	os.chdir(cwd)

	# Copy options files to cwd
	shutil.copy('/home/ac833/daophot-options-files/daophot.opt', 'daophot.opt')
	shutil.copy('/home/ac833/daophot-options-files/photo.opt', 'photo.opt')
	shutil.copy('/home/ac833/daophot-options-files/allstar.opt', 'allstar.opt')

	################################################################################################
	# 							MAKE FILE OF COORDINATE TRANSFORMATIONS
	################################################################################################
	
	daomatch = pexpect.spawn('daomatch')

	# Set up log file
	fout = file('daomatch_log.txt','w')
	daomatch.logfile = fout

	# Give DAOMATCH all the 5 dithers to be used in making the medianed image
	# The '/' after the file names lets DAOMATCH know that the scale is the same and it is only the rotation and shifts that might have changed
	daomatch.expect("Master input file:")
	daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither) + '_cbcd_dn.ap') # Give it the first BCD phot file
	daomatch.expect("Output file name")
	daomatch.sendline(star_name + '_' + wavelength +'_f1.mch')
	daomatch.expect("Next input file:")
	daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+1) + '_cbcd_dn.ap/') # Give it the second BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+2) + '_cbcd_dn.ap/') # Give it the third BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+3) + '_cbcd_dn.ap/') # Give it the fourth BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+4) + '_cbcd_dn.ap/') # Give it the fifth BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline("") # exit

	# Refine transformations with DAOMASTER
	daomaster = pexpect.spawn('daomaster')

	# Set up log file
	fout = file('daomaster_log.txt','w')
	daomaster.logfile = fout

	daomaster.expect("File with list of input files:")
	daomaster.sendline(star_name + '_' + wavelength + '_f1.mch')
	daomaster.expect("Minimum number, minimum fraction, enough frames:")
	daomaster.sendline("1, 0.5, 5") 
	daomaster.expect("Maximum sigma:")
	daomaster.sendline("99") 
	daomaster.expect("Your choice:")
	daomaster.sendline("6") # solve for 6 degrees of freedom
	daomaster.expect("Critical match-up radius:")
	daomaster.sendline("7") 
	
	for dither in range(start_dither+1,start_dither+5):
		daomaster.expect(star_name + '_' + wavelength + '_e' + epoch_number + '_d' + str(dither) + '_cbcd_dn.ap')
		daomaster.sendline("")

	# Repeat with decreasing match up size
	for match_up in range(7,-1,-1):
		daomaster.expect("New match-up radius")
		daomaster.sendline(str(match_up))

	# Options for different output files; we only require the refined transformations
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
	daomaster.sendline(star_name + '_' + wavelength +'_f1.mch')
	daomaster.expect("New output file name")
	daomaster.sendline("")
	daomaster.expect("A file with the transfer table?")
	daomaster.sendline("e") # exits rest of options

	################################################################################################
	# 								MAKE MEDIANED IMAGE
	################################################################################################

	montage2 = pexpect.spawn('montage2')

	# Set up log file - need this to obtain the offsets in X and Y
	fout = file('montage_log.txt','w')
	montage2.logfile = fout

	montage2.expect("File with transformations:")
	montage2.sendline(star_name + '_' + wavelength +'_f1.mch')
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
	montage2.sendline(star_name + '_' + wavelength +'_f1.fits')


	################################################################################################
	# 									OBTAIN OFFSETS
	################################################################################################

	log = open('montage_log.txt', 'r')
	lines = log.readlines()

	offsets = []
	
	for line in lines:
		if "Offsets" in line:

			offsets.append(line.split(' ')[-3])
			offsets.append(line.split(' ')[-2])		


	################################################################################################
	# 								CREATE MASTER STAR LIST
	################################################################################################

	daophot = pexpect.spawn('daophot') 

	# Set up logfile
	fout = file('daophot_log.txt','w')
	daophot.logfile = fout

	# Attach medianed image
	daophot.expect("Command:")
	daophot.sendline("at " + star_name + '_' + wavelength + '_f1.fits')
	daophot.expect("Command:")
	daophot.sendline("opt")
	daophot.expect("File with parameters")
	daophot.sendline("")
	daophot.expect("OPT>")
	daophot.sendline("th=20") # set to a higher threshold to account for higher S/N 
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
	daophot.sendline(star_name + '_' + wavelength + '_f1.coo')
	daophot.expect("Output file")
	daophot.sendline(star_name + '_' + wavelength + '_f1.ap')


	################################################################################################
	# 									CREATE PSF MODEL
	################################################################################################


	# Choose 20 brightest stars in image to be candidate PSF stars
	daophot.expect("Command:")
	daophot.sendline("pi")
	daophot.expect("Input file name")
	daophot.sendline(star_name + '_' + wavelength + '_f1.ap')
	daophot.expect("Desired number of stars, faintest magnitude:")
	daophot.sendline("20,99")
	daophot.expect("Output file name")
	daophot.sendline(star_name + '_' + wavelength + '_f1.lst') 
	daophot.expect("Command:")
	daophot.sendline("ex")
	daophot.close(force=True)

	# Run these candidate stars through tests to remove bad stars

	hdulist = fits.open(star_name + '_' + wavelength + '_f1.fits')

	# Access the primary header-data unit (HDU)
	hdu = hdulist[0]
	data = hdu.data

	# Obtain the length of the x and y axis of the image
	x_axis = hdulist[0].header['NAXIS1']
	y_axis = hdulist[0].header['NAXIS2']

	centre = [x_axis/2, y_axis/2] # centre of frame

	# Obtain lower and upper x and y limits for Test 1
	x_lo = centre[0] - (9*centre[0])/10
	x_up = centre[0] + (9*centre[0])/10
	y_lo = centre[1] - (9*centre[1])/10
	y_up = centre[1] + (9*centre[1])/10

	psf_stars = pd.read_csv(star_name + '_' + wavelength + '_f1.lst', delim_whitespace=True, skiprows=3, header=None, names=['ID', 'X', 'Y', 'Mag', 'Error'], index_col=0)

	# Carry out all the tests on each star in the df
	for index, row in psf_stars.iterrows():

		execute = 1

		# TEST 1 : TOO CLOSE TO EDGE OF FRAME

		# If X < x_lo or X > x_up, drop row
		if row['X'] < x_lo or row['X'] > x_up:
			psf_stars.drop(index, inplace=True)
			#print "Deleting star %d because it is too close to edge of frame" % index
			execute = 0 # don't need to carry out rest of tests

		if execute == 1:

			# If Y < y_lo or Y > y_up, drop row
			if row['Y'] < y_lo or row['Y'] > y_up:
				psf_stars.drop(index, inplace=True)
				#print "Deleting star %d because it is too close to edge of frame" % index
				execute = 0 # don't need to carry out rest of tests

		# TEST 2 : NOT BRIGHT ENOUGH

		# Get x and y coords of the star in question
		x_coord = int(round(row['X'] - 1)) # zero-indexed in data and must be rounded to nearest integer
		y_coord = int(round(row['Y'] - 1)) # zero-indexed in data and must be rounded to nearest integer

		if execute == 1:
			if data[y_coord, x_coord] < 100:
				psf_stars.drop(index, inplace=True)
				#print "Deleting star %d because it is not bright enough" % index
				execute = 0 # don't need to carry out rest of tests

		# Write out final list of stars to the lst file in the correct format

		# Get header of lst file
		f = open(star_name + '_' + wavelength + '_f1.lst', 'r')
		header = f.read().splitlines()[0:3]
		f.close()

		# Now overwrite this file
		f = open(star_name + '_' + wavelength + '_f1.lst', 'w')
		f.writelines(header[0] + '\n' + header[1] + '\n' + header[2] + '\n')

		# Send stars to lst file
		psf_stars.to_csv(f, sep=' ', mode='a', header=None)

		f.close()

	# Make PSF model
	daophot = pexpect.spawn('daophot')

	daophot.expect("Command:")
	daophot.sendline("at " + star_name + '_' + wavelength + '_f1.fits')
	daophot.expect("Command:")
	daophot.sendline("psf")
	daophot.expect("File with aperture results")
	daophot.sendline(star_name + '_' + wavelength + '_f1.ap')
	daophot.expect("File with PSF stars")
	daophot.sendline(star_name + '_' + wavelength + '_f1.lst')
	daophot.expect("File for the PSF")
	daophot.sendline(star_name + '_' + wavelength + '_f1.psf')
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
	allstar.sendline(star_name + '_' + wavelength + '_f1.fits')
	allstar.expect("File with the PSF")
	allstar.sendline(star_name + '_' + wavelength + '_f1.psf')
	allstar.expect("Input file")
	allstar.sendline(star_name + '_' + wavelength + '_f1.ap')
	allstar.expect("File for results")
	allstar.sendline(star_name + '_' + wavelength + '_f1.als')
	allstar.expect("Name for subtracted image")
	allstar.sendline(star_name + '_' + wavelength + '_f1_dns.fits')
	allstar.expect("Good bye")
	allstar.close(force=True)	

	################################################################################################
	# 							ADD OFFSETS BACK IN TO STAR LIST
	################################################################################################	

	daophot = pexpect.spawn('daophot')

	daophot.expect("Command:")
	daophot.sendline("off") # offsets to put x and y back in
	daophot.expect("Input file name:")
	daophot.sendline(star_name + '_' + wavelength + '_f1.als')
	daophot.expect("Additive offsets ID, DX, DY, DMAG:")
	daophot.sendline("0," + offsets[0] + "," + offsets[1] + ",0")
	daophot.expect("Output file name")
	daophot.sendline(star_name + '_' + wavelength + '_f1.mag') # this is the final file for star positions that we want 
	daophot.expect("Command:")
	daophot.sendline("ex")
	daophot.close(force=True)

	################################################################################################
	# 							TIDYING UP AND RENAMING FILES 
	################################################################################################	

	# Change .ap file names in .mch file to .als ready for next script
	f = open(star_name +'_' + wavelength + '_f1.mch', 'r')
	filedata = f.read()
	f.close()

	newdata = filedata.replace(".ap",".als")

	f = open(star_name +'_' + wavelength + '_f1.mch','w')
	f.write(newdata)
	f.close()	

	# Rename files we want for the next scripts to be consistent with the on-target version of this code
	os.rename(star_name + '_' + wavelength + '_f1.fits', star_name + '_' + wavelength + '_f1_master.fits')
	os.rename(star_name + '_' + wavelength + '_f1.mch', star_name + '_' + wavelength + '_f1_master.mch')
	os.rename(star_name + '_' + wavelength + '_f1.mag', star_name + '_' + wavelength + '_f1_master.mag')
	os.rename(star_name + '_' + wavelength + '_f1.psf', star_name + '_' + wavelength + '_f1_master.psf')

	return(0)

# PSF photometry
def psf_phot(star_name, galaxy, channel, wavelength, epoch_number):

	for j in range(1,11):

		dither = str(j)

		# Use f1 psf for dithers 1-5 and f2 psf for dithers 6-10
		if j <= 5:
			field = '1'
		else: field = '2'

		# Copy ALLSTAR options file
		shutil.copy('/home/ac833/daophot-options-files/allstar.opt', 'allstar.opt')

		# Spawn ALLSTAR and perform PSF photometry using PSF model made from master image
		allstar = pexpect.spawn('allstar')

		fout = file('allstar_log.txt', 'w')
		allstar.logfile = fout

		allstar.expect("OPT>")
		allstar.sendline("")
		allstar.expect("Input image name:")
		allstar.sendline(star_name + '_' + wavelength + '_e'+ epoch_number + '_d' + str(dither) + '_cbcd_dn.fits')
		allstar.expect("File with the PSF")
		allstar.sendline(star_name + '_' + wavelength + '_f' + str(field) + '_master.psf')

		allstar.expect("Input file")
		allstar.sendline(star_name + '_' + wavelength + '_e'+ epoch_number + '_d' + str(dither) + '_cbcd_dn.ap')
		allstar.expect("File for results")
		allstar.sendline(star_name + '_' + wavelength + '_e'+ epoch_number + '_d' + str(dither) + '_cbcd_dn.als')
		allstar.expect("Name for subtracted image")
		allstar.sendline(star_name + '_' + wavelength + '_e'+ epoch_number + '_d' + str(dither) + '_cbcd_dns.fits')
		allstar.expect("Good bye")
		allstar.close(force=True)	

	return(0)

# ALLFRAME
def allframe(star_name, galaxy, channel, wavelength, epoch_number, field):

	if field == '1':
		start_dither = 1
	else: start_dither = 6

	epoch_string = 'e' + epoch_number

	################################################################################################
	# 							MODIFY TRANSFORMATIONS FILE
	################################################################################################

	# The transformation file currently has transformations for ALL epochs and dithers
	# We want to remove frames not relevant to the current epoch
	# However, we want to keep e01 d1 (for field 1) or d6 (for field 2) as this is the 'master' input file

	# Open mch file i.e. the transformation file and delete rows not relevant to current epoch

	mch_file = star_name + '_' + wavelength + '_f' + field + '_master.mch'
	mch_df = pd.read_csv(mch_file, delim_whitespace=True, header=None, names=['Filename', 'Apostrophe', 'A', 'B', 'C', 'D', 'E', 'F', 'Mag_offset', 'Scatter'])

	indices_to_remove = []

	for index, row in mch_df.iterrows():

		if index != 0: # we don't get rid of first line as this is the 'master' input file

			# Obtain list of indices that contain files we want to remove as not relevant to current epoch
			if epoch_string not in row['Filename']:
				indices_to_remove.append(index)

	# Drop rows that are not relevant to current epoch
	mch_df.drop(indices_to_remove, inplace=True)

	# Write out to new mch file - overwrites file with ALL epochs ins
	mch_df.to_csv(mch_file, header=None, sep=' ', index=False)


	################################################################################################
	# 							COPY NECESSARY EPOCH 1 FILES
	################################################################################################

	# If epoch_number not '01' then copy across necessary files from epoch 1
	if epoch_number != '01':
		shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e01/' + star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.ap', star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.ap')
		shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e01/' + star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.als', star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.als')
		shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e01/' + star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.psf', star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.psf')
		shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e01/' + star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.fits', star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.fits')


	################################################################################################
	# 					COPY MASTER PSF MODEL TO INDIVIDUAL DITHER NAMES
	################################################################################################

	# This step is required by ALLFRAME as it requires a PSF model for each dither
	for dither in range(start_dither, start_dither+5):
		shutil.copy(star_name + '_' + wavelength + '_f' + field + '_master.psf', star_name + '_' + wavelength + '_' + epoch_string + '_d' + str(dither) + '_cbcd_dn.psf')	


	################################################################################################
	# 									RUN ALLFRAME
	################################################################################################

	# Copy ALLFRAME option file to cwd
	shutil.copy('/home/ac833/daophot-options-files/allframe.opt', 'allframe.opt')	

	allframe = pexpect.spawn('allframe')

	fout = file('allframe_log.txt', 'w')
	allframe.logfile = fout

	allframe.expect("OPT>")
	allframe.sendline("")
	allframe.expect("File with list of images:")
	allframe.sendline(mch_file)
	allframe.expect("File with list of stars")
	allframe.sendline(star_name + '_' + wavelength + '_f' + field + '_master.mag')
	allframe.expect("Good bye")
	allframe.close(force=True)

	# Finally, copy e01 .alf file across as will need it later on
	if epoch_number != '01':
		shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e01/' + star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.alf', star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.alf')

	return(0)

# Calibration procedure onto IRAC Vega system
def calibration_procedure(star_name, galaxy, channel, wavelength, epoch_number, field):

	if field == '1':
		start_dither = 1
	else: start_dither = 6

	################################################################################################
	# 							SET UP LISTS OF FILES TO CALIBRATE
	################################################################################################

	# Set up lists of alf files to loop over later
	alf_files = []

	for dither in range(start_dither, start_dither+5):
		alf_files.append(star_name + '_' + wavelength + '_e' + epoch_number + '_d' + str(dither) + '_cbcd_dn.alf')


	################################################################################################
	# 									APERTURE CORRECTION
	################################################################################################

	# The aperture correction is applied to the PSF magnitudes from ALLFRAME (i.e. the .alf files) 
	# To calculate the aperture correction, the aperture magnitudes (.ap) and PSF magnitudes (.als) for a 
	# bright list of stars is used to calculate the average offset.
	# This average offset is then applied to the PSF magnitudes.

	for alf in alf_files:

		################################################################################################
		# 									SET UP DATAFRAMES
		################################################################################################

		# Get filenames for the two files that need comparing 
		als = alf.replace('.alf', '.als')
		ap = alf.replace('.alf', '.ap')

		# Image name
		image = alf.replace('.alf', '.fits')

		# Bright stars that are used to compare magnitude difference
		lst = alf.replace('.alf', '.lst')

		# Open files into dataframes
		# The lst dataframe contains the aperture magnitudes and ALL stars will be used from this file
		# The als dataframe contains the psf magnitudes and only a SUBSET of these stars with IDs matching those of the lst df will be used
		df_lst = pd.read_csv(lst, skiprows=3, header=None, names=['ID_lst', 'X_lst', 'Y_lst', 'ap_mag', 'ap_error'], usecols=(0,1,2,3,4), delim_whitespace=True)
		df_als = pd.read_csv(als, skiprows=3, header=None, names=['ID_als', 'X_als', 'Y_als', 'psf_mag', 'psf_error'], usecols=(0,1,2,3,4), delim_whitespace=True)

		# Sort rows by ID number
		df_lst.sort_values('ID_lst', axis=0, inplace=True)
		df_als.sort_values('ID_als', axis=0, inplace=True)

		# Reset index to start at 0
		df_lst.reset_index(drop=True, inplace=True)
		df_als.reset_index(drop=True, inplace=True)	

		# Crossmatch the two dataframes to only keep stars that appear in both dataframes
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


		################################################################################################
		# 							CALCULATE APERTURE CORRECTION
		################################################################################################

		# Calculate the difference between aperture and PSF magnitudes
		df_combined['mag_difference'] = df_combined['ap_mag'] - df_combined['psf_mag']	

		# Calculate aperture correction value by taking average
		apc = round(df_combined['mag_difference'].mean(), 3)


		################################################################################################
		# 					APPLY APERTURE CORRECTION TO ALLFRAME MAGNITUDES
		################################################################################################

		# This step converts the PSF magnitudes onto the (3,12,20) aperture system

		df_alf = pd.read_csv(alf, header=None, skiprows=3, usecols=[0,1,2,3,4], names=['ID', 'x', 'y', 'mag', 'error'], delim_whitespace=True)
		df_alf['mag'] += apc # ADD apc value to magnitudes


		################################################################################################
		# 						WRITE CORRECTED MAGNITUDES TO NEW FILE
		################################################################################################

		new_filename = alf.replace('.alf', '.alf_apc')
		df_alf.to_csv(new_filename, sep=' ', index=False)


	################################################################################################
	# 							STANDARD APERTURE AND ZERO POINT
	################################################################################################

	# This correction takes the (3,12,20) magnitudes and puts them onto the standard (10,12,20) 
	# system used by Reach et al. (2005) and is the standard IRAC Vega system.
	# It also corrects for the zero point i.e. the flux of a m=0 star.

	for alf in alf_files:

		# Working on the .alf_apc magnitudes which have been aperture corrected already
		alf = alf.replace('.alf', '.alf_apc')

		################################################################################################
		# 								CALCULATE ZMAG VALUE
		################################################################################################

		# Calculate zmag value
		image = alf.replace('.alf_apc', '.fits')
		fits_image = fits.open(image)
		hdr = fits_image[0].header
		fluxconv = hdr['fluxconv'] * (10 ** 6) 
		px_ste = 3.3845 * (10 ** -11) # px size in steradians 

		# Zero magnitude flux values different for each channel
		# Values are taken from Spitzer IRAC handbook
		if channel == '1':
			F0 = 280.9
		elif channel == '2':
			F0 = 179.7
		else: F0 = 'Invalid'

		zmag = round(2.5 * log10(F0/(fluxconv*px_ste)), 2)

		# To go from (3,12,20) to (10,12,20) we have a correction value from IRAC handbook
		# For this aperture and annulus size, it is the same value for [3.6] and [4.5]
		std_corr = 1.112

		# Put alf_apc magnitudes into dataframe
		df_alf = pd.read_csv(alf, header=0, delim_whitespace=True)

		################################################################################################
		# 								APPLY CORRECTIONS
		################################################################################################

		for index in range(0, len(df_alf)):

			# Convert magnitude to flux
			flux =  F0 * (10 ** (-df_alf['mag'][index]/2.5))

			# Apply std aperture correction value
			flux *= std_corr

			# Convert back to magnitude
			df_alf.ix[index, 'mag'] = -2.5 * log10(flux/F0)

			# Round magnitude to 4 dp
			df_alf.ix[index, 'mag'] = round(df_alf['mag'][index], 4)

		# Apply zmag to magnitudes and subtract the zp of 25 that was used by DAOPHOT
		df_alf['mag'] = df_alf['mag'] + zmag - 25

		################################################################################################
		# 						WRITE CORRECTED MAGNITUDES TO NEW FILE
		################################################################################################

		# Write out to file with ext .alf_zp 
		new_filename = alf.replace('.alf_apc', '.alf_zp')
		df_alf.to_csv(new_filename, index=False, sep=' ')


	################################################################################################
	# 									LOCATION CORRECTION
	################################################################################################

	# This correction accounts for the fact that where the star is in the image affects the flux recorded
	# This step corrects for this using the correction images provided by Spitzer handbook
	# The flux of the star is multiplied by this correction value based on where the star is

	################################################################################################
	# 								  LOAD CORRECTION IMAGE
	################################################################################################

	# Open correction file corresponding to correct channel
	if channel == '1':
		corr_file = '/home/ac833/Spitzer-corr-images/ch1_al_s192.fits'
		corr_image = fits.open(corr_file)
		corr_data = corr_image[0].data
	elif channel == '2':
		corr_file = '/home/ac833/Spitzer-corr-images/ch2_al_s192.fits'
		corr_image = fits.open(corr_file)
		corr_data = corr_image[0].data
	else: corr_file = 'Corr file not found'	

	# Loop over all '.alf' files and correct
	for alf in alf_files:	

		# Get files that have been aperture corrected, std aperture and zmag corrected
		alf = alf.replace('.alf', '.alf_zp')

		# Load .alf_zp file into df
		df = pd.read_csv(alf, header=0, delim_whitespace=True)	

		################################################################################################
		# 								 APPLY LOCATION CORRECTION
		################################################################################################

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


		################################################################################################
		# 								  WRITE OUT TO NEW FILE
		################################################################################################

		# Make new filename by replacing .alf_zp with .alf_loc
		alf = alf.replace('.alf_zp', '.alf_loc')
		df.to_csv(alf, sep=' ', index=False)


	################################################################################################
	# 								PIXEL PHASE CORRECTION
	################################################################################################

	################################################################################################
	# 									SET UP PARAMETERS
	################################################################################################

	# Set parameters for aperture = 3, inner sky = 12, outer sky = 20
	if channel == '1':
		x0 = 0.190
		y0 = 0.022
		sig_x = 0.189
		sig_y = 0.164
		delf_x = 0.0327
		delf_y = 0.0544
		f0 = 0.962
	elif channel == '2':
		x0 = 0.070
		y0 = 0.098
		sig_x = 0.203
		sig_y = 0.216
		delf_x = 0.0175
		delf_y = 0.0186
		f0 = 0.981
	else: x0 = 'break'		

	################################################################################################
	# 								APPLY PIXEL PHASE CORRECTION
	################################################################################################

	# Loop over all '.alf' files and correct
	for alf in alf_files:	

		# Get files that have been aperture corrected, std aperture, zmag corrected and location corrected
		alf = alf.replace('.alf', '.alf_loc')

		# Load .alf_zp file into df
		df = pd.read_csv(alf, header=0, delim_whitespace=True)

		# Loop over every star in file
		for index in range(0, len(df)):

			if df['mag'][index] != 99.999:

				# Get x and y coords
				x = float(df['x'][index])
				y = float(df['y'][index])

				# Convert magnitude to flux
				flux = 10 ** (df['mag'][index]/-2.5)

				# Convert x and y coords of star to a phase
				x_phase = x - int(x) - 0.5
				y_phase = y - int(y) - 0.5

				# Calculate distances dx and dy between centre of star and most responsive part of pixel
				dx = x_phase - x0
				dy = y_phase - y0

				# Compute relative flux
				relative_flux = delf_x * exp((-(dx**2))/(2*sig_x**2)) + delf_y * exp(-(dy**2)/(2*sig_y**2)) + f0

				# Compute corrected flux
				corrected_flux = flux / relative_flux

				# Convert corrected flux back to magnitude
				df.loc[index, 'mag'] = -2.5 * log10(corrected_flux)


		################################################################################################
		# 								  WRITE OUT TO NEW FILE
		################################################################################################

		# Make new filename by replacing .alf_loc with .alf_cal
		# This is the final correction
		alf = alf.replace('.alf_loc', '.alf_cal')
		df.to_csv(alf, sep=' ', index=False)

	return(0)

# Combine magnitudes for each epoch
def combine_dithers(star_name, galaxy, channel, wavelength, epoch_number, field):

	if field == '1':
		start_dither = 1
	else: start_dither = 6

	################################################################################################
	# 								FORMAT .ALF_CAL FILES 
	################################################################################################

	# This step takes the calibrated PSF mag files (.alf_cal) and adds the header to the top of the file
	# Output are .alf_all files
	# This step is required so that the files can be run through DAOMATCH/DAOMASTER

	# Format all 5 dithers for this star, epoch and field
	for i in range(start_dither, start_dither+5):

		# Get filenames
		alf = star_name + '_' + wavelength + '_e' + epoch_number + '_d' + str(i) + '_cbcd_dn.alf' # has the header
		alf_apc = star_name + '_' + wavelength + '_e' + epoch_number + '_d' + str(i) + '_cbcd_dn.alf_cal' # has the calibrated magnitudes

		f = open(alf)

		# Write the 3 line header to a new file
		new_filename = star_name + '_' + wavelength + '_e' + epoch_number + '_d' + str(i) + '_cbcd_dn.alf_all'
		header = f.read().splitlines()[0:3]
		g = open(new_filename, 'w')
		g.writelines(header[0] + '\n' + header[1] + '\n' + header[2] + '\n')

		# Open alf and alf_apc
		df = pd.read_csv(alf, delim_whitespace=True, header=None, skiprows=3, names=['ID', 'X', 'Y', 'Uncorr_Mag', 'Error', 'Sky', 'Iters', 'Chi', 'Sharp'])
		df2 = pd.read_csv(alf_apc, skiprows=1, delim_whitespace=True, header=None, names=['ID_nn', 'X_nn', 'Y_nn', 'mag', 'error'])

		# Concat the two df's and drop unnecessary columns
		data = pd.concat((df,df2), axis=1)
		data.drop(['Uncorr_Mag', 'ID_nn','X_nn', 'Y_nn', 'error'], axis=1, inplace=True)
		data = data[['ID', 'X', 'Y', 'mag', 'Error', 'Sky', 'Iters', 'Chi', 'Sharp']]

		g.close()
		f.close()

		# Append to the new file with the header at the top
		data.to_csv(new_filename, mode='a', header=None, sep=' ', index=False)	


	################################################################################################
	# 							MATCH APERTURE MAGNITUDE FILES 
	################################################################################################

	# The next step is to use DAOMATCH to get transformations for the 5 dithers that make up the field
	# This is differenf from the .mch file we already have as the current file includes the EPOCH 1 FILE
	# We do not want this when we create our final magnitudes file for this epoch (unless of course it is epoch 1)

	daomatch = pexpect.spawn('daomatch')

	fout = file('ap_log.txt','w')
	daomatch.logfile = fout

	# Input first dither of field as master input file
	daomatch.expect("Master input file:")
	daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither) + '_cbcd_dn.ap') 
	daomatch.expect("Output file name")
	daomatch.sendline(star_name + '_f' + field + '_ap.mch')
	daomatch.expect("Next input file:")
	daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number+'_d'+str(start_dither + 1)+'_cbcd_dn.ap/') 
	daomatch.expect("Next input file")
	daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number+'_d'+str(start_dither + 2)+'_cbcd_dn.ap/')
	daomatch.expect("Next input file")
	daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number+'_d'+str(start_dither + 3)+'_cbcd_dn.ap/') 
	daomatch.expect("Next input file")
	daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number+'_d'+str(start_dither + 4)+'_cbcd_dn.ap/') 
	daomatch.expect("Next input file")
	daomatch.sendline("") # exit
	daomatch.expect("Good bye.")
	daomatch.close(force=True)	


	################################################################################################
	# 							OBTAIN FINAL MAGNITUDES FILE
	################################################################################################

	# The final step is to use DAOMASTER with the _ap.mch file of coordinate transformations to actually make final file
	# Output is .mag file

	# File with the coordinate transformations in
	match_file = star_name + '_f' + field + '_ap.mch'

	# Change the file extensions in the .mch file from .ap to .alf_all as we now want to match the PSF mags
	f = open(match_file, 'r')
	filedata = f.read()
	f.close()

	newdata = filedata.replace(".ap",".alf_all")

	f = open(match_file,'w')
	f.write(newdata)
	f.close()	

	# Check alf_all files actually contain at least one star.
	# If one doesn't, then add a fake star just so DAOMASTER will run
	for dither in range(start_dither, start_dither+5):

		num_lines = sum(1 for line in open(star_name + '_' + wavelength + '_e' + epoch_number + '_d' + str(dither) + '_cbcd_dn.alf_all'))

		if num_lines <= 3:
			# Open file and append a fake star to it
			with open(star_name + '_' + wavelength + '_e' + epoch_number + '_d' + str(dither) + '_cbcd_dn.alf_all', "a") as myfile:
				myfile.write("1 100.00 100.00 10.000 1.000 0.01 -0.1 5.0 1.00 0.01") # fake data that shouldn't be matched to anything	

	# Run DAOMASTER
	daomaster = pexpect.spawn('daomaster')

	fout = file('daomaster_log.txt', 'w')
	daomaster.logfile = fout

	daomaster.expect("File with list of input files:")
	daomaster.sendline(match_file)
	daomaster.expect("Minimum number, minimum fraction, enough frames:")
	daomaster.sendline("1, 0.5, 5") # play around with these values
	daomaster.expect("Maximum sigma:")
	daomaster.sendline("0.5") # play around with this value
	daomaster.expect("Your choice:")
	daomaster.sendline("6") # solve for 6 degrees of freedom
	daomaster.expect("Critical match-up radius:")
	daomaster.sendline("7") # play around with this

	for dither in range(start_dither,start_dither+5):
		daomaster.sendline("")

	for dither in range(start_dither+1,start_dither+5):
		daomaster.expect(star_name + '_' + wavelength + '_e' + epoch_number +'_d'+str(dither)+'_cbcd_dn.alf_all')
		daomaster.sendline("")

	# Repeat with decreasing match up size
	for match_up in range(7,-1,-1):
		daomaster.expect("New match-up radius")
		daomaster.sendline(str(match_up))

	# Options for different output files - only want transformations according to cookbook
	daomaster.expect("Assign new star IDs?")
	daomaster.sendline("y") # assign new ids so all frames have same ids - already done previously
	daomaster.expect("A file with mean magnitudes and scatter?")
	daomaster.sendline("n")
	daomaster.expect("A file with corrected magnitudes and errors?")
	daomaster.sendline("n")
	daomaster.expect("A file with raw magnitudes and errors?")
	daomaster.sendline("y")
	daomaster.expect("Output file name")
	daomaster.sendline(star_name + '_' + wavelength + '_f' + str(field) + '.cal')
	daomaster.expect("A file with the new transformations?")
	daomaster.sendline("e") # exits rest of options

	# Remove header from .cal file as don't need it anymore
	cal_data = pd.read_csv(star_name + '_' + wavelength + '_f' + field + '.cal', skiprows=3, header=None, delim_whitespace=True, names=['ID', 'X', 'Y', 'M1', 'E1', 'M2', 'E2', 'M3', 'E3', 'M4', 'E4', 'M5', 'E5'], usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12])
	cal_data.to_csv(star_name + '_' + wavelength + '_f' + field + '.cal', sep=' ', index=False)

	return(0)

# Average magnitudes
def ave_mag(star_name, galaxy, channel, wavelength, epoch_number, field):

	################################################################################################
	# 								SET UP OF DATAFRAME 
	################################################################################################

	# Import data into df
	filename = star_name + '_' + wavelength + '_f' + field + '.cal'
	data = pd.read_csv(filename, header=0, delim_whitespace=True)

	# Get right f0 for the channel
	if wavelength == '3p6um':
		f0 = 280.9
	if wavelength == '4p5um':
		f0 = 179.7		

	# New columns for fluxes
	data['F1'] = 0
	data['F2'] = 0
	data['F3'] = 0
	data['F4'] = 0
	data['F5'] = 0

	# Create new columns in dataframe for the average magnitude, average error and standard error on the mean for each star
	data['f_ave'] = 0
	data['m_ave'] = 0
	data['e_ave'] = 0
	data['std_err'] = 0	


	################################################################################################
	# 							CONVERT MAGNITUDES TO FLUXES 
	################################################################################################

	# Loop over all rows of dataframe to convert any non 99.9999 mags to a flux
	for index in range(0, len(data)):
		if data['M1'][index] != 99.9999:
			f = f0 * (10 ** (-data['M1'][index]/2.5))
			data.loc[index, 'F1'] = f
		if data['M2'][index] != 99.9999:
			f = f0 * (10 ** (-data['M2'][index]/2.5))
			data.loc[index, 'F2'] = f
		if data['M3'][index] != 99.9999:
			f = f0 * (10 ** (-data['M3'][index]/2.5))
			data.loc[index, 'F3'] = f
		if data['M4'][index] != 99.9999:
			f = f0 * (10 ** (-data['M4'][index]/2.5))
			data.loc[index, 'F4'] = f
		if data['M5'][index] != 99.9999:
			f = f0 * (10 ** (-data['M5'][index]/2.5))
			data.loc[index, 'F5'] = f


	################################################################################################
	# 						COMPUTE AVERAGE MAGNITUDES AND ERRORS 
	################################################################################################

	for index in range(0, len(data)):
		total = 0
		error = 0
		count = 0

		if data['M1'][index] != 99.9999:
			total += data['F1'][index]
			error += data['E1'][index]
			count += 1
		if data['M2'][index] != 99.9999:
			total += data['F2'][index]
			error += data['E2'][index]
			count += 1
		if data['M3'][index] != 99.9999:
			total += data['F3'][index]
			error += data['E3'][index]
			count += 1
		if data['M4'][index] != 99.9999:
			total += data['F4'][index]
			error += data['E4'][index]
			count += 1
		if data['M5'][index] != 99.9999:
			total += data['F5'][index]
			error += data['E5'][index]
			count += 1

		if count != 0:

			data.loc[index, 'f_ave'] = total/count
			data.loc[index, 'e_ave'] = round(error/count, 4)

	# Convert average flux to an average magnitude and put in m_ave column
	for index in range(0, len(data)):

		if data['f_ave'][index] != 0:

			m = -2.5 * log10(data['f_ave'][index]/f0)
			data.loc[index, 'm_ave'] = round(m, 4)


	# Calculate standard error on the mean for each star
	for index in range(0, len(data)):

		mags = []

		if data['M1'][index] != 99.9999:
			mags.append(data['M1'][index])
		if data['M2'][index] != 99.9999:
			mags.append(data['M2'][index])
		if data['M3'][index] != 99.9999:
			mags.append(data['M3'][index])
		if data['M4'][index] != 99.9999:
			mags.append(data['M4'][index])
		if data['M5'][index] != 99.9999:
			mags.append(data['M5'][index])

		data.loc[index, 'std_err'] = round(np.std(mags)/sqrt(len(mags)),4)


	################################################################################################
	# 						WRITE AVERAGE MAGNITUDES TO FILE 
	################################################################################################

	filename = filename.replace('.cal', '.ave')
	data.to_csv(filename, columns=['ID', 'X', 'Y', 'm_ave', 'e_ave', 'std_err'], header=True, sep=" ", index=False)

	return(0)

# Get Cepheid magnitude 
def get_mag(star_name, galaxy, channel, wavelength, num_epochs, ra, dec):

	################################################################################################
	# 								SET UP OF MAGNITUDE FILE 
	################################################################################################

	# Open file to write mag and error to 
	filename = '/home/ac833/Magnitudes/' + galaxy + '/' + star_name + '_' + wavelength +'.txt'
	f = open(filename, 'w')
	f.write("Epoch   ID   X   Y   Mag   Error   Std_err \n") # THIS NEEDS MODIFYING I THINK

	count = 0 # counts how many epochs star was found in --> should equal num_epochs

	# Cepheid is in field 2 files for ch1 and field 1 files for ch2
	if channel == '1':
		field = '2'
		start_dither = 6
	else: 
		field = '1'
		start_dither = 1

	################################################################################################
	# 						FIND CEPHEID AT EACH EPOCH AND WRITE TO FILE
	################################################################################################	

	for epoch in range(1,num_epochs+1):

		# Get epoch in correct form
		if epoch < 10:
			epoch = '0' + str(epoch)
		else: epoch = str(epoch)

		# Change cwd to folder with data in
		cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch + '/'
		os.chdir(cwd)

		# Convert RA and Dec to pixel coordinates using astropy.wcs
		hdulist = fits.open(star_name + '_' + wavelength + '_e' + epoch + '_d' + str(start_dither) + '_cbcd_dn.fits')		
		w = wcs.WCS(hdulist[0].header) # Parse the WCS keywords in the primary HDU
		world = np.array([[ra, dec]], np.float_) # the world coordinates of the target Cepheid
		pix = w.wcs_world2pix(world,1) # convert world coordinates to pixel coordinates in an array. For just list of coords want pix[0]

		# Coordinates we want to compare .ave file with
		pix_coord = pix[0] # just as a list. x = pix_coord[0] and y = pix_coord[1]

		# Open file containing average magnitudes
		ave_file = star_name + '_' + wavelength + '_f' + field + '.ave'

		if (os.path.isfile(ave_file)):

			ave = pd.read_csv(ave_file, header=0, delim_whitespace=True)

			for i in range(0, len(ave)):

				# Find stars that could potentially be the target star because they are in the right area of the frame 
				if abs(ave['X'][i]-pix_coord[0]) <= 3 and abs(ave['Y'][i]-pix_coord[1]) <= 3:

					f.writelines("%s %f %f %f %f %f %f \n" % (epoch, float(ave['ID'][i]), float(ave['X'][i]), float(ave['Y'][i]), float(ave['m_ave'][i]), float(ave['e_ave'][i]), float(ave['std_err'][i])))
					count += 1

	print "Epochs found = " + str(count)

	f.close()

	return(count)

# Format to GLOESS-style file
def format_gloess(star_name, galaxy, channel, wavelength, num_epochs, period):

	################################################################################################
	# 								GET HMJD FOR EACH EPOCH
	################################################################################################

	hmjd = []

	for j in range(1,num_epochs+1):

		if j < 10:
			filename = '/home/ac833/Data/'+str(galaxy)+'/BCD/'+str(star_name)+'/ch'+str(channel)+'/e0'+str(j)+'/'+str(star_name)+'_'+str(wavelength)+'_e0'+str(j)+"_d1_cbcd_dn.fits"
		else: filename = '/home/ac833/Data/'+str(galaxy)+'/BCD/'+str(star_name)+'/ch'+str(channel)+'/e'+str(j)+'/'+str(star_name)+'_'+str(wavelength)+'_e'+str(j)+"_d1_cbcd_dn.fits"

		image = fits.open(filename)
		header = image[0].header
		hmjd.append(header['HMJD_OBS'])


	################################################################################################
	# 								GET MAGNITUDES AND ERRORS
	################################################################################################

	mag_file = '/home/ac833/Magnitudes/'+str(galaxy)+'/'+str(star_name)+'_'+wavelength+'.txt'
	mags = pd.read_csv(mag_file, header=0, delim_whitespace=True)
	magnitudes = list(mags['Mag'])
	errors = list(mags['Std_err'])


	################################################################################################
	# 									WRITE OUT NEW FILE
	################################################################################################

	output_filename = '/home/ac833/GLOESS_files/' + str(star_name) + '_' + wavelength + '_gloess.txt'

	f = open(output_filename, 'w')
	f.write(str(star_name) + '\n')
	f.write(str(period) + '\n')
	f.write(str(1) + '\n') # doesn't matter what goes here as it's not used by GLOESS code

	# If short period then slightly relaxed fitting parameter of 0.25
	# Longer period Cepheids generally have better light curves so can use a tighter fitting parameter of 0.20
	if period <= 12:
		f.write(str(0.25) + '\n') 
	else: f.write(str(0.2) + '\n')

	for x in range(0, num_epochs):
		f.write(str(hmjd[x]) + ' ' + str(magnitudes[x]) + ' ' + str(errors[x]) + '\n')

	f.close()

	return(0)

# Fit GLOESS
def gloess_single_band(star_name, galaxy, channel, wavelength):

	################################################################################################
	# 									READ IN DATA
	################################################################################################

	# This file has all the magnitudes
	input_filename = '/home/ac833/GLOESS_files/' + star_name + '_' + wavelength + '_gloess.txt'

	dir1 = [] # magnitudes
	deir1 = [] # errors
	dmjd = [] # hmjd of observations

	counter = 0

	## Want to know whether the IRAC data is phased or not. 
	## If it is phased, must reduce the uncertainty by another factor of sqrt(N)
	## if phased == 1 then true. if phased == 0, false
	## The data is phased if Period > 12 days

	for line in open(input_filename):

		data = line.split()

		# Counter always starts equal to 0 so this will always be the first line
		if counter == 0:
			# Cephied name	
			cepname = data[0]
		if counter == 1:
			# if counter = 1 then the first element in this line of the input file is the period of the cepheid
			period = float(data[0])
			if period > 12:
				# for Cepheids with periods > 12 days then the observations were taken equally spaced
				# therefore, error scales with 1/N
				phased = 1
			else:
				# for Cepheids with periods < 12 days they weren't equally spaced
				# therefore, error scales with 1/sqrt(N)
				phased = 0
		if counter == 2:
			# if counter = 2, the first element in this line of the input file is the number of lines
			# this part of the code is redundant since don't use num_lines anymore
			nlines = float(data[0])
		if counter == 3:
			# if counter = 3, this line of data contains the x values in all the 12 bands
			# only doing one band here, so will only be one number
			xir1 = float(data[0])
		if counter > 3:
			# if counter>3, then this line of the input consists of many values that get appended to the lists at the beginning of the code
			# I think the d in front of du say is just so that when convert to numpy arrays the u isn't already taken
			dmjd.append(float(data[0])) # hmjd
			dir1.append(float(data[1])) # magnitude
			deir1.append(float(data[2])) # error in magnitude

		# Counter incrememnts so that each line from the input file is read 
		counter  = counter + 1	
			
	## Read in all the data from the file and filled the arrays. Need to convert these to numpy arrays.

	# Number of lines of data as opposed to the Name, period etc
	number = counter - 4 # Number data lines in the file

	# Convert to numpy arrays
	ir1 = np.array(dir1)
	eir1 = np.array(deir1)
	mjd = np.array(dmjd)	

	# Number of observations
	nir1 = sum(ir1<50)


	################################################################################################
	# 									PHASE THE DATA
	################################################################################################

	# Phases don't need to be done individually by band - only depends on P
	phase = (mjd / period) - np.floor(mjd / period)
	multiple_phase = np.concatenate((phase,(phase+1.0),(phase+2.0),(phase+3.0),(phase+4.0)))


	################################################################################################
	# 								GET MAX AND MIN VALUES
	################################################################################################

	# These steps are more necessary when plotting multiple bands

	maxvals = []
	minvals = []

	if nir1 > 0:
		maxvals.append(np.amax(ir1[ir1<50]))
		minvals.append(np.amin(ir1[ir1<50]))

	# Convert the max and min values lists to np arrays
	maxvals = np.array(maxvals)
	minvals = np.array(minvals)

	# Max and min across all bands are the largest values of the max and min arrays
	max = np.max(maxvals)
	min = np.min(minvals)

	# Set max and min limits
	maxlim = max + 0.5
	minlim = min - 0.5

	#print cepname, ' ---- Period =', period, 'days'
	#print '------------------------------------------------------'

	################################################################################################
	# 								FIT GLOESS CURVE TO DATA
	################################################################################################

	# Clear figure
	plt.clf()

	# GridSpec specifies the geometry of the grid that a subplot will be placed. 
	# The number of rows and number of columns of the grid need to be set.
	# So here, we have 3 rows and 4 columns
	gs = gridspec.GridSpec(3, 4)
	ax1 = plt.subplot(gs[:,:]) # plot just takes up entire grid

	# Values in the [] are xmin,xmax,ymin,ymax
	ax1.axis([1,3.5,(maxlim),(minlim)])

	titlestring = cepname + ', P = ' + str(period) + ' days'
	plt.suptitle(titlestring, fontsize=20)

	ax1.set_ylabel('Magnitude')
	ax1.set_xlabel('Phase $\phi$')

	## Fitting and plotting for each band
	if nir1 > 0:

		ir11, ir1x, yir1, yeir1, xphaseir1 = gf.fit_one_band(ir1,eir1,phase,nir1,xir1)

		ax1.plot(ir1x,ir11,'k-')
	 	ax1.plot(xphaseir1,yir1,color='MediumVioletRed',marker='o',ls='None', label='mag')

		aveir1, adevir1, sdevir1, varir1, skewir1, kurtosisir1, ampir1 = gf.moment(ir11[200:300],100)

		if phased == 1:
			factor = sqrt(nir1)
		if phased == 0:
			factor = 1 
		if nir1 > 1:
			print  '<mag> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f}' .format(aveir1, sdevir1/factor, ampir1)
		if nir1 == 1:
			print  'mag = {0:.3f} --- single point'.format(aveir1)

	handles, labels = ax1.get_legend_handles_labels() 

	################################################################################################
	# 						SAVE LIGHT CURVE AND MAGNITUDE TO FILE
	################################################################################################

	# Save light curve
	lc_name = '/home/ac833/Light_Curves/' + galaxy +'/' + star_name + '_' + wavelength +'.png'
	plt.savefig(lc_name)

	# Output all relevant information to a file
	# Galaxy, star, channel, period, mean mag, std_dev, amplitude

	master_mags_file = '/home/ac833/master_stars.txt'

	f = open(master_mags_file, 'a')

	if phased == 1:
		std_error = sdevir1/(factor ** 2) # error scales with 1/N
	else:
		std_error = sdevir1/factor # error scales with 1/sqrt(N)


	f.write(str(galaxy) + ' ' + str(star_name) + ' ' + str(channel) + ' ' + str(period) + ' ' + str(round(aveir1,2)) + ' ' + str(round(std_error, 4)) + ' ' + str(round(ampir1,4)) + '\n' ) 
	f.close()

	return(0)