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
		daophot.sendline("30,99") # used for the aperture correction later
		daophot.expect("Output file name")
		daophot.sendline(file_stem + '.lst')   

		# Exit daophot
		daophot.expect("Command:")
		daophot.sendline("exit")
		daophot.close(force=True)

	return(0)

# Master image on target
def master_on_target(star_name, galaxy, channel, wavelength, epoch_number, num_epochs):

	if channel == '1':
		field = '2' # the on-target images are the second 5 dithers d6-10
		start_dither = 6
	else:
		field = '1' # the on-target images are the first 5 dithers d1-5
		start_dither = 1

	num_images = num_epochs * 5 # 5 dithers per epoch 

	# List of files to use
	files = []

	for epoch in range(1,num_epochs+1):
		for dither in range(start_dither,start_dither+5): 

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
	daomatch.sendline(files[0]) 
	daomatch.expect("Output file name")
	daomatch.sendline(star_name + '_' + wavelength + '_f' + field + '.mch')	

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
	daomaster.sendline(star_name + '_' + wavelength + '_f' + field + '.mch')
	daomaster.expect("Minimum number, minimum fraction, enough frames:")
	daomaster.sendline("1, 0.5, " + str(num_images)) # play around with these values
	daomaster.expect("Maximum sigma:")
	daomaster.sendline("99") 
	daomaster.expect("Your choice:")
	daomaster.sendline("6") # solve for 6 degrees of freedom
	daomaster.expect("Critical match-up radius:")
	daomaster.sendline("7") 

	for j in range(1,num_images):

		#daomaster.expect(files[j])
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
	daomaster.sendline(star_name + '_' + wavelength + '_f' + field + '.mch') # these are more refined transformations
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

	# # Now run these candidate PSF stars through series of tests to get rid of bad stars

	# # Read in FITS image
	# hdulist = fits.open(star_name + '_' + wavelength + '_f' + field +'.fits')

	# # Access the primary header-data unit (HDU)
	# hdu = hdulist[0]
	# data = hdu.data

	# # Obtain the length of the x and y axis of the image
	# x_axis = hdulist[0].header['NAXIS1']
	# y_axis = hdulist[0].header['NAXIS2']

	# centre = [x_axis/2, y_axis/2] # centre of frame

	# # Obtain lower and upper x and y limits for Test 1
	# x_lo = centre[0] - (3*centre[0])/4
	# x_up = centre[0] + (3*centre[0])/4
	# y_lo = centre[1] - (3*centre[1])/4
	# y_up = centre[1] + (3*centre[1])/4

	# psf_stars = pd.read_csv(star_name + '_' + wavelength + '_f' + field + '.lst', delim_whitespace=True, skiprows=3, header=None, names=['ID', 'X', 'Y', 'Mag', 'Error'], index_col=0)

	# deleted_stars = 0 

	# # Carry out all the tests on each star in the df 'psf_stars'
	# for index, row in psf_stars.iterrows():

	# 	execute = 1

	# 	# TEST 1 : TOO CLOSE TO EDGE OF FRAME

	# 	# If X < x_lo or X > x_up, drop row
	# 	if row['X'] < x_lo or row['X'] > x_up:
	# 		psf_stars.drop(index, inplace=True)
	# 		deleted_stars += 1
	# 		#print "Deleting star %d because it is too close to edge of frame" % index
	# 		execute = 0 # don't need to carry out rest of tests

	# 	if execute == 1:

	# 		# If Y < y_lo or Y > y_up, drop row
	# 		if row['Y'] < y_lo or row['Y'] > y_up:
	# 			psf_stars.drop(index, inplace=True)
	# 			deleted_stars += 1
	# 			#print "Deleting star %d because it is too close to edge of frame" % index
	# 			execute = 0 # don't need to carry out rest of tests

	# 	# TEST 2 : NOT BRIGHT ENOUGH

	# 	# Get x and y coords of the star in question
	# 	x_coord = int(round(row['X'] - 1)) # zero-indexed in data and must be rounded to nearest integer
	# 	y_coord = int(round(row['Y'] - 1)) # zero-indexed in data and must be rounded to nearest integer

	# 	if execute == 1:
	# 		if data[y_coord, x_coord] < 150:
	# 			psf_stars.drop(index, inplace=True)
	# 			deleted_stars += 1
	# 			#print "Deleting star %d because it is not bright enough" % index
	# 			execute = 0 # don't need to carry out rest of tests

	# print "Stars remaining: " + str(20-deleted_stars)

	# # Write out final list of stars to the lst file in the correct format

	# # Get header of lst file
	# f = open(star_name + '_' + wavelength + '_f' + field + '.lst', 'r')
	# header = f.read().splitlines()[0:3]
	# f.close()

	# # Now overwrite this file
	# f = open(star_name + '_' + wavelength + '_f' + field + '.lst', 'w')
	# f.writelines(header[0] + '\n' + header[1] + '\n' + header[2] + '\n')

	# # Send stars to lst file
	# psf_stars.to_csv(f, sep=' ', mode='a', header=None)

	# f.close()

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
	f = open(star_name +'_' + wavelength + '_f' + field + '_full.mch', 'r')
	filedata = f.read()
	f.close()

	newdata = filedata.replace(".ap",".als")

	f = open(star_name +'_' + wavelength + '_f' + field + '_full.mch','w')
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

	num_images = 5 # for off-target field, medianed images are made for each epoch rather than a combined one across all epochs

	if channel == '1':
		field = '1' # the off-target images are the first 5 dithers d1-5
		start_dither = 1
	else:
		field = '2' # the on-target images are the first 5 dithers d6-10
		start_dither = 6

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
	daomatch.sendline(star_name + '_' + wavelength +'_f' + field + '.mch')
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
	daomaster.sendline(star_name + '_' + wavelength + '_f' + field + '.mch')
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
	daomaster.sendline(star_name + '_' + wavelength +'_f' + field + '.mch')
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
	montage2.sendline(star_name + '_' + wavelength +'_f' + field + '.mch')
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
	daophot.sendline("at " + star_name + '_' + wavelength + '_f' + field + '.fits')
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
	daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.coo')
	daophot.expect("Output file")
	daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.ap')


	################################################################################################
	# 									CREATE PSF MODEL
	################################################################################################


	# Choose 20 brightest stars in image to be candidate PSF stars
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

	# Run these candidate stars through tests to remove bad stars

	hdulist = fits.open(star_name + '_' + wavelength + '_f' + field + '.fits')

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

	psf_stars = pd.read_csv(star_name + '_' + wavelength + '_f' + field + '.lst', delim_whitespace=True, skiprows=3, header=None, names=['ID', 'X', 'Y', 'Mag', 'Error'], index_col=0)

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
		f = open(star_name + '_' + wavelength + '_f' + field + '.lst', 'r')
		header = f.read().splitlines()[0:3]
		f.close()

		# Now overwrite this file
		f = open(star_name + '_' + wavelength + '_f' + field + '.lst', 'w')
		f.writelines(header[0] + '\n' + header[1] + '\n' + header[2] + '\n')

		# Send stars to lst file
		psf_stars.to_csv(f, sep=' ', mode='a', header=None)

		f.close()

	# Make PSF model
	daophot = pexpect.spawn('daophot')

	daophot.expect("Command:")
	daophot.sendline("at " + star_name + '_' + wavelength + '_f' + field + '.fits')
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
	allstar.sendline(star_name + '_' + wavelength + '_f' + field + '.fits')
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
	# 							ADD OFFSETS BACK IN TO STAR LIST
	################################################################################################	

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

	################################################################################################
	# 							TIDYING UP AND RENAMING FILES 
	################################################################################################	

	# Change .ap file names in .mch file to .als ready for next script
	f = open(star_name +'_' + wavelength + '_f' + field + '.mch', 'r')
	filedata = f.read()
	f.close()

	newdata = filedata.replace(".ap",".als")

	f = open(star_name +'_' + wavelength + '_f' + field + '.mch','w')
	f.write(newdata)
	f.close()	

	# Rename files we want for the next scripts to be consistent with the on-target version of this code
	os.rename(star_name + '_' + wavelength + '_f' + field + '.fits', star_name + '_' + wavelength + '_f' + field + '_master.fits')
	os.rename(star_name + '_' + wavelength + '_f' + field + '.mch', star_name + '_' + wavelength + '_f' + field + '_master.mch')
	os.rename(star_name + '_' + wavelength + '_f' + field + '.mag', star_name + '_' + wavelength + '_f' + field + '_master.mag')
	os.rename(star_name + '_' + wavelength + '_f' + field + '.psf', star_name + '_' + wavelength + '_f' + field + '_master.psf')

	return(0)

# PSF photometry
def psf_phot(star_name, galaxy, channel, wavelength, epoch_number):

	for j in range(1,11):

		dither = str(j)

		# Copy ALLSTAR options file
		shutil.copy('/home/ac833/daophot-options-files/allstar.opt', 'allstar.opt')

		# Spawn ALLSTAR and perform PSF photometry using PSF model made from master on target image
		allstar = pexpect.spawn('allstar')

		fout = file('allstar_log.txt', 'w')
		allstar.logfile = fout

		allstar.expect("OPT>")
		allstar.sendline("")
		allstar.expect("Input image name:")
		allstar.sendline(star_name + '_' + wavelength + '_e'+ epoch_number + '_d' + str(dither) + '_cbcd_dn.fits')
		allstar.expect("File with the PSF")
		allstar.sendline(star_name + '_on_target_master.psf') # this psf was made from the medianed image of all 120 on target 3p6 and 120 on target 4p5 frames
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

	# This list contains all magnitude differences for a given field and will be averaged over to get aperture correction value
	mag_differences = []

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

		# Delete any rows which have a 99.9999 psf mag
		df_combined.drop(df_combined[df_combined.psf_mag == 99.9999].index, inplace=True)

		if len(df_combined) == 0:
			print "%s - No PSF stars with non-99.9999 values" % alf

		################################################################################################
		# 							CALCULATE APERTURE CORRECTION
		################################################################################################

		# Calculate the difference between aperture and PSF magnitudes
		df_combined['mag_difference'] = df_combined['ap_mag'] - df_combined['psf_mag']	

		# # Calculate aperture correction value by taking average
		# apc = round(df_combined['mag_difference'].mean(), 3)

		# Append every value in df_combined['mag_difference'] to list mag_differences
		for i in range(0,len(df_combined)):

			mag_differences.append(df_combined['mag_difference'][i])



	################################################################################################
	# 							CALCULATE AVERAEGE APERTURE CORRECTION
	################################################################################################

	mag_differences = np.array(mag_differences) # convert to np array

	if len(mag_differences) != 0:

		# Calculate aperture correction value that will be used for this whole field
		apc = round(np.average(mag_differences),3)

	else:

		print "Wavelength = %s  Epoch = %s  Field = %s has no stars for aperture correction." % (str(wavelength), str(epoch_number), str(field))
		apc = 0


	# # Save aperture correction value to file
	# apc_file = open('/home/ac833/apc.txt', 'a')
	# apc_file.write("Aperture correction for Star " + str(star_name) + " Wavelength " + str(wavelength) + " Epoch " + str(epoch_number) + " Field " + str(field) + " is " + str(apc) + '\n')
	# apc_file.close()



	for alf in alf_files:


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
			#std_corr = 1.060
		elif channel == '2':
			F0 = 179.7
			#std_corr = 1.063
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

	# # Set parameters for aperture = 5, inner sky = 5, outer sky = 10
	# if channel == '1':
	# 	x0 = 0.191
	# 	y0 = 0.027
	# 	sig_x = 0.190
	# 	sig_y = 0.163
	# 	delf_x = 0.0298
	# 	delf_y = 0.0479
	# 	f0 = 0.966
	# elif channel == '2':
	# 	x0 = 0.091
	# 	y0 = 0.091
	# 	sig_x = 0.209
	# 	sig_y = 0.209
	# 	delf_x = 0.0176
	# 	delf_y = 0.0183
	# 	f0 = 0.981
	# else: x0 = 'break'	


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

# Fornat the calibrated magnitude files
def format_cal_files(star_name, galaxy, channel, wavelength, epoch_number, field):

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

	return(0)

# Match all on target aperture files 
def match_all_frames(star_name, galaxy, num_epochs):

	num_images = 5 * num_epochs * 2

	#####################################################################################################
	# 									CREATE TEMP FOLDER 
	#####################################################################################################

	# Create temporary folder to work in
	temp = '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/temp/'

	# Delete temp folder if it already exists
	if (os.path.isdir(temp)):
		shutil.rmtree(temp)

	# Make temp folder
	os.mkdir(temp)	

	#####################################################################################################
	# 								CREATE LIST OF FILES TO USE  
	#####################################################################################################

	# List of files to use
	files = []

	# The on target frames for [3.6] are dithers 6 - 10
	for epoch in range(1,num_epochs+1):

		for dither in range(6,11): 

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			dither = str(dither)

			filename = star_name + '_3p6um_e' + epoch + '_d' + dither + '_cbcd_dn.ap'
			files.append(filename)	


	# The on target frames for [4.5] are dithers 1 - 5
	for epoch in range(1,num_epochs+1):
		
		for dither in range(1,6): 

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			dither = str(dither)

			filename = star_name + '_4p5um_e' + epoch + '_d' + dither + '_cbcd_dn.ap'
			files.append(filename)

	#####################################################################################################
	# 								COPY FILES TO TEMP FOLDER  
	#####################################################################################################


	# Copy all FITS images and aperture photometry files to temp folder
	for epoch in range(1,num_epochs+1):

		# Get epoch in correct format
		if epoch < 10:
			epoch = '0' + str(epoch)
		else: epoch = str(epoch)

		for channel in [1,2]:

			if channel == 1:
				wavelength = '3p6um'
				start_dither = 6
			else: 
				wavelength = '4p5um'
				start_dither = 1

			cwd = '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+str(channel)+'/e'+epoch+'/'

			for dither in range(start_dither,start_dither+5):

				dither = str(dither)

				shutil.copyfile(cwd + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.ap', temp + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.ap')
				shutil.copyfile(cwd + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.alf_all', temp + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.alf_all') 
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
	daomatch.sendline(files[0]) 
	daomatch.expect("Output file name")
	daomatch.sendline(star_name + '_on_target.mch')	

	index = 0
	master = 0

	# Give it the rest of the images
	for j in range(1,len(files)):

		#print files[j]

		daomatch.sendline(files[j]+'/')

		index = daomatch.expect(["Next input file", "Write this transformation?"])
		#index = daomatch.expect(["Next input file", "Write this transformation?", pexpect.exceptions.TIMEOUT])


		if index == 1:

			master = 1

			daomatch.sendline('Y')
			daomatch.sendline(files[j]+'/')

		# if index == 2:

		# 	daomatch.sendline('Y')

		# This makes sure the last file is read properly
		if j == len(files)-1:

			daomatch.sendline("")		

		# print "Index = " + str(index)

		# # If index = 0, then we got "Next input file" so matching worked
		# # Do nothing
		# if index == 1:

		# 	print index

		# 	print files[j] + " could not match properly"

		# 	#daomatch.expect('Write this transformation')
		# 	daomatch.sendline('No')

		# # If index = 1, then we got "Write this transformation" so matching bad
		# else:

		# 	print index

		# 	print files[j] + " matched correctly"

		# 	if j == len(files)-1:

		# 		daomatch.sendline("")


		# # Index will tell you whether DAOMATCH could macth the file or not
		# index = daomatch.expect(["Next input file", pexpect.exceptions.TIMEOUT])


		#daomatch.expect("Next input file")
		#daomatch.sendline(files[j]+'/')

	#daomatch.expect("Next input file")
	daomatch.sendline("") # exit

	################################################################################################
	# 							REMOVE DUPLICATE LINES IN MATCH FILE
	################################################################################################

	lines_seen = set() # holds lines already seen
	outfile = open('temp.mch', "w")
	for line in open(star_name+'_on_target.mch', "r"):
	    if line not in lines_seen: # not a duplicate
	        outfile.write(line)
	        lines_seen.add(line)
	outfile.close()

	shutil.copy('temp.mch', star_name + '_on_target.mch')


	################################################################################################
	# 						REPLACE .AP WITH .ALF_ALL IN MCH FILE
	################################################################################################

	contents = ""
	with open(star_name + '_on_target.mch') as f:
		for line in f.readlines():
			contents += line

	newdata = contents.replace(".ap",".alf_all")

	f = open(star_name + '_on_target.mch','w')
	f.write(newdata)
	f.close()	
		

	################################################################################################
	# 						DAOMASTER TO MAKE GIANT FILE OF MAGNITUDES
	################################################################################################

	daomaster = pexpect.spawn('daomaster')

	fout = file('daomaster_log.txt','w')
	daomaster.logfile = fout

	daomaster.expect("File with list of input files:")
	daomaster.sendline(star_name + '_on_target.mch')
	daomaster.expect("Minimum number, minimum fraction, enough frames:")
	daomaster.sendline("1, 0.5, " + str(num_images)) # play around with these values
	daomaster.expect("Maximum sigma:")
	daomaster.sendline("99") 
	daomaster.expect("Your choice:")
	daomaster.sendline("6") # solve for 6 degrees of freedom
	daomaster.expect("Critical match-up radius:")
	daomaster.sendline("7") 

	for j in range(1,len(files)):
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
	daomaster.sendline("y")
	daomaster.expect("Output file name")
	daomaster.sendline(star_name + '_on_target.mag')
	daomaster.expect("A file with the new transformations?")
	daomaster.sendline("y")
	daomaster.expect("Output file name")
	daomaster.sendline(star_name + '_on_target.mch') # these are more refined transformations
	daomaster.expect("New output file name")
	daomaster.sendline("")
	daomaster.expect("A file with the transfer table?")
	daomaster.sendline("e") # exits rest of options

	shutil.copy(star_name + '_on_target.mag', '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/' + star_name + '_on_target.mag')

	# # Move up one folder to star folder
	# os.chdir(os.path.expanduser('../'))
	# shutil.rmtree('temp')

	return(0)

# Average magnitudes in giant magnitudes file
def average_mags(star_name, galaxy, num_epochs):

	os.chdir('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/')

	mag_file = star_name + '_on_target.mag'

	all_df = pd.read_csv(mag_file, skiprows=3, delim_whitespace=True, header=None)

	##########################################################################
	# 					SORT DATAFRAME
	##########################################################################

	# Number of rows of data each star has depends on number of epochs which in turn depends on galaxy
	if galaxy == 'LMC':
		no_rows = 41
	else: no_rows = 21

	# Master dataframe
	df = pd.DataFrame()

	# Put data into dataframe
	for i in range(0, no_rows):

		df_temp = all_df[i::no_rows]
		df_temp.reset_index(drop=True, inplace=True)

		df = pd.concat([df, df_temp], axis=1)

	# Drop rows
	df.dropna(axis=1, inplace=True)

	##########################################################################
	# 					ADD COLUMN NAMES
	##########################################################################

	column_names = ['ID','X','Y']

	for i in range(1,num_epochs+1):

		if i < 10:
			i = '0' + str(i)
		else: i = str(i)

		for j in range(6,11):
			column_names.append('3p6_e' + i + '_d' + str(j) + '_mag')
			column_names.append('3p6_e' + i + '_d' + str(j) + '_err')

	for i in range(1,num_epochs+1):

		if i < 10:
			i = '0' + str(i)
		else: i = str(i)

		for j in range(1,6):
			column_names.append('4p5_e' + i + '_d' + str(j) + '_mag')
			column_names.append('4p5_e' + i + '_d' + str(j) + '_err') 


	# Determine number of Del columns to add
	min_columns = 3 + (num_epochs*2*5*2) # num_epochs x (mag + err) x dithers x channels
	no_del_cols = len(df.columns) - min_columns

	del_cols = []

	for i in range(1, no_del_cols+1):

		column_names.append('Del'+str(i))
		del_cols.append('Del'+str(i))

	# rename column names
	df.columns = column_names

	# Drop 'Del' columns
	df.drop(del_cols, axis=1, inplace=True)	


	# # Each star has 40 rows of data

	# df_new1 = df[0::41]
	# df_new1.reset_index(drop=True, inplace=True)

	# df_new2 = df[1::41]
	# df_new2.reset_index(drop=True, inplace=True)

	# df_new3 = df[2::41]
	# df_new3.reset_index(drop=True, inplace=True)

	# df_new4 = df[3::41]
	# df_new4.reset_index(drop=True, inplace=True)

	# df_new5 = df[4::41]
	# df_new5.reset_index(drop=True, inplace=True)

	# df_new6 = df[5::41]
	# df_new6.reset_index(drop=True, inplace=True)

	# df_new7 = df[6::41]
	# df_new7.reset_index(drop=True, inplace=True)

	# df_new8 = df[7::41]
	# df_new8.reset_index(drop=True, inplace=True)

	# df_new9 = df[8::41]
	# df_new9.reset_index(drop=True, inplace=True)

	# df_new10 = df[9::41]
	# df_new10.reset_index(drop=True, inplace=True)

	# df_new11 = df[10::41]
	# df_new11.reset_index(drop=True, inplace=True)

	# df_new12 = df[11::41]
	# df_new12.reset_index(drop=True, inplace=True)

	# df_new13 = df[12::41]
	# df_new13.reset_index(drop=True, inplace=True)

	# df_new14 = df[13::41]
	# df_new14.reset_index(drop=True, inplace=True)

	# df_new15 = df[14::41]
	# df_new15.reset_index(drop=True, inplace=True)

	# df_new16 = df[15::41]
	# df_new16.reset_index(drop=True, inplace=True)

	# df_new17 = df[16::41]
	# df_new17.reset_index(drop=True, inplace=True)

	# df_new18 = df[17::41]
	# df_new18.reset_index(drop=True, inplace=True)

	# df_new19 = df[18::41]
	# df_new19.reset_index(drop=True, inplace=True)

	# df_new20 = df[19::41]
	# df_new20.reset_index(drop=True, inplace=True)

	# df_new21 = df[20::41]
	# df_new21.reset_index(drop=True, inplace=True)

	# df_new22 = df[21::41]
	# df_new22.reset_index(drop=True, inplace=True)

	# df_new23 = df[22::41]
	# df_new23.reset_index(drop=True, inplace=True)

	# df_new24 = df[23::41]
	# df_new24.reset_index(drop=True, inplace=True)

	# df_new25 = df[24::41]
	# df_new25.reset_index(drop=True, inplace=True)

	# df_new26 = df[25::41]
	# df_new26.reset_index(drop=True, inplace=True)

	# df_new27 = df[26::41]
	# df_new27.reset_index(drop=True, inplace=True)

	# df_new28 = df[27::41]
	# df_new28.reset_index(drop=True, inplace=True)

	# df_new29 = df[28::41]
	# df_new29.reset_index(drop=True, inplace=True)

	# df_new30 = df[29::41]
	# df_new30.reset_index(drop=True, inplace=True)

	# df_new31 = df[30::41]
	# df_new31.reset_index(drop=True, inplace=True)

	# df_new32 = df[31::41]
	# df_new32.reset_index(drop=True, inplace=True)

	# df_new33 = df[32::41]
	# df_new33.reset_index(drop=True, inplace=True)

	# df_new34 = df[33::41]
	# df_new34.reset_index(drop=True, inplace=True)

	# df_new35 = df[34::41]
	# df_new35.reset_index(drop=True, inplace=True)

	# df_new36 = df[35::41]
	# df_new36.reset_index(drop=True, inplace=True)

	# df_new37 = df[36::41]
	# df_new37.reset_index(drop=True, inplace=True)

	# df_new38 = df[37::41]
	# df_new38.reset_index(drop=True, inplace=True)

	# df_new39 = df[38::41]
	# df_new39.reset_index(drop=True, inplace=True)

	# df_new40 = df[39::41]
	# df_new40.reset_index(drop=True, inplace=True)

	# df_new41 = df[40::41]
	# df_new41.reset_index(drop=True, inplace=True)

	# dfs = [df_new1, df_new2, df_new3, df_new4, df_new5, df_new6, df_new7, df_new8, df_new9, df_new10,
	#       df_new11, df_new12, df_new13, df_new14, df_new15, df_new16, df_new17, df_new18, df_new19, df_new20,
	#       df_new21, df_new22, df_new23, df_new24, df_new25, df_new26, df_new27, df_new28, df_new29, df_new30,
	#       df_new31, df_new32, df_new33, df_new34, df_new35, df_new36, df_new37, df_new38, df_new39, df_new40, df_new41]


	# df = pd.concat(dfs, axis=1)
	# df.dropna(axis=1, inplace=True)

	# #######################################################################################################
	# #										ADD COLUMN NAMES
	# #######################################################################################################

	# column_names = ['ID', 'X', 'Y']

	# for i in range(1,25):
	    
	#     if i < 10:
	#         i = '0' + str(i)
	#     else: i = str(i)
	        
	#     for j in range(6,11):
	#         column_names.append('3p6_e' + i + '_d' + str(j) + '_mag')
	#         column_names.append('3p6_e' + i + '_d' + str(j) + '_err')

	# for i in range(1,25):
	    
	#     if i < 10:
	#         i = '0' + str(i)
	#     else: i = str(i)
	        
	#     for j in range(1,6):
	#         column_names.append('4p5_e' + i + '_d' + str(j) + '_mag')
	#         column_names.append('4p5_e' + i + '_d' + str(j) + '_err') 

	
	# # Determine number of Del columns to add
	# no_del_cols = len(df.columns) - 483 # 483 columns always 

	# del_cols = []

	# for i in range(1, no_del_cols+1):

	# 	column_names.append('Del'+str(i))
	# 	del_cols.append('Del'+str(i))

	# # rename column names
	# df.columns = column_names

	# # Drop 'Del' columns
	# df.drop(del_cols, axis=1, inplace=True)

	#######################################################################################################
	#								AVERAGE MAGNITUDES AND ERRORS
	#######################################################################################################	

	# Create empty average mag and average error columns for each wavelength and epoch.
	# Total new cols = 24 epochs x 2 wavelengths x 2 cols = 96
	for epoch in range(1,num_epochs+1):

	    if epoch < 10:
	        epoch_string = 'e0' + str(epoch)
	    else: epoch_string = 'e' + str(epoch)
	    
	    df['3p6_' + epoch_string + '_ave_mag'] = np.nan
	    df['3p6_' + epoch_string + '_ave_err'] = np.nan
	    df['4p5_' + epoch_string + '_ave_mag'] = np.nan
	    df['4p5_' + epoch_string + '_ave_err'] = np.nan

	# Compute intensity-averaged magnitudes for [3.6] and [4.5] at each epoch and for each star
	# So every 10 (5 mags, 5 errs) columns need to be converted to flux and averaged
	# Ignore cols ID, X and Y when averaging

	for channel in [1,2]:
	    
	    if channel == 1:
	        wavelength = '3p6'
	        start_dither = 6
	        F0 = 280.9
	    else: 
	        wavelength = '4p5'
	        start_dither = 1
	        F0 = 179.7
	    
	    for epoch in range(1,num_epochs+1): 
	        
	        if epoch < 10:
	            epoch_string = 'e0' + str(epoch)
	        else: epoch_string = 'e' + str(epoch)
	        
	        # Get list of columns to average over
	        
	        list_of_mags = []
	        list_of_errs = []
	        
	        for dither in range(start_dither, start_dither+5):
	            
	            list_of_mags.append(wavelength + '_' + epoch_string + '_d' + str(dither) + '_mag')
	            list_of_errs.append(wavelength + '_' + epoch_string + '_d' + str(dither) + '_err')
	    
	        
	        # Now convert each magnitude to a flux and average over        
	        for index, row in df.iterrows():
	            
	            fluxes = []
	            errors = []
	            
	            for col in list_of_mags:
	                
	                if df[col][index] != 99.9999:
	                    
	                    # Convert to flux and append to list
	                    flux =  F0 * (10 ** (-float(df[col][index])/2.5))
	                    fluxes.append(flux)
	            
	            # Calculate average flux
	            if len(fluxes) != 0:
	                ave_flux = sum(fluxes)/len(fluxes)
	            else: ave_flux = 0
	            
	            # Convert flux back to magnitude and append to new column
	            if ave_flux != 0:
	                ave_mag = -2.5 * log10(ave_flux/F0)
	            else: ave_mag = 99.9999
	                
	            # Write average magnitude for that star at that epoch to new column
	            new_col = wavelength + '_' + epoch_string + '_ave_mag'
	            df.ix[index, new_col] = ave_mag
	            
	            for col in list_of_errs:
	                
	                if df[col][index] != 9.9999:
	                    
	                    # Convert error to flux and append to list
	                    error = F0 * (10 ** (-float(df[col][index])/2.5))
	                    errors.append(error)
	           
	            # Calculate average error in flux
	            if len(errors) != 0:
	                ave_err_flux = sum(errors)/len(errors)
	            else: ave_err_flux = 0
	            
	            # Convert flux back to magnitude and append to new column
	            if ave_err_flux != 0:
	                ave_err = -2.5 * log10(ave_err_flux/F0)
	            else: ave_err = 9.9999 
	            
	            # Write average magnitude for that star at that epoch to new column
	            new_col = wavelength + '_' + epoch_string + '_ave_err'
	            df.ix[index, new_col] = ave_err 


	#######################################################################################################
	#								DELETE STARS WITHOUT ALL EPOCH DATA
	#######################################################################################################	

	# Delete stars that don't have data for every epoch
	# This is fine as they definitely won't be our V* 

	keep_cols = [col for col in df.columns if 'ave_mag' in col or 'ave_err' in col or 'X' in col or 'Y' in col or 'ID' in col]

	df = df[keep_cols]
	df = df.replace(to_replace=99.9999, value=np.nan)
	df.dropna(inplace=True)

	df.to_csv(star_name + '_ave_mag.txt', index=False)

	return(df)

# Search for Cepheid
def find_target(star_name, galaxy, ra, dec, num_epochs, df):

	################################################################################################
	# 						SET UP OF OUTPUT MAGNITUDE FILE 
	################################################################################################

	# Open file to write mag and error to 
	filename = '/home/ac833/Magnitudes/' + galaxy + '/' + star_name + '.txt'
	f = open(filename, 'w')
	f.write("Epoch   ID   X   Y   3p6_mag  3p6_err  4p5_mag   4p5_err \n") 

	count = 0 # counts how many epochs star was found in --> should equal num_epochs


	################################################################################################
	# 						GET COORDINATES FROM FITS HEADER  
	################################################################################################

	# The coordinates in the dataframe that need to be compared to the header are from the 'master'
	# input file which is always 3p6 e01 d6 

	master_input = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch1/e01/' + star_name + '_3p6um_e01_d6_cbcd_dn.fits'

	# Convert RA and Dec to pixel coordinates using astropy.wcs
	hdulist = fits.open(master_input)		
	w = wcs.WCS(hdulist[0].header) # Parse the WCS keywords in the primary HDU
	world = np.array([[ra, dec]], np.float_) # the world coordinates of the target Cepheid
	pix = w.wcs_world2pix(world,1) # convert world coordinates to pixel coordinates in an array. For just list of coords want pix[0]

	# Coordinates we want to compare dataframe with with
	true_x = pix[0][0] 
	true_y = pix[0][1]

	print "NED coordinates X = " + str(true_x) + " Y = " + str(true_y)

	################################################################################################
	# 				COMPARE TRUE COORDS WITH DF COORDS TO FIND TARGET V*  
	################################################################################################

	# Want to search for stars within an increasing search box until star found

	star_found = False

	# Initial increment of 2 from NED coordinates
	increment = 2

	while star_found == False:

		increment += 1 # increase search radius

		# Set up of initial search box size
		x_low = true_x - increment
		x_high = true_x + increment
		y_low = true_y - increment
		y_high = true_y + increment	

		# Filter df to only keep stars within this box
		box_df = df[(df['X']>x_low) & (df['X']<x_high) & (df['Y']>y_low) & (df['Y']<y_high)]
		box_df.reset_index(drop=True, inplace=True)

		# Print number of stars in this filter df
		no_potentials = len(box_df)

		# If only one star found, assume this is the target cepheid and write their data to file
		if no_potentials == 1:

			print "Only one star found. Writing this data to file."
			print "Star found at X = %f  Y = %f" % (box_df['X'], box_df['Y'])

			star_found = True # stop the loop

			# Write to file
			for epoch in range(1, num_epochs+1):

			    if epoch < 10:
			        epoch_string = 'e0' + str(epoch)
			    else: epoch_string = 'e' + str(epoch)

			    f.write("%d %d %.2f %.2f %.4f %.4f %.4f %.4f \n" % (epoch, box_df['ID'], box_df['X'], box_df['Y'], box_df['3p6_' + epoch_string + '_ave_mag'], box_df['3p6_' + epoch_string + '_ave_err'], box_df['4p5_' + epoch_string + '_ave_mag'], box_df['4p5_' + epoch_string + '_ave_err']))			


		elif no_potentials == 0:

			print "No stars found within the search box."

			# fail_file = '/home/ac833/Magnitudes/failed_stars.txt'

			# g = open(fail_file, 'a')
			# g.write(star_name + ': No stars found within the search box.')
			# g.close()

		else:

			print "%d stars found in this search box. Now computing variability index." % no_potentials

			# Do variability index and pick star with largest variability index

			################################################################################################
			# 						1. COMPUTE WEIGHTED MEAN MAGNITUDES  
			################################################################################################

			# Add empty columns for the weighted mean magnitudes
			df['3p6_weighted_mean_mag'] = np.nan
			df['4p5_weighted_mean_mag'] = np.nan

			# Now want to calculate the weighted mean magnitude per star at each wavelength across all epochs

			for channel in [1,2]:
			    
			    if channel == 1:
			        wavelength = '3p6'
			        start_dither = 6
			        F0 = 280.9
			    else: 
			        wavelength = '4p5'
			        start_dither = 1
			        F0 = 179.7
			        
			    for index, row in box_df.iterrows():

			        numerator = 0
			        denominator = 0

			        for epoch in range(1,num_epochs+1):

			            if epoch < 10:
			                epoch_string = 'e0' + str(epoch)
			            else: epoch_string = 'e' + str(epoch) 

			            mag_col = wavelength + '_' + epoch_string + '_ave_mag'
			            err_col = wavelength + '_' + epoch_string + '_ave_err'
			           
		                numerator += box_df[mag_col][index]/(box_df[err_col][index]**2)
		                denominator += 1/(box_df[err_col][index]**2)
			        
			        # Compute weighted mean mag as in Welch & Stetson (1993)
		            weighted_mean_mag = numerator/denominator
		            column = wavelength + '_weighted_mean_mag'
		            box_df.ix[index, column] = weighted_mean_mag


			################################################################################################
			# 						2. COMPUTE NORMALISED MAGNITUDE RESIDUALS  
			################################################################################################

			# Add extra columns in df for normalised magnitude residuals
			for epoch in range(1,num_epochs+1):

			    if epoch < 10:
			        epoch_string = 'e0' + str(epoch)
			    else: epoch_string = 'e' + str(epoch)
			    
			    box_df['3p6_' + epoch_string + '_res'] = np.nan
			    box_df['4p5_' + epoch_string + '_res'] = np.nan

			for channel in [1,2]:
	    
			    if channel == 1:
			        wavelength = '3p6'

			    else: 
			        wavelength = '4p5'
			        
			    for index, row in box_df.iterrows():

			        for epoch in range(1,num_epochs+1):

			            if epoch < 10:
			                epoch_string = 'e0' + str(epoch)
			            else: epoch_string = 'e' + str(epoch)
			            
			            res_col = wavelength + '_' + epoch_string + '_res'
			            mag_col = wavelength + '_' + epoch_string + '_ave_mag'
			            weight_col = wavelength + '_weighted_mean_mag'
			            err_col = wavelength + '_' + epoch_string + '_ave_err'
			            
			            epoch_mag = float(box_df[mag_col][index])
			            weighted_mag = float(box_df[weight_col][index])
			            error = float(box_df[err_col][index])
			            
			            box_df.ix[index, res_col] = (epoch_mag - weighted_mag)/error 


			################################################################################################
			# 						3. COMPUTE VARIABILITY INDEX, I
			################################################################################################

			# Add new I column
			box_df['I'] = np.nan

			# This factor is the first term in the calculation of I from Welch & Stetson
			epoch_factor = sqrt(1/(float(num_epochs)*(float(num_epochs-1))))

			for index, row in box_df.iterrows():
			    
			    summation = 0

			    for epoch in range(1,num_epochs+1):

			        if epoch < 10:
			            epoch_string = 'e0' + str(epoch)
			        else: epoch_string = 'e' + str(epoch)
			        
			        col_3p6 = '3p6_' + epoch_string + '_res'
			        col_4p5 = '4p5_' + epoch_string + '_res'
			        
			        res_3p6 = box_df[col_3p6][index]
			        res_4p5 = box_df[col_4p5][index]
			        
			        summation += res_3p6 * res_4p5
			    
			    box_df.ix[index, 'I'] = epoch_factor * summation

			# Index with largest value of I
			star_index = box_df['I'].idxmax()

			print box_df 

			# Write to file the star at this index
			for epoch in range(1, num_epochs+1):

			    if epoch < 10:
			        epoch_string = 'e0' + str(epoch)
			    else: epoch_string = 'e' + str(epoch)

			    f.write("%d %d %.2f %.2f %.4f %.4f %.4f %.4f \n" % (int(epoch), int(box_df['ID'][star_index]), box_df['X'][star_index], box_df['Y'][star_index], box_df['3p6_' + epoch_string + '_ave_mag'][star_index], box_df['3p6_' + epoch_string + '_ave_err'][star_index], box_df['4p5_' + epoch_string + '_ave_mag'][star_index], box_df['4p5_' + epoch_string + '_ave_err'][star_index]))			

			star_found = True

	f.close()

	return(0)

# Format GLOESS file 
def gloess_files(star_name, galaxy, num_epochs, period):

	################################################################################################
	# 								GET HMJD FOR EACH EPOCH
	################################################################################################

	hmjd_3p6 = []
	hmjd_4p5 = []

	for j in range(1,num_epochs+1):

		if j < 10:
			filename = '/home/ac833/Data/'+str(galaxy)+'/BCD/'+str(star_name)+'/ch1/e0'+str(j)+'/'+str(star_name)+'_3p6um_e0'+str(j)+"_d1_cbcd_dn.fits"
		else: filename = '/home/ac833/Data/'+str(galaxy)+'/BCD/'+str(star_name)+'/ch1/e'+str(j)+'/'+str(star_name)+'_3p6um_e'+str(j)+"_d1_cbcd_dn.fits"

		image = fits.open(filename)
		header = image[0].header
		hmjd_3p6.append(header['HMJD_OBS'])

	for j in range(1,num_epochs+1):

		if j < 10:
			filename = '/home/ac833/Data/'+str(galaxy)+'/BCD/'+str(star_name)+'/ch2/e0'+str(j)+'/'+str(star_name)+'_4p5um_e0'+str(j)+"_d1_cbcd_dn.fits"
		else: filename = '/home/ac833/Data/'+str(galaxy)+'/BCD/'+str(star_name)+'/ch2/e'+str(j)+'/'+str(star_name)+'_4p5um_e'+str(j)+"_d1_cbcd_dn.fits"

		image = fits.open(filename)
		header = image[0].header
		hmjd_4p5.append(header['HMJD_OBS'])


	################################################################################################
	# 								OPEN MAGNITUDES FILE INTO DF
	################################################################################################

	df = pd.read_csv('/home/ac833/Magnitudes/' + galaxy + '/' + star_name + '.txt', header=0, delim_whitespace=True)


	################################################################################################
	# 								WRITE OUT TO FILE
	################################################################################################

	f = open('/home/ac833/GLOESS_files/' + star_name + '_gloess.txt', 'w')

	f.write(str(star_name) + '\n')
	f.write(str(period) + '\n')
	f.write(str(1) + '\n') # doesn't matter what goes here as it's not used by GLOESS code

	# If short period then slightly relaxed fitting parameter of 0.25
	# Longer period Cepheids generally have better light curves so can use a tighter fitting parameter of 0.20
	if period <= 12:
		f.write(str(0.25) + ' ' + str(0.25) + '\n') 
	else: f.write(str(0.2) + ' ' + str(0.2) + '\n')

	# Write out the 3p6 data first
	for epoch in range(0, num_epochs):

		f.write("%f %f %f 99.9999 9.9999 \n" % (hmjd_3p6[epoch], df['3p6_mag'][epoch], df['3p6_err'][epoch]))

	# Then write out 4p5 data
	for epoch in range(0, num_epochs):

		f.write("%f 99.9999 9.9999 %f %f \n" % (hmjd_4p5[epoch], df['4p5_mag'][epoch], df['4p5_err'][epoch]))

	f.close()

# Fit multiband GLOESS curve
def gloess_multiband(star_name, galaxy):

	################################################################################################
	# 									READ IN DATA
	################################################################################################

	# This file has all the magnitudes
	input_filename = '/home/ac833/GLOESS_files/' + star_name + '_gloess.txt'

	dir1 = [] # 3p6 magnitudes
	dir2 = [] # 4p5 magnitudes
	deir1 = [] # 3p6 errors
	deir2 = [] # 4p5 errors
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
			xir1 = float(data[0])
			xir2 = float(data[1])
		if counter > 3:
			# if counter>3, then this line of the input consists of many values that get appended to the lists at the beginning of the code
			# I think the d in front of du say is just so that when convert to numpy arrays the u isn't already taken
			dmjd.append(float(data[0])) # hmjd
			dir1.append(float(data[1])) 
			deir1.append(float(data[2])) 
			dir2.append(float(data[3]))
			deir2.append(float(data[4]))

		# Counter incrememnts so that each line from the input file is read 
		counter  = counter + 1	
			
	## Read in all the data from the file and filled the arrays. Need to convert these to numpy arrays.

	# Number of lines of data as opposed to the Name, period etc
	number = counter - 4 # Number data lines in the file

	# Convert to numpy arrays
	ir1 = np.array(dir1)
	ir2 = np.array(dir2)
	eir1 = np.array(deir1)
	eir2 = np.array(deir2)
	mjd = np.array(dmjd)	

	# Number of observations - doesn't count 99.9999 values
	nir1 = sum(ir1<50)
	nir2 = sum(ir2<50)

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

	if nir2 > 0:
		maxvals.append(np.amax(ir2[ir2<50]))
		minvals.append(np.amin(ir2[ir2<50]))

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
	plt.figure(figsize=[5,7])

	# GridSpec specifies the geometry of the grid that a subplot will be placed. 
	# The number of rows and number of columns of the grid need to be set.
	# So here, we have 3 rows and 4 columns
	gs = gridspec.GridSpec(3, 4)
	ax2 = plt.subplot(gs[0, :])
	ax3 = plt.subplot(gs[1, :])
	ax4 = plt.subplot(gs[2, :])

	# Values in the [] are xmin,xmax,ymin,ymax
	#ax1.axis([1,3.5,(maxlim),(minlim)])

	titlestring = cepname + ', P = ' + str(period) + ' days'
	plt.suptitle(titlestring, fontsize=20)

	#ax1.set_ylabel('Magnitude')
	#ax1.set_xlabel('Phase $\phi$')

	## Fitting and plotting for each band
	if nir1 > 0:

		ir11, ir1x, yir1, yeir1, xphaseir1 = gf.fit_one_band(ir1,eir1,phase,nir1,xir1)

		# ax1.plot(ir1x,ir11,'k-')
	 # 	ax1.plot(xphaseir1,yir1,color='MediumVioletRed',marker='o',ls='None', label='mag')

		aveir1, adevir1, sdevir1, varir1, skewir1, kurtosisir1, ampir1 = gf.moment(ir11[200:300],100)

		if phased == 1:
			factor = sqrt(nir1)
		if phased == 0:
			factor = 1 
		if nir1 > 1:
			print  '<mag> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f}' .format(aveir1, sdevir1/factor, ampir1)
		if nir1 == 1:
			print  'mag = {0:.3f} --- single point'.format(aveir1)


	if nir2 > 0:

		ir21, ir2x, yir2, yeir2, xphaseir2 = gf.fit_one_band(ir2,eir2,phase,nir2,xir2)

		# ax1.plot(ir2x,ir21,'k-')
	 # 	ax1.plot(xphaseir2,yir2,color='DeepPink',marker='o',ls='None', label='mag')

		aveir2, adevir2, sdevir2, varir2, skewir2, kurtosisir2, ampir2 = gf.moment(ir21[200:300],100)

		if phased == 1:
			factor = sqrt(nir2)
		if phased == 0:
			factor = 1 
		if nir2 > 1:
			print  '<mag> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f}' .format(aveir2, sdevir2/factor, ampir2)
		if nir2 == 1:
			print  'mag = {0:.3f} --- single point'.format(aveir1)

	#handles, labels = ax1.get_legend_handles_labels() 
	#ax1.legend(handles[::-1],labels[::-1],loc=4, numpoints=1,prop={'size':10})

	### Define the colour curve
	colour_curve = ir11 - ir21
	## Define the colour points
	ch1_points = yir1[yir1<99]
	ch2_points = yir2[yir2<99]

	colour_points = ch1_points - ch2_points
	colour_phases = xphaseir1[yir1<99]

	colour_points = np.concatenate((colour_points,colour_points,colour_points,colour_points,colour_points))
	colour_phases = np.concatenate((colour_phases,(colour_phases+1.),(colour_phases+2.),(colour_phases+3.),(colour_phases+4.)))


	avecol, adevcol, sdevcol, varcol, skewcol, kurtosiscol, ampcol = gf.moment(colour_curve[200:300],100)

	#print >> avsout, '<[3.6] - [4.5]> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f}' .format(avecol, sdevcol/factor, ampcol)
	print  '<[3.6] - [4.5]> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f}' .format(avecol, sdevcol/factor, ampcol)

	ax2.axis([1,3.5,(np.average(ir11[200:300]) + 0.4),(np.average(ir11[200:300]) - 0.4)])
	ax2.yaxis.tick_right()
	ax2.plot(ir1x,ir11,'k-')
	ax2.plot(xphaseir1,yir1,color='MediumVioletRed',marker='o',ls='None', alpha=0.6, label='$[3.6]$')
	ax2.annotate('$[3.6]$', xy=(0.04, 0.8375), xycoords='axes fraction', fontsize=16)

	ax3.axis([1,3.5,(np.average(ir21[200:300]) + 0.4),(np.average(ir21[200:300]) - 0.4)])
	ax3.yaxis.tick_right()
	ax3.plot(ir2x,ir21,'k-')
	ax3.plot(xphaseir2,yir2,color='DeepPink',marker='o',ls='None', alpha=0.6, label='$[3.6]$')
	ax3.annotate('$[4.5]$', xy=(0.04, 0.8375), xycoords='axes fraction',fontsize=16)

	myaxis2 = [1,3.5, 0.4,-0.4]
	ax4.axis(myaxis2)
	ax4.yaxis.tick_right()
	ax4.yaxis.set_major_locator(plt.FixedLocator([0.2, 0.1, 0, -0.1, -0.2]))
	ax4.plot(ir1x,colour_curve,'k-')
	ax4.plot(colour_phases,colour_points,color='Black',marker='o',ls='None', alpha=0.3, label='$[3.6]-[4.5]$')

	ax4.set_xlabel('Phase $\phi$')
	ax4.annotate('$[3.6] - [4.5]$', xy=(0.04, 0.8375), xycoords='axes fraction',fontsize=16)

	ax4.hlines(0,1,3.5,'k','dashdot')

	plt.setp(ax2.get_xticklabels(),visible=False)
	plt.setp(ax3.get_xticklabels(),visible=False)

	plotname = cepname+'.eps'
	plt.savefig(plotname, transparent='True')

	#avsout.close()

	# Save light curve
	lc_name = '/home/ac833/Light_Curves/' + galaxy +'/' + star_name + '.png'
	plt.savefig(lc_name, dpi=100)
	#plt.savefig(lc_name, transparent=True)

	# Output all relevant information to a file
	# Galaxy, star, channel, period, mean mag, std_dev, amplitude

	master_mags_file = '/home/ac833/master_stars.txt'

	f = open(master_mags_file, 'a')

	if phased == 1:
		std_error_ir1 = sdevir1/(factor ** 2) # error scales with 1/N
		std_error_ir2 = sdevir2/(factor ** 2)
	else:
		std_error_ir1 = sdevir1/factor # error scales with 1/sqrt(N)
		std_error_ir2 = sdevir2/factor # error scales with 1/sqrt(N)


	f.write("%s %s %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f \n" % (galaxy, star_name, period, aveir1, std_error_ir1, ampir1, aveir2, std_error_ir2, ampir2, avecol, sdevcol/factor, ampcol))

	f.close()

# Master image off target
def master_off_target_no_psf(star_name, galaxy, channel, wavelength, epoch_number, num_epochs):

	num_images = 5 # for off-target field, medianed images are made for each epoch rather than a combined one across all epochs

	if channel == '1':
		field = '1' # the off-target images are the first 5 dithers d1-5
		start_dither = 1
	else:
		field = '2' # the off-target images are the second 5 dithers d6-10
		start_dither = 6

	shutil.copy(star_name + '_on_target_master.psf', star_name + '_off_target_master.psf')

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

	# Make initial file manually and then let DAOMASTER refine the transformations later

	mch_filename = star_name + '_off_target.mch'

	# Open file
	f = open(mch_filename, 'w')

	file1 = star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither) + '_cbcd_dn.ap'
	file2 = star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+1) + '_cbcd_dn.ap'
	file3 = star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+2) + '_cbcd_dn.ap'
	file4 = star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+3) + '_cbcd_dn.ap'
	file5 = star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+4) + '_cbcd_dn.ap'

	# Write transfomration lines to file
	f.write("'%s' 0.0000 0.0000 1.0000 0.0000 0.0000 1.0000 0.0000 0.0000 \n" % file1)
	f.write("'%s' -10.0000 33.0000 1.0000 0.0000 0.0000 1.0000 0.0000 0.0000 \n" % file2)
	f.write("'%s' -42.0000 -18.0000 1.0000 0.0000 0.0000 1.0000 0.0000 0.0000 \n" % file3)	
	f.write("'%s' 4.0000 -9.5000 1.0000 0.0000 0.0000 1.0000 0.0000 0.0000 \n" % file4)
	f.write("'%s' 47.0000 -26.0000 1.0000 0.0000 0.0000 1.0000 0.0000 0.0000 \n" % file5)

	f.close()
	
	# daomatch = pexpect.spawn('daomatch')

	# # Set up log file
	# fout = file('daomatch_log.txt','w')
	# daomatch.logfile = fout

	# # Give DAOMATCH all the 5 dithers to be used in making the medianed image
	# # The '/' after the file names lets DAOMATCH know that the scale is the same and it is only the rotation and shifts that might have changed
	# daomatch.expect("Master input file:")
	# daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither) + '_cbcd_dn.ap') # Give it the first BCD phot file
	# daomatch.expect("Output file name")
	# daomatch.sendline(star_name + '_off_target.mch')
	# daomatch.expect("Next input file:")
	# daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+1) + '_cbcd_dn.ap/') # Give it the second BCD phot file
	# daomatch.expect("Next input file")
	# daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+2) + '_cbcd_dn.ap/') # Give it the third BCD phot file
	# daomatch.expect("Next input file")
	# daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+3) + '_cbcd_dn.ap/') # Give it the fourth BCD phot file
	# daomatch.expect("Next input file")
	# daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+4) + '_cbcd_dn.ap/') # Give it the fifth BCD phot file
	# daomatch.expect("Next input file")
	# daomatch.sendline("") # exit

	# Refine transformations with DAOMASTER
	daomaster = pexpect.spawn('daomaster')

	# Set up log file
	fout = file('daomaster_log.txt','w')
	daomaster.logfile = fout

	daomaster.expect("File with list of input files:")
	daomaster.sendline(star_name + '_off_target.mch')
	daomaster.expect("Minimum number, minimum fraction, enough frames:")
	daomaster.sendline("1, 0.5, 5") 
	daomaster.expect("Maximum sigma:")
	daomaster.sendline("99") 
	daomaster.expect("Your choice:")
	daomaster.sendline("6") # solve for 6 degrees of freedom
	daomaster.expect("Critical match-up radius:")
	daomaster.sendline("7") 
	
	for dither in range(start_dither+1,start_dither+5):
		#daomaster.expect(star_name + '_' + wavelength + '_e' + epoch_number + '_d' + str(dither) + '_cbcd_dn.ap')
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
	daomaster.sendline(star_name + '_off_target.mch')
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
	montage2.sendline(star_name + '_off_target.mch')
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
	montage2.sendline(star_name + '_off_target.fits')


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
	daophot.sendline("at " + star_name + '_off_target.fits')
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
	daophot.sendline(star_name + '_off_target.coo')
	daophot.expect("Output file")
	daophot.sendline(star_name + '_off_target.ap')

	# Run ALLSTAR
	allstar = pexpect.spawn('allstar')

	fout = file('allstar_log.txt', 'w')
	allstar.logfile = fout

	allstar.expect("OPT>")
	allstar.sendline("")
	allstar.expect("Input image name:")
	allstar.sendline(star_name + '_off_target.fits')
	allstar.expect("File with the PSF")
	allstar.sendline(star_name + '_off_target_master.psf')
	allstar.expect("Input file")
	allstar.sendline(star_name + '_off_target.ap')
	allstar.expect("File for results")
	allstar.sendline(star_name + '_off_target.als')
	allstar.expect("Name for subtracted image")
	allstar.sendline(star_name + '_off_target_dns.fits')
	allstar.expect("Good bye")
	allstar.close(force=True)	

	################################################################################################
	# 							ADD OFFSETS BACK IN TO STAR LIST
	################################################################################################	

	daophot = pexpect.spawn('daophot')

	daophot.expect("Command:")
	daophot.sendline("off") # offsets to put x and y back in
	daophot.expect("Input file name:")
	daophot.sendline(star_name + '_off_target.als')
	daophot.expect("Additive offsets ID, DX, DY, DMAG:")
	daophot.sendline("0," + offsets[0] + "," + offsets[1] + ",0")
	daophot.expect("Output file name")
	daophot.sendline(star_name + '_off_target.mag') # this is the final file for star positions that we want 
	daophot.expect("Command:")
	daophot.sendline("ex")
	daophot.close(force=True)

	################################################################################################
	# 							TIDYING UP AND RENAMING FILES 
	################################################################################################	

	# Change .ap file names in .mch file to .als ready for next script
	f = open(star_name +'_off_target.mch', 'r')
	filedata = f.read()
	f.close()

	newdata = filedata.replace(".ap",".als")

	f = open(star_name +'_off_target_master.mch','w')
	f.write(newdata)
	f.close()	

	# Rename files we want for the next scripts to be consistent with the on-target version of this code
	os.rename(star_name + '_off_target.fits', star_name + '_off_target_master.fits')
	os.rename(star_name + '_off_target.mag', star_name + '_off_target_master.mag')

	return(0)

# Master image off target
def master_off_target_update(star_name, galaxy, channel, wavelength, epoch_number, num_epochs):

	num_images = 5 # for off-target field, medianed images are made for each epoch rather than a combined one across all epochs

	if channel == '1':
		field = '1' # the off-target images are the first 5 dithers d1-5
		start_dither = 1
	else:
		field = '2' # the off-target images are the second 5 dithers d6-10
		start_dither = 6

	# Copy on-target PSF to off-target PSF 
	if channel == '1':
		on_field = '2'
	else: on_field = '1' 

	shutil.copy(star_name + '_' + wavelength + '_f' + on_field + '_master.psf', star_name + '_' + wavelength + '_f' + field + '_master.psf')

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
	
	daomatch = pexpect.spawn('daomatch', timeout=10)

	# Set up log file
	fout = file('daomatch_log.txt','w')
	daomatch.logfile = fout

	# Give DAOMATCH all the 5 dithers to be used in making the medianed image
	# The '/' after the file names lets DAOMATCH know that the scale is the same and it is only the rotation and shifts that might have changed
	daomatch.expect("Master input file:")
	daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither) + '_cbcd_dn.ap') # Give it the first BCD phot file
	daomatch.expect("Output file name")
	daomatch.sendline(star_name + '_' + wavelength +'_f' + field + '.mch')

	# contents = ""
	# with open('daomatch_log.txt') as f:
	# 	for line in f.readlines():
	# 		contents += line

	# if "Write this transformation?" in contents:
	# 	print "Fail here"

	# f = open(star_name + '_on_target.mch','w')
	# f.write(newdata)
	# f.close()

	print star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+1) + '_cbcd_dn.ap/'

	try: 
		print "Attempting to match file"
		daomatch.expect("Next input file:")
		daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+1) + '_cbcd_dn.ap/') # Give it the second BCD phot file

	except pexpect.TIMEOUT :

		print "Can't match. Session timed out."

	# contents = ""
	# with open('daomatch_log.txt') as f:
	# 	for line in f.readlines():
	# 		contents += line

	# if "Write this transformation?" in contents:
	# 	print "Fail here"

	print star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+2) + '_cbcd_dn.ap/'

	daomatch.expect("Next input file")
	daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+2) + '_cbcd_dn.ap/') # Give it the third BCD phot file

	try: 
		print "Attempting to match file"
		daomatch.expect("Next input file:")
		daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+2) + '_cbcd_dn.ap/') # Give it the second BCD phot file

	except pexpect.TIMEOUT :

		print "Can't match. Session timed out."


	# contents = ""
	# with open('daomatch_log.txt') as f:
	# 	for line in f.readlines():
	# 		contents += line

	# if "Write this transformation?" in contents:
	# 	print "Fail here"

	print star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+3) + '_cbcd_dn.ap/'

	daomatch.expect("Next input file")
	daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+3) + '_cbcd_dn.ap/') # Give it the fourth BCD phot file

	try: 
		print "Attempting to match file"
		daomatch.expect("Next input file:")
		daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+3) + '_cbcd_dn.ap/') # Give it the second BCD phot file

	except pexpect.TIMEOUT :

		print "Can't match. Session timed out."

	# contents = ""
	# with open('daomatch_log.txt') as f:
	# 	for line in f.readlines():
	# 		contents += line

	# if "Write this transformation?" in contents:
	# 	print "Fail here"

	print star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+4) + '_cbcd_dn.ap/'

	daomatch.expect("Next input file")
	daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+4) + '_cbcd_dn.ap/') # Give it the fifth BCD phot file

	try: 
		print "Attempting to match file"
		daomatch.expect("Next input file:")
		daomatch.sendline(star_name + '_' + wavelength + '_e' + epoch_number +'_d' + str(start_dither+4) + '_cbcd_dn.ap/') # Give it the second BCD phot file

	except pexpect.TIMEOUT :

		print "Can't match. Session timed out."


	# contents = ""
	# with open('daomatch_log.txt') as f:
	# 	for line in f.readlines():
	# 		contents += line

	# print contents

	# if "Write this transformation?" in contents:
	# 	print "Fail here"

	daomatch.expect("Next input file")
	daomatch.sendline("") # exit


	# # Refine transformations with DAOMASTER
	# daomaster = pexpect.spawn('daomaster')

	# # Set up log file
	# fout = file('daomaster_log.txt','w')
	# daomaster.logfile = fout

	# daomaster.expect("File with list of input files:")
	# daomaster.sendline(star_name + '_' + wavelength + '_f' + field + '.mch')
	# daomaster.expect("Minimum number, minimum fraction, enough frames:")
	# daomaster.sendline("1, 0.5, 5") 
	# daomaster.expect("Maximum sigma:")
	# daomaster.sendline("99") 
	# daomaster.expect("Your choice:")
	# daomaster.sendline("6") # solve for 6 degrees of freedom
	# daomaster.expect("Critical match-up radius:")
	# daomaster.sendline("7") 
	
	# for dither in range(start_dither+1,start_dither+5):
	# 	daomaster.expect(star_name + '_' + wavelength + '_e' + epoch_number + '_d' + str(dither) + '_cbcd_dn.ap')
	# 	daomaster.sendline("")

	# # Repeat with decreasing match up size
	# for match_up in range(7,-1,-1):
	# 	daomaster.expect("New match-up radius")
	# 	daomaster.sendline(str(match_up))

	# # Options for different output files; we only require the refined transformations
	# daomaster.expect("Assign new star IDs?")
	# daomaster.sendline("y") # assign new ids so all frames have same ids
	# daomaster.expect("A file with mean magnitudes and scatter?")
	# daomaster.sendline("n")
	# daomaster.expect("A file with corrected magnitudes and errors?")
	# daomaster.sendline("n")
	# daomaster.expect("A file with raw magnitudes and errors?")
	# daomaster.sendline("n")
	# daomaster.expect("A file with the new transformations?")
	# daomaster.sendline("y")
	# daomaster.expect("Output file name")
	# daomaster.sendline(star_name + '_' + wavelength +'_f' + field + '.mch')
	# daomaster.expect("New output file name")
	# daomaster.sendline("")
	# daomaster.expect("A file with the transfer table?")
	# daomaster.sendline("e") # exits rest of options

	# ################################################################################################
	# # 								MAKE MEDIANED IMAGE
	# ################################################################################################

	# montage2 = pexpect.spawn('montage2')

	# # Set up log file - need this to obtain the offsets in X and Y
	# fout = file('montage_log.txt','w')
	# montage2.logfile = fout

	# montage2.expect("File with transformations:")
	# montage2.sendline(star_name + '_' + wavelength +'_f' + field + '.mch')
	# montage2.expect("Image-name suffix:")
	# montage2.sendline("")
	# montage2.expect("Minimum number of frames, percentile:")
	# montage2.sendline("1,0.5") # play around with minimum number of frames
	# montage2.expect("X limits of output image:")
	# montage2.sendline("e")
	# montage2.expect("Y limits of output image:")
	# montage2.sendline("e")
	# montage2.expect("Expansion factor:")
	# montage2.sendline("1") # creates image with same scale as bcd images
	# montage2.expect("Determine sky from overlap region?")
	# montage2.sendline("y")
	# montage2.expect("Name for output image")
	# montage2.sendline(star_name + '_' + wavelength +'_f' + field + '.fits')


	# ################################################################################################
	# # 									OBTAIN OFFSETS
	# ################################################################################################

	# log = open('montage_log.txt', 'r')
	# lines = log.readlines()

	# offsets = []
	
	# for line in lines:
	# 	if "Offsets" in line:

	# 		offsets.append(line.split(' ')[-3])
	# 		offsets.append(line.split(' ')[-2])		


	# ################################################################################################
	# # 								CREATE MASTER STAR LIST
	# ################################################################################################

	# daophot = pexpect.spawn('daophot') 

	# # Set up logfile
	# fout = file('daophot_log.txt','w')
	# daophot.logfile = fout

	# # Attach medianed image
	# daophot.expect("Command:")
	# daophot.sendline("at " + star_name + '_' + wavelength + '_f' + field + '.fits')
	# daophot.expect("Command:")
	# daophot.sendline("opt")
	# daophot.expect("File with parameters")
	# daophot.sendline("")
	# daophot.expect("OPT>")
	# daophot.sendline("th=20") # set to a higher threshold to account for higher S/N 
	# daophot.expect("OPT>")
	# daophot.sendline("")

	# daophot.expect("Command:")
	# daophot.sendline("fi")
	# daophot.expect("Number of frames averaged, summed:")
	# daophot.sendline("5,1") 
	# daophot.expect("File for positions")
	# daophot.sendline("")
	# daophot.expect("Are you happy with this?")
	# daophot.sendline("y")

	# daophot.expect("Command:")
	# daophot.sendline("ph")
	# daophot.expect("File with aperture radii")
	# daophot.sendline("")
	# daophot.expect("PHO>")
	# daophot.sendline("")
	# daophot.expect("Input position file")
	# daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.coo')
	# daophot.expect("Output file")
	# daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.ap')

	# # Run ALLSTAR
	# allstar = pexpect.spawn('allstar')

	# fout = file('allstar_log.txt', 'w')
	# allstar.logfile = fout

	# allstar.expect("OPT>")
	# allstar.sendline("")
	# allstar.expect("Input image name:")
	# allstar.sendline(star_name + '_' + wavelength + '_f' + field + '.fits')
	# allstar.expect("File with the PSF")
	# allstar.sendline(star_name + '_' + wavelength + '_f' + field + '_master.psf')
	# allstar.expect("Input file")
	# allstar.sendline(star_name + '_' + wavelength + '_f' + field + '.ap')
	# allstar.expect("File for results")
	# allstar.sendline(star_name + '_' + wavelength + '_f' + field + '.als')
	# allstar.expect("Name for subtracted image")
	# allstar.sendline(star_name + '_' + wavelength + '_f' + field + '_dns.fits')
	# allstar.expect("Good bye")
	# allstar.close(force=True)	

	# ################################################################################################
	# # 							ADD OFFSETS BACK IN TO STAR LIST
	# ################################################################################################	

	# daophot = pexpect.spawn('daophot')

	# daophot.expect("Command:")
	# daophot.sendline("off") # offsets to put x and y back in
	# daophot.expect("Input file name:")
	# daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.als')
	# daophot.expect("Additive offsets ID, DX, DY, DMAG:")
	# daophot.sendline("0," + offsets[0] + "," + offsets[1] + ",0")
	# daophot.expect("Output file name")
	# daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.mag') # this is the final file for star positions that we want 
	# daophot.expect("Command:")
	# daophot.sendline("ex")
	# daophot.close(force=True)

	# ################################################################################################
	# # 							TIDYING UP AND RENAMING FILES 
	# ################################################################################################	

	# # Change .ap file names in .mch file to .als ready for next script
	# f = open(star_name +'_' + wavelength + '_f' + field + '.mch', 'r')
	# filedata = f.read()
	# f.close()

	# newdata = filedata.replace(".ap",".als")

	# f = open(star_name +'_' + wavelength + '_f' + field + '.mch','w')
	# f.write(newdata)
	# f.close()	

	# # Rename files we want for the next scripts to be consistent with the on-target version of this code
	# os.rename(star_name + '_' + wavelength + '_f' + field + '.fits', star_name + '_' + wavelength + '_f' + field + '_master.fits')
	# os.rename(star_name + '_' + wavelength + '_f' + field + '.mch', star_name + '_' + wavelength + '_f' + field + '_master.mch')
	# os.rename(star_name + '_' + wavelength + '_f' + field + '.mag', star_name + '_' + wavelength + '_f' + field + '_master.mag')



	return(0)

# Master image on target
def master_on_target_both_channels(star_name, galaxy, num_epochs):

	# if channel == '1':
	# 	field = '2' # the on-target images are the second 5 dithers d6-10
	# 	start_dither = 6
	# else:
	# 	field = '1' # the on-target images are the first 5 dithers d1-5
	# 	start_dither = 1

	num_images = num_epochs * 5 * 2 # 5 dithers per epoch for both channels

	# List of files to use
	files = []

	# All 3p6 on target files
	for epoch in range(1,num_epochs+1):
		for dither in range(6,11): 

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			dither = str(dither)

			filename = star_name + '_3p6um_e' + epoch + '_d' + dither + '_cbcd_dn.ap'
			files.append(filename)	

	# All 4p5 on target files
	for epoch in range(1,num_epochs+1):
		for dither in range(1,6): 

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			dither = str(dither)

			filename = star_name + '_4p5um_e' + epoch + '_d' + dither + '_cbcd_dn.ap'
			files.append(filename)	

	# Create temporary folder to make image - will then be copied to each epoch folder
	temp = '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/temp/'

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

		cwd = '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch1/e'+epoch+'/'

		# Copy 3p6 files
		for dither in range(6,11):

			dither = str(dither)

			shutil.copyfile(cwd + star_name + '_3p6um_e' + epoch + '_d' + dither + '_cbcd_dn.ap', temp + star_name + '_3p6um_e' + epoch + '_d' + dither + '_cbcd_dn.ap') #.ap
			shutil.copyfile(cwd + star_name + '_3p6um_e' + epoch + '_d' + dither + '_cbcd_dn.fits', temp + star_name + '_3p6um_e' + epoch + '_d' + dither + '_cbcd_dn.fits')

	# Copy all FITS images and aperture photometry files to temp folder
	for epoch in range(1,num_epochs+1):

		# Get epoch in correct format
		if epoch < 10:
			epoch = '0' + str(epoch)
		else: epoch = str(epoch)

		cwd = '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch2/e'+epoch+'/'

		# Copy 3p6 files
		for dither in range(1,6):

			dither = str(dither)

			shutil.copyfile(cwd + star_name + '_4p5um_e' + epoch + '_d' + dither + '_cbcd_dn.ap', temp + star_name + '_4p5um_e' + epoch + '_d' + dither + '_cbcd_dn.ap') #.ap
			shutil.copyfile(cwd + star_name + '_4p5um_e' + epoch + '_d' + dither + '_cbcd_dn.fits', temp + star_name + '_4p5um_e' + epoch + '_d' + dither + '_cbcd_dn.fits')

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
	daomatch.sendline(files[0]) 
	daomatch.expect("Output file name")
	daomatch.sendline(star_name + '_on_target.mch')	

	# Give it the rest of the images
	for j in range(1,len(files)):
	#for j in range(1,num_images):

		daomatch.expect("Next input file")
		daomatch.sendline(files[j]+'/')

	daomatch.expect("Next input file")
	daomatch.sendline("") # exit

	# Use DAOMASTER to refine the transformations
	daomaster = pexpect.spawn('daomaster')

	fout = file('daomaster_log.txt','w')
	daomaster.logfile = fout

	daomaster.expect("File with list of input files:")
	daomaster.sendline(star_name + '_on_target.mch')
	daomaster.expect("Minimum number, minimum fraction, enough frames:")
	daomaster.sendline("1, 0.5, " + str(num_images)) # play around with these values
	daomaster.expect("Maximum sigma:")
	daomaster.sendline("99") 
	daomaster.expect("Your choice:")
	daomaster.sendline("6") # solve for 6 degrees of freedom
	daomaster.expect("Critical match-up radius:")
	daomaster.sendline("7") 

	for j in range(1,len(files)):
	#for j in range(1,num_images):

		#daomaster.expect(files[j])
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
	daomaster.sendline(star_name + '_on_target.mch') # these are more refined transformations
	daomaster.expect("New output file name")
	daomaster.sendline("")
	daomaster.expect("A file with the transfer table?")
	daomaster.sendline("e") # exits rest of options


	################################################################################################
	# 					CUT NUMBER OF IMAGES TO MAKE MEDIANED IMAGE
	################################################################################################

	medianed_files = []

	for wavelength in ['3p6um', '4p5um']:

		if wavelength == '3p6um':
			dither = 6
		else: dither = 1

		for epoch in range(1,6):

			for dith in range(dither, dither+5):

				filename = star_name + '_' + wavelength + '_e0' + str(epoch) + '_d' + str(dith) + '_cbcd_dn.ap'
				medianed_files.append(filename)

	# Now create copy of transformation file
	shutil.copy(star_name + '_on_target.mch', star_name + '_on_target_cut.mch')

	# Go through transformation file and delete files that do not match any filenames in medianed_files
	cut_file = star_name + '_on_target_cut.mch'
	cut_df = pd.read_csv(cut_file, delim_whitespace=True, header=None, names=['Filename', 'Apostrophe', 'A', 'B', 'C', 'D', 'E', 'F', 'Mag_offset', 'Scatter'])

	indices_to_remove = []

	for index, row in cut_df.iterrows():

		n = 0

		for filename in medianed_files:

			if filename in row['Filename']:
				n = 1 # we want to keep this row

		# If it didn't match any files then want to remove
		if n == 0:
			indices_to_remove.append(index)

	# Drop rows that are not relevant to current epoch
	cut_df.drop(indices_to_remove, inplace=True)

	# Write out to new mch file - overwrites file with ALL epochs ins
	cut_df.to_csv(cut_file, header=None, sep=' ', index=False)


	################################################################################################
	# 					MAKE MEDIANED IMAGE FROM CUT TRANSFORMATION FILE
	################################################################################################

	montage2 = pexpect.spawn('montage2')

	# Set up log file
	fout = file('montage_log.txt','w')
	montage2.logfile = fout

	montage2.expect("File with transformations:")
	montage2.sendline(star_name + '_on_target_cut.mch')
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
	montage2.sendline(star_name + '_on_target.fits')
	montage2.expect("Good bye")
	montage2.close(force=True)	

	# Now want to check whether a good medianed image was made
	# To do this, examine the weightings in the montage log file 
	# They should all be quite high. If there exists one which isn't then bad image made
	# In this case, the overly-weighted image gets discarded and a new image made
	# This repeats until all the weightings are acceptable

	# Want to keep a copy of the transformation file for all frames as this will be required later on
	# This step is just in case weightings are bad 
	mch_file = star_name + '_on_target.mch'
	shutil.copy(mch_file, star_name + '_on_target_full.mch') 

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
			with open(cut_file) as oldfile, open(cut_file+'_new', 'w') as newfile:
				for line in oldfile:
					if bad_frame not in line:
						newfile.write(line)

			shutil.copy(cut_file+'_new',  cut_file)

			# Re-run MONTAGE2
			montage2 = pexpect.spawn('montage2')

			# Set up log file
			fout = file('montage_log.txt','w')
			montage2.logfile = fout

			montage2.expect("File with transformations:")
			montage2.sendline(cut_file)
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
			montage2.sendline(star_name + '_on_target.fits')
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
	daophot.sendline("at " + star_name + '_on_target.fits')
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
	daophot.sendline(str(50)+",1") 
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
	daophot.sendline(star_name + '_on_target.coo')
	daophot.expect("Output file")
	daophot.sendline(star_name + '_on_target.ap')

	################################################################################################
	# 							MAKE PSF MODEL
	################################################################################################

	# Want to make the PSF model from a medianed image of all 240 frames as it will have a higher
	# S/N ratio and so the stars it has for PSF stars should be better

	# Make medianed image from master transformation file

	montage2 = pexpect.spawn('montage2')

	# Set up log file
	fout = file('montage_log.txt','w')
	montage2.logfile = fout

	montage2.expect("File with transformations:")
	montage2.sendline(star_name + '_on_target.mch')
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
	montage2.sendline(star_name + '_on_target_full.fits')
	montage2.expect("Good bye")
	montage2.close(force=True)	

	# Get full star list from this medianed image
	daophot = pexpect.spawn('daophot') 

	# Set up logfile
	fout = file('daophot_log.txt','w')
	daophot.logfile = fout

	# Attach medianed image and obtain star list
	daophot.expect("Command:")
	daophot.sendline("at " + star_name + '_on_target_full.fits')
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
	daophot.sendline(star_name + '_on_target_full.coo')
	daophot.expect("Are you happy with this?")
	daophot.sendline("y")

	daophot.expect("Command:")
	daophot.sendline("ph")
	daophot.expect("File with aperture radii")
	daophot.sendline("")
	daophot.expect("PHO>")
	daophot.sendline("")
	daophot.expect("Input position file")
	daophot.sendline(star_name + '_on_target_full.coo')
	daophot.expect("Output file")
	daophot.sendline(star_name + '_on_target_full.ap')	

	# Choose 20 brightest stars as candidate PSF stars
	daophot.expect("Command:")
	daophot.sendline("pi")
	daophot.expect("Input file name")
	daophot.sendline(star_name + '_on_target_full.ap')
	daophot.expect("Desired number of stars, faintest magnitude:")
	daophot.sendline("20,99")
	daophot.expect("Output file name")
	daophot.sendline(star_name + '_on_target_full.lst') 
	daophot.expect("Command:")
	daophot.sendline("ex")
	daophot.close(force=True)

	# Now run these candidate PSF stars through series of tests to get rid of bad stars

	# Read in FITS image
	hdulist = fits.open(star_name + '_on_target_full.fits')

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

	psf_stars = pd.read_csv(star_name + '_on_target_full.lst', delim_whitespace=True, skiprows=3, header=None, names=['ID', 'X', 'Y', 'Mag', 'Error'], index_col=0)

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
	f = open(star_name + '_on_target_full.lst', 'r')
	header = f.read().splitlines()[0:3]
	f.close()

	# Now overwrite this file
	f = open(star_name + '_on_target_full.lst', 'w')
	f.writelines(header[0] + '\n' + header[1] + '\n' + header[2] + '\n')

	# Send stars to lst file
	psf_stars.to_csv(f, sep=' ', mode='a', header=None)

	f.close()

	# Use these stars to create PSF model and run ALLSTAR

	daophot = pexpect.spawn('daophot')
	daophot.expect("Command:")
	daophot.sendline("at " + star_name + '_on_target_full.fits')

	daophot.expect("Command:")
	daophot.sendline("psf")
	daophot.expect("File with aperture results")
	daophot.sendline(star_name + '_on_target_full.ap')
	daophot.expect("File with PSF stars")
	daophot.sendline(star_name + '_on_target_full.lst')
	daophot.expect("File for the PSF")
	daophot.sendline(star_name + '_on_target.psf')

	daophot.expect("Command:")
	daophot.sendline("ex")
	daophot.close(force=True)

	allstar = pexpect.spawn('allstar')

	allstar.expect("OPT>")
	allstar.sendline("")
	allstar.expect("Input image name:")
	allstar.sendline(star_name + '_on_target.fits') # Run on cut medianed image as this one better at not smearing stars together
	allstar.expect("File with the PSF")
	allstar.sendline(star_name + '_on_target.psf')
	allstar.expect("Input file")
	allstar.sendline(star_name + '_on_target.ap')
	allstar.expect("File for results")
	allstar.sendline(star_name + '_on_target.als')
	allstar.expect("Name for subtracted image")
	allstar.sendline(star_name + '_on_target_dns.fits')
	allstar.expect("Good bye")
	allstar.close(force=True)

	################################################################################################
	# 						ADD OFFSETS BACK IN TO MASTER STAR LIST
	################################################################################################

	daophot = pexpect.spawn('daophot')

	daophot.expect("Command:")
	daophot.sendline("off") # offsets to put x and y back in
	daophot.expect("Input file name:")
	daophot.sendline(star_name + '_on_target.als')
	daophot.expect("Additive offsets ID, DX, DY, DMAG:")
	daophot.sendline("0," + offsets[0] + "," + offsets[1] + ",0")
	daophot.expect("Output file name")
	daophot.sendline(star_name + '_on_target.mag') # this file is the master star list 
	daophot.expect("Command:")
	daophot.sendline("ex")
	daophot.close(force=True)

	################################################################################################
	# 					TIDYING UP AND COPYING FILES TO INDIVIDUAL EPOCHS
	################################################################################################

	# Change .ap file names in .mch file to .als
	f = open(star_name +'_on_target_full.mch', 'r')
	filedata = f.read()
	f.close()

	newdata = filedata.replace(".ap",".als")

	f = open(star_name +'_on_target_full.mch','w')
	f.write(newdata)
	f.close()	

	for channel in ['1','2']:

		# Now copy the master image and the master star list to each of the epochs 
		for epoch in range(1,num_epochs+1):

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			shutil.copy(star_name + '_on_target_full.fits', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch' + channel + '/e'+epoch+'/'+star_name + '_on_target_master.fits')
			shutil.copy(star_name + '_on_target_full.mch', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch' + channel + '/e'+epoch+'/'+star_name + '_on_target_master.mch')
			shutil.copy(star_name + '_on_target.mag', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'+star_name + '_on_target_master.mag')
			shutil.copy(star_name + '_on_target.psf', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'+star_name + '_on_target_master.psf')

	# Delete temp folder
	shutil.rmtree(temp)

	return(0)

# ALLFRAME
def allframe_new(star_name, galaxy, channel, wavelength, epoch_number, field):

	if field == '1':

		start_dither = 1

		if wavelength == '3p6um':
			on_off = 'off'
		else: on_off = 'on'

	else: 

		start_dither = 6

		if wavelength == '3p6um':
			on_off = 'on'
		else: on_off = 'off'

	epoch_string = 'e' + epoch_number

	################################################################################################
	# 							MODIFY TRANSFORMATIONS FILE
	################################################################################################

	# The transformation file currently has transformations for ALL epochs and dithers and both channels for on target
	# For off target field, transformation file only contains the 5 frames from that epoch for that channel
	# We want to remove frames not relevant to the current epoch
	# However, we want to keep e01 d1 (for field 1) or d6 (for field 2) as this is the 'master' input file

	# Open mch file i.e. the transformation file and delete rows not relevant to current epoch

	mch_file = star_name + '_' + on_off + '_target_master.mch'
	mch_df = pd.read_csv(mch_file, delim_whitespace=True, header=None, names=['Filename', 'Apostrophe', 'A', 'B', 'C', 'D', 'E', 'F', 'Mag_offset', 'Scatter'])

	indices_to_remove = []

	for index, row in mch_df.iterrows():

		if index != 0: # we don't get rid of first line as this is the 'master' input file

			# Obtain list of indices that contain files we want to remove as not relevant to current epoch
			if epoch_string not in row['Filename'] or wavelength not in row['Filename']:
				indices_to_remove.append(index)

	# Drop rows that are not relevant to current epoch
	mch_df.drop(indices_to_remove, inplace=True)

	# Write out to new mch file - overwrites file with ALL epochs ins
	mch_df.to_csv(mch_file, header=None, sep=' ', index=False)


	################################################################################################
	# 					COPY MASTER PSF MODEL TO INDIVIDUAL DITHER NAMES
	################################################################################################

	# This step is required by ALLFRAME as it requires a PSF model for each dither
	for dither in range(start_dither, start_dither+5):
		shutil.copy(star_name + '_' + on_off + '_target_master.psf', star_name + '_' + wavelength + '_' + epoch_string + '_d' + str(dither) + '_cbcd_dn.psf')	


	################################################################################################
	# 							COPY NECESSARY EPOCH 1 FILES
	################################################################################################

	# If on target field then copy across necessary files from 3p6 e01 d6 as this is the master input file
	if epoch_string != 'e01' and on_off == 'on' and wavelength == '3p6um':
		shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch1/e01/' + star_name + '_3p6um_e01_d6_cbcd_dn.ap', star_name + '_3p6um_e01_d6_cbcd_dn.ap')
		shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch1/e01/' + star_name + '_3p6um_e01_d6_cbcd_dn.als', star_name + '_3p6um_e01_d6_cbcd_dn.als')
		shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch1/e01/' + star_name + '_3p6um_e01_d6_cbcd_dn.psf', star_name + '_3p6um_e01_d6_cbcd_dn.psf')
		shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch1/e01/' + star_name + '_3p6um_e01_d6_cbcd_dn.fits', star_name + '_3p6um_e01_d6_cbcd_dn.fits')

	if on_off == 'on' and wavelength == '4p5um':
		shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch1/e01/' + star_name + '_3p6um_e01_d6_cbcd_dn.ap', star_name + '_3p6um_e01_d6_cbcd_dn.ap')
		shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch1/e01/' + star_name + '_3p6um_e01_d6_cbcd_dn.als', star_name + '_3p6um_e01_d6_cbcd_dn.als')
		shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch1/e01/' + star_name + '_3p6um_e01_d6_cbcd_dn.psf', star_name + '_3p6um_e01_d6_cbcd_dn.psf')
		shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch1/e01/' + star_name + '_3p6um_e01_d6_cbcd_dn.fits', star_name + '_3p6um_e01_d6_cbcd_dn.fits')


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
	allframe.sendline(star_name + '_' + on_off + '_target_master.mag')
	allframe.expect("Good bye")
	allframe.close(force=True)

	# Finally, copy e01 3p6 d6 .alf file across as will need it later on
	if epoch_number != '01':
		shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch1/e01/' + star_name + '_3p6um_e01_d6_cbcd_dn.alf', star_name + '_3p6um_e01_d6_cbcd_dn.alf')


	################################################################################################
	# 								CHECK ALLFRAME OUTPUTS
	################################################################################################

	# Check there are no empty .alf files and if there is, add one fake star so that subsequent scripts work

	for dither in range(start_dither, start_dither+5):

		alf_file = star_name + '_' + wavelength + '_' + epoch_string + '_d' + str(dither) + '_cbcd_dn.alf'

		num_lines = sum(1 for line in open(alf_file))

		if num_lines <= 3:

			# Add fake star to alf file
			f = open(alf_file, 'a')
			f.write("10000   50.000  50.000   19.000   0.999     0.26       5.     0.87   -0.021 \n")
			f.close()

	return(0)

# Remove PSF stars from .lst which have stars in annulus
def remove_neighbours(star_name, galaxy, wavelength, epoch_number, dither):

	# Get filenames 
	lst_file = star_name + '_' + wavelength + '_e' + epoch_number + '_d' + str(dither) + '_cbcd_dn.lst'
	coo_file = star_name + '_' + wavelength + '_e' + epoch_number + '_d' + str(dither) + '_cbcd_dn.coo'

	# Open data into df's
	lst_df = pd.read_csv(lst_file, delim_whitespace=True, header=None, skiprows=3, names=['ID', 'X', 'Y', 'mag', 'err'] )
	coo_df = pd.read_csv(coo_file, delim_whitespace=True, header=None, skiprows=3, names=['ID', 'X', 'Y', 'mag', 'err'], usecols=[0,1,2,3,4] )

	# Create empty list for PSF stars to be appended to that need to be removed
	PSF_stars_to_remove = []

	# Iterate over each PSF star, determining whether it has a star in its annulus
	for index, row in lst_df.iterrows():
	    
	    #print "PSF star ID = " + str(lst_df['ID'][index])
	    
	    neighbour = 0 # set initial counter to 0 i.e. no neighbour

	    # Get x and y coordinates of PSF star
	    PSF_x = lst_df['X'][index]
	    PSF_y = lst_df['Y'][index]
	    
	    # Iterate over each star in .coo file
	    for neighbour_index, neighbour_row in coo_df.iterrows():
	        
	        # Get x and y coordinates of potential neighbour star
	        Nei_x = coo_df['X'][neighbour_index]
	        Nei_y = coo_df['Y'][neighbour_index]
	        
	        # Calculate distance between the two stars
	        d = calc_distance(PSF_x, PSF_y, Nei_x, Nei_y)
	        
	        # If distance between is=12 and os=20 then star is in annulus
	        # neighbour gets set to 1 and code loops out as there is at least one neighbour
	        if d>=12 and d<=20:
	            
	            #print "Neighbour with ID = %d inside annulus" % coo_df['ID'][neighbour_index]
	            neighbour = 1 # update counter
	            break # don't need to test anymore stars
	    
	    # Either keep or remove PSF star based on whether there is a neighbour
	    if neighbour == 1:    
	        PSF_stars_to_remove.append(index)

	# Drop rows with indices in PSF_stars_to_remove
	rem_df = lst_df.drop(PSF_stars_to_remove, axis=0)		    

	# # In the case that there exists at least one star without a neighbour, write to file
	# if len(rem_df) != 0:

	# 	# Output remaining stars into .lst file format
	# 	f = open(lst_file)

	# 	# Write the 3 line header to a new file
	# 	header = f.read().splitlines()[0:3]
	# 	g = open(lst_file, 'w')
	# 	g.writelines(header[0] + '\n' + header[1] + '\n' + header[2] + '\n')
		  
	# 	g.close()
	# 	f.close()

	# 	# Write remaining PSF stars to file
	# 	rem_df.to_csv(lst_file, mode='a', header=None, sep=' ', index=False)	
	

	# else:

	# 	print "Wavelength = %s  Epoch = %s  Dither = %s  has no PSF stars without neighbours. Using the brightest 5." % (str(wavelength), str(epoch_number), str(dither))

	# 	# Drop all rows except the first 5 as these are the brightest stars
	# 	lst_df.drop(lst_df.index[5:], axis=0, inplace=True)	


	# 	# Output remaining stars into .lst file format
	# 	f = open(lst_file)

	# 	# Write the 3 line header to a new file
	# 	header = f.read().splitlines()[0:3]
	# 	g = open(lst_file, 'w')
	# 	g.writelines(header[0] + '\n' + header[1] + '\n' + header[2] + '\n')
		  
	# 	g.close()
	# 	f.close()

	# 	# Write remaining PSF stars to file
	# 	lst_df.to_csv(lst_file, mode='a', header=None, sep=' ', index=False)


	# Output remaining stars into .lst file format
	f = open(lst_file)

	# Write the 3 line header to a new file
	header = f.read().splitlines()[0:3]
	g = open(lst_file, 'w')
	g.writelines(header[0] + '\n' + header[1] + '\n' + header[2] + '\n')
	  
	g.close()
	f.close()

	# Write remaining PSF stars to file
	rem_df.to_csv(lst_file, mode='a', header=None, sep=' ', index=False)	

	return(0)


# Calculate distance between two points
def calc_distance(x1,y1,x2,y2):
	return (sqrt((x2-x1)**2 + (y2-y1)**2))

# Calibration procedure onto IRAC Vega system
def calibration_procedure_3_3_7(star_name, galaxy, channel, wavelength, epoch_number, field):

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

	# This list contains all magnitude differences for a given field and will be averaged over to get aperture correction value
	mag_differences = []

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

		# Delete any rows which have a 99.9999 psf mag
		df_combined.drop(df_combined[df_combined.psf_mag == 99.9999].index, inplace=True)

		if len(df_combined) == 0:
			print "%s - No PSF stars with non-99.9999 values" % alf

		################################################################################################
		# 							CALCULATE APERTURE CORRECTION
		################################################################################################

		# Calculate the difference between aperture and PSF magnitudes
		df_combined['mag_difference'] = df_combined['ap_mag'] - df_combined['psf_mag']	

		# # Calculate aperture correction value by taking average
		# apc = round(df_combined['mag_difference'].mean(), 3)

		# Append every value in df_combined['mag_difference'] to list mag_differences
		for i in range(0,len(df_combined)):

			mag_differences.append(df_combined['mag_difference'][i])



	################################################################################################
	# 							CALCULATE AVERAEGE APERTURE CORRECTION
	################################################################################################

	mag_differences = np.array(mag_differences) # convert to np array

	# Calculate aperture correction value that will be used for this whole field
	apc = round(np.average(mag_differences),3)


	for alf in alf_files:


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
	# sy'/home/ac833/Data/' used by Reach et al. (2005) and is the standard IRAC Vega system.
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

		# To go from (3,3,7) to (10,12,20) we have a correction value from IRAC handbook
		# The values are different for ch1 and ch2
		std_corr_3p6 = 1.125
		std_corr_4p5 = 1.120

		# Put alf_apc magnitudes into dataframe
		df_alf = pd.read_csv(alf, header=0, delim_whitespace=True)

		################################################################################################
		# 								APPLY CORRECTIONS
		################################################################################################

		for index in range(0, len(df_alf)):

			# Convert magnitude to flux
			flux =  F0 * (10 ** (-df_alf['mag'][index]/2.5))

			# Apply std aperture correction value
			if wavelength == '3p6um':
				flux *= std_corr_3p6
			else: flux *= std_corr_4p5

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

	# Set parameters for aperture = 3, inner sky = 3, outer sky = 7
	if channel == '1':
		x0 = 0.186
		y0 = 0.021
		sig_x = 0.191
		sig_y = 0.164
		delf_x = 0.0330
		delf_y = 0.0540
		f0 = 0.962
	elif channel == '2':
		x0 = 0.071
		y0 = 0.099
		sig_x = 0.201
		sig_y = 0.218
		delf_x = 0.0179
		delf_y = 0.0192
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

# Master image on target
def master_on_target_no_timeout(star_name, galaxy, num_epochs):

	# if channel == '1':
	# 	field = '2' # the on-target images are the second 5 dithers d6-10
	# 	start_dither = 6
	# else:
	# 	field = '1' # the on-target images are the first 5 dithers d1-5
	# 	start_dither = 1

	num_images = num_epochs * 5 * 2 # 5 dithers per epoch for both channels

	# List of files to use
	files = []

	# All 3p6 on target files
	for epoch in range(1,num_epochs+1):
		for dither in range(6,11): 

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			dither = str(dither)

			filename = star_name + '_3p6um_e' + epoch + '_d' + dither + '_cbcd_dn.ap'
			files.append(filename)	

	# All 4p5 on target files
	for epoch in range(1,num_epochs+1):
		for dither in range(1,6): 

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			dither = str(dither)

			filename = star_name + '_4p5um_e' + epoch + '_d' + dither + '_cbcd_dn.ap'
			files.append(filename)	

	# Create temporary folder to make image - will then be copied to each epoch folder
	temp = '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/temp/'

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

		cwd = '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch1/e'+epoch+'/'

		# Copy 3p6 files
		for dither in range(6,11):

			dither = str(dither)

			shutil.copyfile(cwd + star_name + '_3p6um_e' + epoch + '_d' + dither + '_cbcd_dn.ap', temp + star_name + '_3p6um_e' + epoch + '_d' + dither + '_cbcd_dn.ap') #.ap
			shutil.copyfile(cwd + star_name + '_3p6um_e' + epoch + '_d' + dither + '_cbcd_dn.fits', temp + star_name + '_3p6um_e' + epoch + '_d' + dither + '_cbcd_dn.fits')

	# Copy all FITS images and aperture photometry files to temp folder
	for epoch in range(1,num_epochs+1):

		# Get epoch in correct format
		if epoch < 10:
			epoch = '0' + str(epoch)
		else: epoch = str(epoch)

		cwd = '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch2/e'+epoch+'/'

		# Copy 3p6 files
		for dither in range(1,6):

			dither = str(dither)

			shutil.copyfile(cwd + star_name + '_4p5um_e' + epoch + '_d' + dither + '_cbcd_dn.ap', temp + star_name + '_4p5um_e' + epoch + '_d' + dither + '_cbcd_dn.ap') #.ap
			shutil.copyfile(cwd + star_name + '_4p5um_e' + epoch + '_d' + dither + '_cbcd_dn.fits', temp + star_name + '_4p5um_e' + epoch + '_d' + dither + '_cbcd_dn.fits')

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
	daomatch.sendline(files[0]) 
	daomatch.expect("Output file name")
	daomatch.sendline(star_name + '_on_target.mch')	

	index = 0
	master = 0

	# Give it the rest of the images
	for j in range(1,len(files)):

		print '----------------------------------'
		print files[j]

		daomatch.sendline(files[j]+'/')

		index = daomatch.expect(["Next input file", "Write this transformation?"])
		#index = daomatch.expect(["Next input file", "Write this transformation?", pexpect.exceptions.TIMEOUT])
		#print "Index = " + str(index)

		print index

		if index == 1:

			master = 1

			daomatch.sendline('Y')
			#daomatch.sendline(files[j]+'/')

		# if index == 2:

		# 	daomatch.sendline('Y')

		# This makes sure the last file is read properly
		if j == len(files)-1:

			daomatch.sendline("")		


		#reset index
		index = 0
		print index

		# print "Index = " + str(index)

		# # If index = 0, then we got "Next input file" so matching worked
		# # Do nothing
		# if index == 1:

		# 	print index

		# 	print files[j] + " could not match properly"

		# 	#daomatch.expect('Write this transformation')
		# 	daomatch.sendline('No')

		# # If index = 1, then we got "Write this transformation" so matching bad
		# else:

		# 	print index

		# 	print files[j] + " matched correctly"

		# 	if j == len(files)-1:

		# 		daomatch.sendline("")


		# # Index will tell you whether DAOMATCH could macth the file or not
		# index = daomatch.expect(["Next input file", pexpect.exceptions.TIMEOUT])


		#daomatch.expect("Next input file")
		#daomatch.sendline(files[j]+'/')

	#daomatch.expect("Next input file")
	daomatch.sendline("") # exit

	################################################################################################
	# 							REMOVE DUPLICATE LINES IN MATCH FILE
	################################################################################################

	lines_seen = set() # holds lines already seen
	outfile = open('temp.mch', "w")
	for line in open(star_name+'_on_target.mch', "r"):
	    if line not in lines_seen: # not a duplicate
	        outfile.write(line)
	        lines_seen.add(line)
	outfile.close()

	shutil.copy('temp.mch', star_name + '_on_target.mch')


	################################################################################################
	# 									RUN DAOMASTER
	################################################################################################


	if master == 0:

		# Use DAOMASTER to refine the transformations
		daomaster = pexpect.spawn('daomaster')

		fout = file('daomaster_log.txt','w')
		daomaster.logfile = fout

		daomaster.expect("File with list of input files:")
		daomaster.sendline(star_name + '_on_target.mch')
		daomaster.expect("Minimum number, minimum fraction, enough frames:")
		daomaster.sendline("1, 0.5, " + str(num_images)) # play around with these values
		daomaster.expect("Maximum sigma:")
		daomaster.sendline("99") 
		daomaster.expect("Your choice:")
		daomaster.sendline("6") # solve for 6 degrees of freedom
		daomaster.expect("Critical match-up radius:")
		daomaster.sendline("7") 

		for j in range(1,len(files)):
		#for j in range(1,num_images):

			#daomaster.expect(files[j])
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
		daomaster.sendline(star_name + '_on_target.mch') # these are more refined transformations
		daomaster.expect("New output file name")
		daomaster.sendline("")
		daomaster.expect("A file with the transfer table?")
		daomaster.sendline("e") # exits rest of options


	################################################################################################
	# 					CUT NUMBER OF IMAGES TO MAKE MEDIANED IMAGE
	################################################################################################

	medianed_files = []

	for wavelength in ['3p6um', '4p5um']:

		if wavelength == '3p6um':
			dither = 6
		else: dither = 1

		for epoch in range(1,6):

			for dith in range(dither, dither+5):

				filename = star_name + '_' + wavelength + '_e0' + str(epoch) + '_d' + str(dith) + '_cbcd_dn.ap'
				medianed_files.append(filename)

	# Now create copy of transformation file
	shutil.copy(star_name + '_on_target.mch', star_name + '_on_target_cut.mch')

	# Go through transformation file and delete files that do not match any filenames in medianed_files
	cut_file = star_name + '_on_target_cut.mch'
	cut_df = pd.read_csv(cut_file, delim_whitespace=True, header=None, names=['Filename', 'Apostrophe', 'A', 'B', 'C', 'D', 'E', 'F', 'Mag_offset', 'Scatter'])

	indices_to_remove = []

	for index, row in cut_df.iterrows():

		n = 0

		for filename in medianed_files:

			if filename in row['Filename']:
				n = 1 # we want to keep this row

		# If it didn't match any files then want to remove
		if n == 0:
			indices_to_remove.append(index)

	# Drop rows that are not relevant to current epoch
	cut_df.drop(indices_to_remove, inplace=True)

	# Write out to new mch file - overwrites file with ALL epochs ins
	cut_df.to_csv(cut_file, header=None, sep=' ', index=False)


	################################################################################################
	# 					MAKE MEDIANED IMAGE FROM CUT TRANSFORMATION FILE
	################################################################################################

	montage2 = pexpect.spawn('montage2')

	# Set up log file
	fout = file('montage_log.txt','w')
	montage2.logfile = fout

	montage2.expect("File with transformations:")
	montage2.sendline(star_name + '_on_target_cut.mch')
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
	montage2.sendline(star_name + '_on_target.fits')
	montage2.expect("Good bye")
	montage2.close(force=True)	

	# Now want to check whether a good medianed image was made
	# To do this, examine the weightings in the montage log file 
	# They should all be quite high. If there exists one which isn't then bad image made
	# In this case, the overly-weighted image gets discarded and a new image made
	# This repeats until all the weightings are acceptable

	# Want to keep a copy of the transformation file for all frames as this will be required later on
	# This step is just in case weightings are bad 
	mch_file = star_name + '_on_target.mch'
	shutil.copy(mch_file, star_name + '_on_target_full.mch') 

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
			with open(cut_file) as oldfile, open(cut_file+'_new', 'w') as newfile:
				for line in oldfile:
					if bad_frame not in line:
						newfile.write(line)

			shutil.copy(cut_file+'_new',  cut_file)

			# Re-run MONTAGE2
			montage2 = pexpect.spawn('montage2')

			# Set up log file
			fout = file('montage_log.txt','w')
			montage2.logfile = fout

			montage2.expect("File with transformations:")
			montage2.sendline(cut_file)
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
			montage2.sendline(star_name + '_on_target.fits')
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
	daophot.sendline("at " + star_name + '_on_target.fits')
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
	daophot.sendline(str(50)+",1") 
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
	daophot.sendline(star_name + '_on_target.coo')
	daophot.expect("Output file")
	daophot.sendline(star_name + '_on_target.ap')

	################################################################################################
	# 							MAKE PSF MODEL
	################################################################################################

	# Want to make the PSF model from a medianed image of all 240 frames as it will have a higher
	# S/N ratio and so the stars it has for PSF stars should be better

	# Make medianed image from master transformation file

	montage2 = pexpect.spawn('montage2')

	# Set up log file
	fout = file('montage_log.txt','w')
	montage2.logfile = fout

	montage2.expect("File with transformations:")
	montage2.sendline(star_name + '_on_target.mch')
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
	montage2.sendline(star_name + '_on_target_full.fits')
	montage2.expect("Good bye")
	montage2.close(force=True)	

	# Get full star list from this medianed image
	daophot = pexpect.spawn('daophot') 

	# Set up logfile
	fout = file('daophot_log.txt','w')
	daophot.logfile = fout

	# Attach medianed image and obtain star list
	daophot.expect("Command:")
	daophot.sendline("at " + star_name + '_on_target_full.fits')
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
	daophot.sendline(star_name + '_on_target_full.coo')
	daophot.expect("Are you happy with this?")
	daophot.sendline("y")

	daophot.expect("Command:")
	daophot.sendline("ph")
	daophot.expect("File with aperture radii")
	daophot.sendline("")
	daophot.expect("PHO>")
	daophot.sendline("")
	daophot.expect("Input position file")
	daophot.sendline(star_name + '_on_target_full.coo')
	daophot.expect("Output file")
	daophot.sendline(star_name + '_on_target_full.ap')	

	# Choose 20 brightest stars as candidate PSF stars
	daophot.expect("Command:")
	daophot.sendline("pi")
	daophot.expect("Input file name")
	daophot.sendline(star_name + '_on_target_full.ap')
	daophot.expect("Desired number of stars, faintest magnitude:")
	daophot.sendline("20,99")
	daophot.expect("Output file name")
	daophot.sendline(star_name + '_on_target_full.lst') 
	daophot.expect("Command:")
	daophot.sendline("ex")
	daophot.close(force=True)

	# Now run these candidate PSF stars through series of tests to get rid of bad stars

	# Read in FITS image
	hdulist = fits.open(star_name + '_on_target_full.fits')

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

	psf_stars = pd.read_csv(star_name + '_on_target_full.lst', delim_whitespace=True, skiprows=3, header=None, names=['ID', 'X', 'Y', 'Mag', 'Error'], index_col=0)

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
	f = open(star_name + '_on_target_full.lst', 'r')
	header = f.read().splitlines()[0:3]
	f.close()

	# Now overwrite this file
	f = open(star_name + '_on_target_full.lst', 'w')
	f.writelines(header[0] + '\n' + header[1] + '\n' + header[2] + '\n')

	# Send stars to lst file
	psf_stars.to_csv(f, sep=' ', mode='a', header=None)

	f.close()

	# Use these stars to create PSF model and run ALLSTAR

	daophot = pexpect.spawn('daophot')
	daophot.expect("Command:")
	daophot.sendline("at " + star_name + '_on_target_full.fits')

	daophot.expect("Command:")
	daophot.sendline("psf")
	daophot.expect("File with aperture results")
	daophot.sendline(star_name + '_on_target_full.ap')
	daophot.expect("File with PSF stars")
	daophot.sendline(star_name + '_on_target_full.lst')
	daophot.expect("File for the PSF")
	daophot.sendline(star_name + '_on_target.psf')

	daophot.expect("Command:")
	daophot.sendline("ex")
	daophot.close(force=True)

	allstar = pexpect.spawn('allstar')

	allstar.expect("OPT>")
	allstar.sendline("")
	allstar.expect("Input image name:")
	allstar.sendline(star_name + '_on_target.fits') # Run on cut medianed image as this one better at not smearing stars together
	allstar.expect("File with the PSF")
	allstar.sendline(star_name + '_on_target.psf')
	allstar.expect("Input file")
	allstar.sendline(star_name + '_on_target.ap')
	allstar.expect("File for results")
	allstar.sendline(star_name + '_on_target.als')
	allstar.expect("Name for subtracted image")
	allstar.sendline(star_name + '_on_target_dns.fits')
	allstar.expect("Good bye")
	allstar.close(force=True)

	################################################################################################
	# 						ADD OFFSETS BACK IN TO MASTER STAR LIST
	################################################################################################

	daophot = pexpect.spawn('daophot')

	daophot.expect("Command:")
	daophot.sendline("off") # offsets to put x and y back in
	daophot.expect("Input file name:")
	daophot.sendline(star_name + '_on_target.als')
	daophot.expect("Additive offsets ID, DX, DY, DMAG:")
	daophot.sendline("0," + offsets[0] + "," + offsets[1] + ",0")
	daophot.expect("Output file name")
	daophot.sendline(star_name + '_on_target.mag') # this file is the master star list 
	daophot.expect("Command:")
	daophot.sendline("ex")
	daophot.close(force=True)

	################################################################################################
	# 					TIDYING UP AND COPYING FILES TO INDIVIDUAL EPOCHS
	################################################################################################

	# Change .ap file names in .mch file to .als
	f = open(star_name +'_on_target_full.mch', 'r')
	filedata = f.read()
	f.close()

	newdata = filedata.replace(".ap",".als")

	f = open(star_name +'_on_target_full.mch','w')
	f.write(newdata)
	f.close()	

	for channel in ['1','2']:

		# Now copy the master image and the master star list to each of the epochs 
		for epoch in range(1,num_epochs+1):

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			shutil.copy(star_name + '_on_target_full.fits', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch' + channel + '/e'+epoch+'/'+star_name + '_on_target_master.fits')
			shutil.copy(star_name + '_on_target_full.mch', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch' + channel + '/e'+epoch+'/'+star_name + '_on_target_master.mch')
			shutil.copy(star_name + '_on_target.mag', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'+star_name + '_on_target_master.mag')
			shutil.copy(star_name + '_on_target.psf', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'+star_name + '_on_target_master.psf')

	# Delete temp folder
	#shutil.rmtree(temp)

	return(0)


# Aperture correction using all stars for one star from all epochs
def aperture_correction(star_name, galaxy, channel, wavelength):

	if galaxy == 'LMC':
		num_epochs = 24
	else: num_epochs = 12

	# Set up list of files to obtain stars from
	als_files = [] # file containing all stars and PSF mags
	lst_files = [] # file containing PSF stars and aperture mags
	alf_files = [] # the magnitudes that need to be corrected

	for epoch in range(1,num_epochs+1):

		if epoch < 10:
			epoch_number = '0' + str(epoch)
		else: epoch_number = str(epoch)

		for dither in range(1,11):

			als_files.append('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/' + star_name + '_' + wavelength + '_e' + epoch_number + '_d' + str(dither) + '_cbcd_dn.als')
			lst_files.append('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/' + star_name + '_' + wavelength + '_e' + epoch_number + '_d' + str(dither) + '_cbcd_dn.lst')
			alf_files.append('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/' + star_name + '_' + wavelength + '_e' + epoch_number + '_d' + str(dither) + '_cbcd_dn.alf')

	# List of all mag differences between aperture and psf mags
	mag_differences = []

	# Open each lst/als file combo, cross-match stars and calculate aperture correction value
	for index in range(0,len(als_files)):

		als_file = als_files[index]
		lst_file = lst_files[index]
		alf_file = alf_files[index]

		# Open files into dataframes
		# The lst dataframe contains the aperture magnitudes and ALL stars will be used from this file
		# The als dataframe contains the psf magnitudes and only a SUBSET of these stars with IDs matching those of the lst df will be used
		df_lst = pd.read_csv(lst_file, skiprows=3, header=None, names=['ID_lst', 'X_lst', 'Y_lst', 'ap_mag', 'ap_error'], usecols=(0,1,2,3,4), delim_whitespace=True)
		df_als = pd.read_csv(als_file, skiprows=3, header=None, names=['ID_als', 'X_als', 'Y_als', 'psf_mag', 'psf_error'], usecols=(0,1,2,3,4), delim_whitespace=True)

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

		# Delete any rows which have a 99.9999 psf mag
		df_combined.drop(df_combined[df_combined.psf_mag == 99.9999].index, inplace=True)

		if len(df_combined) == 0:
			print "%s - No PSF stars with non-99.9999 values" % alf_file	

		################################################################################################
		# 							CALCULATE APERTURE CORRECTION
		################################################################################################

		# Calculate the difference between aperture and PSF magnitudes
		df_combined['mag_difference'] = df_combined['ap_mag'] - df_combined['psf_mag']	

		# # Calculate aperture correction value by taking average
		# apc = round(df_combined['mag_difference'].mean(), 3)

		# Append every value in df_combined['mag_difference'] to list mag_differences
		for i in range(0,len(df_combined)):

			mag_differences.append(df_combined['mag_difference'][i])


	################################################################################################
	# 							CALCULATE AVERAEGE APERTURE CORRECTION
	################################################################################################

	mag_differences = np.array(mag_differences) # convert to np array

	if len(mag_differences) != 0:

		# Calculate aperture correction value that will be used for this whole field
		apc = round(np.average(mag_differences),3)

	else:

		print "Wavelength = %s  Epoch = %s  Field = %s has no stars for aperture correction." % (str(wavelength), str(epoch_number), str(field))
		apc = 0


	for alf in alf_files:


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


	return(0)

# Calibration procedure onto IRAC Vega system
def calibration_procedure_testing(star_name, galaxy, channel, wavelength, epoch_number, field):

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


	# Copy apc file to cal file extension so next part of code works
	for alf in alf_files:

		alf_apc = alf.replace('.alf', '.alf_apc')
		#alf_cal = alf.replace('.alf', '.alf_cal')

		alf_df = pd.read_csv(alf, skiprows=3, usecols=[0,1,2,3,4], names=['ID', 'x', 'y', 'mag', 'error'], delim_whitespace=True)
		alf_df.to_csv(alf_apc, sep=' ', index=False)


	# ################################################################################################
	# # 									APERTURE CORRECTION
	# ################################################################################################

	# # The aperture correction is applied to the PSF magnitudes from ALLFRAME (i.e. the .alf files) 
	# # To calculate the aperture correction, the aperture magnitudes (.ap) and PSF magnitudes (.als) for a 
	# # bright list of stars is used to calculate the average offset.
	# # This average offset is then applied to the PSF magnitudes.

	# # This list contains all magnitude differences for a given field and will be averaged over to get aperture correction value
	# mag_differences = []

	# for alf in alf_files:

	# 	################################################################################################
	# 	# 									SET UP DATAFRAMES
	# 	################################################################################################

	# 	# Get filenames for the two files that need comparing 
	# 	als = alf.replace('.alf', '.als')
	# 	ap = alf.replace('.alf', '.ap')

	# 	# Image name
	# 	image = alf.replace('.alf', '.fits')

	# 	# Bright stars that are used to compare magnitude difference
	# 	lst = alf.replace('.alf', '.lst')

	# 	# Open files into dataframes
	# 	# The lst dataframe contains the aperture magnitudes and ALL stars will be used from this file
	# 	# The als dataframe contains the psf magnitudes and only a SUBSET of these stars with IDs matching those of the lst df will be used
	# 	df_lst = pd.read_csv(lst, skiprows=3, header=None, names=['ID_lst', 'X_lst', 'Y_lst', 'ap_mag', 'ap_error'], usecols=(0,1,2,3,4), delim_whitespace=True)
	# 	df_als = pd.read_csv(als, skiprows=3, header=None, names=['ID_als', 'X_als', 'Y_als', 'psf_mag', 'psf_error'], usecols=(0,1,2,3,4), delim_whitespace=True)

	# 	# Sort rows by ID number
	# 	df_lst.sort_values('ID_lst', axis=0, inplace=True)
	# 	df_als.sort_values('ID_als', axis=0, inplace=True)

	# 	# Reset index to start at 0
	# 	df_lst.reset_index(drop=True, inplace=True)
	# 	df_als.reset_index(drop=True, inplace=True)	

	# 	# Crossmatch the two dataframes to only keep stars that appear in both dataframes
	# 	for i in range(0, len(df_als)):

	# 		match = 0 

	# 		for j in range(0, len(df_lst)):

	# 			# Check whether IDs match
	# 			if df_als['ID_als'][i] == df_lst['ID_lst'][j]:
	# 				match += 1

	# 		# If they match, do nothing as want to keep this row. Otherwise, delete the row in the als dataframe
	# 		if match == 0:
	# 			df_als.drop(i, axis=0, inplace=True)

	# 	# Reset index to start at 0
	# 	df_lst.reset_index(drop=True, inplace=True)
	# 	df_als.reset_index(drop=True, inplace=True)

	# 	# Now do the same process but removing any stars in the lst file not in the als file as the stars need to have magnitudes in both dataframes
	# 	for i in range(0, len(df_lst)):

	# 		match = 0

	# 		for j in range(0, len(df_als)):

	# 			# Check whether IDs match
	# 			if df_lst['ID_lst'][i] == df_als['ID_als'][j]:
	# 				match += 1

	# 		# If they match, do nothing as want to keep this row. Otherwise, delete the row in the als dataframe
	# 		if match == 0:
	# 			df_lst.drop(i, axis=0, inplace=True)

	# 	# Again, reset index to start at 0
	# 	df_lst.reset_index(drop=True, inplace=True)
	# 	df_als.reset_index(drop=True, inplace=True)	

	# 	# Now df_als and df_lst have the same stars so can concatenate
	# 	df_combined = pd.concat((df_als, df_lst), axis=1)

	# 	# Remove columns not needed
	# 	df_combined.drop(['X_lst', 'Y_lst', 'ID_lst'], axis=1, inplace=True)

	# 	# Delete any rows which have a 99.9999 psf mag
	# 	df_combined.drop(df_combined[df_combined.psf_mag == 99.9999].index, inplace=True)

	# 	if len(df_combined) == 0:
	# 		print "%s - No PSF stars with non-99.9999 values" % alf

	# 	################################################################################################
	# 	# 							CALCULATE APERTURE CORRECTION
	# 	################################################################################################

	# 	# Calculate the difference between aperture and PSF magnitudes
	# 	df_combined['mag_difference'] = df_combined['ap_mag'] - df_combined['psf_mag']	

	# 	# # Calculate aperture correction value by taking average
	# 	# apc = round(df_combined['mag_difference'].mean(), 3)

	# 	# Append every value in df_combined['mag_difference'] to list mag_differences
	# 	for i in range(0,len(df_combined)):

	# 		mag_differences.append(df_combined['mag_difference'][i])



	# ################################################################################################
	# # 							CALCULATE AVERAEGE APERTURE CORRECTION
	# ################################################################################################

	# mag_differences = np.array(mag_differences) # convert to np array

	# if len(mag_differences) != 0:

	# 	# Calculate aperture correction value that will be used for this whole field
	# 	apc = round(np.average(mag_differences),3)

	# else:

	# 	print "Wavelength = %s  Epoch = %s  Field = %s has no stars for aperture correction." % (str(wavelength), str(epoch_number), str(field))
	# 	apc = 0


	# for alf in alf_files:


	# 	################################################################################################
	# 	# 					APPLY APERTURE CORRECTION TO ALLFRAME MAGNITUDES
	# 	################################################################################################

	# 	# This step converts the PSF magnitudes onto the (3,12,20) aperture system

	# 	df_alf = pd.read_csv(alf, header=None, skiprows=3, usecols=[0,1,2,3,4], names=['ID', 'x', 'y', 'mag', 'error'], delim_whitespace=True)
	# 	df_alf['mag'] += apc # ADD apc value to magnitudes


	# 	################################################################################################
	# 	# 						WRITE CORRECTED MAGNITUDES TO NEW FILE
	# 	################################################################################################

	# 	new_filename = alf.replace('.alf', '.alf_apc')
	# 	df_alf.to_csv(new_filename, sep=' ', index=False)


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
			#std_corr = 1.060
		elif channel == '2':
			F0 = 179.7
			#std_corr = 1.063
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

	# # Set parameters for aperture = 5, inner sky = 5, outer sky = 10
	# if channel == '1':
	# 	x0 = 0.191
	# 	y0 = 0.027
	# 	sig_x = 0.190
	# 	sig_y = 0.163
	# 	delf_x = 0.0298
	# 	delf_y = 0.0479
	# 	f0 = 0.966
	# elif channel == '2':
	# 	x0 = 0.091
	# 	y0 = 0.091
	# 	sig_x = 0.209
	# 	sig_y = 0.209
	# 	delf_x = 0.0176
	# 	delf_y = 0.0183
	# 	f0 = 0.981
	# else: x0 = 'break'	


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
