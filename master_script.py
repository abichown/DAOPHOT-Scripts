'''
Purpose: This is the master script that performs the entire photometry process and produces GLOESS curves.
The input for this is a list of stars that contains RA, Dec, Period, Channel and Star Name.
This script assumes that the scripts file_setup.py and convert_to_counts.py have already been executed.
The functions called upon here can be found in the functions.py script.
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

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
#from kneed import DataGenerator, KneeLocator

import matplotlib
import gloess_fits as gf
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
matplotlib.rc('text',usetex=True)
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Garamond']

from functions import aper_phot, master_on_target, master_off_target, psf_phot, allframe, calibration_procedure, combine_dithers, ave_mag, get_mag, format_gloess, gloess_single_band

sys.settrace('line')

start = time.time()

#####################################################################################################
# 										READ IN STAR LIST
#####################################################################################################


# Set up input file into a dataframe
df = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star','Period', 'RA', 'Dec', 'Channel'])


# Iterate over every star in the df to carry out the pipeline
for i in range(0, len(df)):


	#####################################################################################################
	# 									OBTAIN INFORMATION ON STAR 
	#####################################################################################################


	galaxy = str(df['Galaxy'][i])
	star_name = str(df['Star'][i])
	ra = df['RA'][i]
	dec = df['Dec'][i]
	period = df['Period'][i]

	if df['Channel'][i] == 1:
		channel = '1'
		wavelength = '3p6um'
	elif df['Channel'][i] == 2:
		channel = '2'
		wavelength = '4p5um'
	else: wavelength = 'channel not defined'

	if galaxy == 'LMC':
		num_epochs = 24
	elif galaxy == 'SMC':
		num_epochs = 12
	else: num_epochs = 0

	print star_name

	#####################################################################################################
	# 								REMOVE FILES FROM PREVIOUS RUNS 
	#####################################################################################################

	# Basically, remove anything that does not end in _cbcd_dn.fits or _cbcd.fits
	# But double check this

    #####################################################################################################
    # 									PERFORM PHOTOMETRY
    #####################################################################################################

	print "Aperture photometry"

    # Iterate over every epoch for the current star
	for epoch in range(1, num_epochs+1):

    	# Convert epoch into string version
		if epoch < 10:
			epoch_number = '0' + str(epoch)
		else: epoch_number = str(epoch)

		# Change directory to where this epoch is
		cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/'
		os.chdir(cwd)  
		
		# Aperture photometry on all 10 dithers for this epoch
		aper_phot(star_name, galaxy, channel, wavelength, epoch_number)


    #####################################################################################################
    # 						CREATE MEDIANED IMAGES FOR ON AND OFF TARGET FIELDS
    ##################################################################################################### 

	# Then need to come out of the epoch folders and just inside the channel folder so that medianed
	# image can be made and then copied to the epoch folders
	# Make medianed image, master star list and master PSF model for on-target field i.e. field 2

	print "On-target medianed image"

	os.chdir(os.path.expanduser('../')) # move up a level into channel folder
	master_on_target(star_name, galaxy, channel, wavelength, epoch_number, num_epochs)

	print "Off-target medianed image"

	# Then move back into the epoch folders again
	for epoch in range(1, num_epochs+1):

	    # Convert epoch into string version
		if epoch < 10:
			epoch_number = '0' + str(epoch)
		else: epoch_number = str(epoch)

		# Change directory to where this epoch is
		cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/'
		os.chdir(cwd) 

		# Make medianed image, master star list and master PSF model for off-target field i.e. field 1
		master_off_target(star_name, galaxy, channel, wavelength, epoch_number, num_epochs)


    #####################################################################################################
    # 								PSF PHOTOMETRY USING ALLSTAR
    ##################################################################################################### 

	# PSF photometry on all 10 dithers for this epoch
	print "PSF photometry"

    # Iterate over every epoch for the current star
	for epoch in range(1, num_epochs+1):

    	# Convert epoch into string version
		if epoch < 10:
			epoch_number = '0' + str(epoch)
		else: epoch_number = str(epoch)

		# Change directory to where this epoch is
		cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/'
		os.chdir(cwd)  
		
		# PSF photometry on all 10 dithers for this epoch
		psf_phot(star_name, galaxy, channel, wavelength, epoch_number)


    #####################################################################################################
    # 										ALLFRAME
    ##################################################################################################### 

	# Run ALLFRAME
	# Files created are .alf files
	print "ALLFRAME"

	for epoch in range(1, num_epochs+1):

    	# Convert epoch into string version
		if epoch < 10:
			epoch_number = '0' + str(epoch)
		else: epoch_number = str(epoch)

		# Change directory to where this epoch is
		cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/'
		os.chdir(cwd)  

		# Run ALLFRAME on both fields consisting of 5 dithers each
		for field in [1,2]:

			field = str(field)
			allframe(star_name, galaxy, channel, wavelength, epoch_number, field)


    #####################################################################################################
    # 							CALIBRATION TO STANDARD IRAC VEGA SYSTEM
    ##################################################################################################### 


	# Carry out calibration procedure to get magnitudes onto the standard IRAC Vega system of Reach et al. (2005)
	# Steps:
	# 1. Aperture correction 
	# 2. Standard aperture correction
	# 3. Zero point correction
	# 4. Location correction
	# 5. Pixel phase correction 
	# Files created are .alf_cal files

	print "Calibration to standard system"


	for epoch in range(1, num_epochs+1):

    	# Convert epoch into string version
		if epoch < 10:
			epoch_number = '0' + str(epoch)
		else: epoch_number = str(epoch)

		# Change directory to where this epoch is
		cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/'
		os.chdir(cwd)  

		# Run ALLFRAME on both fields consisting of 5 dithers each
		for field in [1,2]:

			field = str(field)
			calibration_procedure(star_name, galaxy, channel, wavelength, epoch_number, field)


    #####################################################################################################
    # 								OBTAIN FINAL MAGNITUDES FILE
    ##################################################################################################### 

	# Obtain final magnitude files that have combined the 5 dithers for each field into one file
	# Files created are .alf_all
	print "Create combined magnitude files for each epoch"

	for epoch in range(1, num_epochs+1):

    	# Convert epoch into string version
		if epoch < 10:
			epoch_number = '0' + str(epoch)
		else: epoch_number = str(epoch)

		# Change directory to where this epoch is
		cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/'
		os.chdir(cwd)  

		# Run ALLFRAME on both fields consisting of 5 dithers each
		for field in [1,2]:

			field = str(field)
			combine_dithers(star_name, galaxy, channel, wavelength, epoch_number, field)


    #####################################################################################################
    # 								GET AVERAGE MAGNITUDE AT EACH EPOCH
    ##################################################################################################### 

	# Obtain an average magnitude for each star at each epoch using the .cal file
	print "Averaging magnitudes for each epoch"

	for epoch in range(1, num_epochs+1):

    	# Convert epoch into string version
		if epoch < 10:
			epoch_number = '0' + str(epoch)
		else: epoch_number = str(epoch)

		# Change directory to where this epoch is
		cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/'
		os.chdir(cwd)  

		# Run ALLFRAME on both fields consisting of 5 dithers each
		for field in [1,2]:

			field = str(field)
			ave_mag(star_name, galaxy, channel, wavelength, epoch_number, field)


    #####################################################################################################
    # 							FIND MAGNITUDE OF CEPHEID AT EACH EPOCH
    ##################################################################################################### 

	# Now for this Cepheid, obtain its magnitude at each epoch and write out to a file
	# If ch1, want f2. If ch2, want f1.
	print "Get magnitudes at each epoch for star " + star_name
	get_mag(star_name, galaxy, channel, wavelength, num_epochs, ra, dec)


	#####################################################################################################
	# 										RUN GLOESS
	#####################################################################################################


	# Format magnitude file into a way in which GLOESS code accpets
	print "Formating file to put through GLOESS"
	format_gloess(star_name, galaxy, channel, wavelength, num_epochs, period)

	# Run GLOESS fitting code to obtain mean magnitude and light curve for Cepheid
	print "Fitting GLOESS curve"
	gloess_single_band(star_name, galaxy, channel, wavelength)

	print '------------------------------------------------------'

	#####################################################################################################
	# 								REPEAT ON ALL STARS IN FILE 
	#####################################################################################################

end = time.time()
print(end - start)