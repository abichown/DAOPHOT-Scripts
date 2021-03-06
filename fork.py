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

#import matplotlib.pyplot as matplotlib.pyplot
# from kneed import DataGenerator, KneeLocator

import matplotlib 
import gloess_fits as gf
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
matplotlib.rc('text',usetex=True)
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Garamond']

from new_functions import aper_phot, master_on_target, master_off_target, psf_phot, allframe, calibration_procedure, combine_dithers, ave_mag, get_mag, format_gloess, gloess_single_band
from new_functions import format_cal_files, match_all_frames, average_mags, find_target, gloess_files, gloess_multiband, master_off_target_no_psf, master_on_target_both_channels, allframe_new, remove_neighbours, calc_distance, calibration_procedure_testing, master_on_target_no_timeout
from new_functions import aperture_correction, calibration_procedure_loc_px, calibration_procedure_10_12_20, aper_phot_10_12_20
from new_functions import calculate_aperture_correction, apply_aperture_correction, calc_ave_epoch_offset, apply_epoch_offset, loc_px_corrections, entire_cal_procedure, delete_intermediate_files, initial_setup, master_image
from new_functions import master_image_test, average_mags_jackknife, average_mags_weighted
from new_functions import master_on_target_new

start = time.time()

#####################################################################################################
# 										READ IN STAR LIST
#####################################################################################################


# Set up input file into a dataframe
#df = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star','Period', 'RA', 'Dec', 'Channel'])
df = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star','Period', 'RA', 'Dec'])
print df

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
	print "Removing files from previous runs"
	initial_setup(star_name, galaxy)


    #####################################################################################################
    # 									PERFORM PHOTOMETRY
    #####################################################################################################

	print "Aperture photometry"

	for channel in [1,2]:

		if channel == 1:
			channel = '1'
			wavelength = '3p6um'
		elif channel == 2:
			channel = '2'
			wavelength = '4p5um'
		else: wavelength = 'channel not defined'

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
    # 				CREATE MEDIANED IMAGE FOR ON TARGET FIELDS FOR BOTH CHANNELS AT ONCE
    ##################################################################################################### 

	print "On-target medianed image from [3.6] AND [4.5] frames"

	os.chdir(os.path.expanduser('../')) # move up a level into channel folder
	#master_on_target_both_channels(star_name, galaxy, num_epochs)
	#master_on_target_no_timeout(star_name, galaxy, num_epochs)
	#master_on_target_new(star_name, galaxy, num_epochs)
	#master_image(star_name, galaxy, num_epochs)
	master_image_test(star_name, galaxy, num_epochs, ra, dec)

    #####################################################################################################
    # 						CREATE MEDIANED IMAGE FOR OFF TARGET FIELDS 
    #####################################################################################################


	for channel in ['1','2']:

		print "Off-target medianed image for channel " + channel

		if channel == '1':
			wavelength = '3p6um'
		else: wavelength = '4p5um'

		# Then move back into the epoch folders again
		for epoch in range(1, num_epochs+1):

		    # Convert epoch into string version
			if epoch < 10:
				epoch_number = '0' + str(epoch)
			else: epoch_number = str(epoch)

			# Change directory to where this epoch is
			cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/'
			os.chdir(cwd) 

			# Make medianed image and star lists for off target fields
			master_off_target_no_psf(star_name, galaxy, channel, wavelength, epoch_number, num_epochs)


    #####################################################################################################
    # 								PSF PHOTOMETRY USING ALLSTAR
    ##################################################################################################### 

	# PSF photometry on all 10 dithers for this epoch
	print "PSF photometry"

	for channel in ['1','2']:

		if channel == '1':
			wavelength = '3p6um'
		else: wavelength = '4p5um'

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

	for channel in ['1','2']:

		if channel == '1':
			wavelength = '3p6um'
		else: wavelength = '4p5um'

		for epoch in range(1, num_epochs+1):

	    	# Convert epoch into string version
			if epoch < 10:
				epoch_number = '0' + str(epoch)
			else: epoch_number = str(epoch)

			# Change directory to where this epoch is
			cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/'
			os.chdir(cwd)  

			# Run ALLFRAME on both fields consisting of 5 dithers each
			# 3p6 - field 1 is off, field 2 is on
			# 4p5 - field 1 is on, field 2 is off
			for field in [1,2]:

				field = str(field)
				allframe_new(star_name, galaxy, channel, wavelength, epoch_number, field)


 #    #####################################################################################################
 #    # 								APERTURE PHOTOMETRY 10 12 20
 #    ##################################################################################################### 

	# for channel in [1,2]:

	# 	if channel == 1:
	# 		channel = '1'
	# 		wavelength = '3p6um'
	# 	elif channel == 2:
	# 		channel = '2'
	# 		wavelength = '4p5um'
	# 	else: wavelength = 'channel not defined'

	# 	print "Aperture photometry"

	#     # Iterate over every epoch for the current star
	# 	for epoch in range(1, num_epochs+1):

	#     	# Convert epoch into string version
	# 		if epoch < 10:
	# 			epoch_number = '0' + str(epoch)
	# 		else: epoch_number = str(epoch)

	# 		# Change directory to where this epoch is
	# 		cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/'
	# 		os.chdir(cwd)  
			
	# 		# Aperture photometry on all 10 dithers for this epoch
	# 	 	aper_phot_10_12_20(star_name, galaxy, channel, wavelength, epoch_number)

 #    #####################################################################################################
 #    # 			PREPROCESSING BEFORE CALIBRATION - REMOVE STARS WITH NEIGHBOURS IN ANNULUS
 #    ##################################################################################################### 

	# # Remove stars with neighbouring stars inside the annulus as this will mess up the sky flux calculation
	# print "Removing PSF stars with contaminating neighbours"

	# for channel in ['1','2']:

	# 	if channel == '1':
	# 		wavelength = '3p6um'
	# 	else: wavelength = '4p5um'

	# 	for epoch in range(1, num_epochs+1):

	#     	# Convert epoch into string version
	# 		if epoch < 10:
	# 			epoch_number = '0' + str(epoch)
	# 		else: epoch_number = str(epoch)

	# 		# Change directory to where this epoch is
	# 		cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/'
	# 		os.chdir(cwd)  

	# 		# Perform on all 10 dithers
	# 		for dither in range(1,11):

	# 			dither = str(dither)

	# 			remove_neighbours(star_name, galaxy, wavelength, epoch_number, dither)

    #####################################################################################################
    # 							NEW CALIBRATION TO STANDARD IRAC VEGA SYSTEM
    ##################################################################################################### 

	print "New calibration "

	for channel in ['1','2']:

		if channel == '1':
			start_dither = 6
			wavelength = '3p6um'
		else:
			wavelength = '4p5um'
			start_dither = 1

		for epoch in range(1,num_epochs+1):

			for dither in range(start_dither, start_dither+5):

				entire_cal_procedure(star_name, galaxy, channel, wavelength, epoch, dither)



		# # Calculate aperture correction value
		# apc = calculate_aperture_correction(star_name, galaxy, channel, wavelength)

		# # Now apply aperture correction to the reference image
		# apply_aperture_correction(star_name, galaxy, channel, wavelength, apc)

		# # Calculate average offset per epoch compared to ref image
		# offsets = calc_ave_epoch_offset(star_name, galaxy, channel, wavelength)

		# # Apply these offsets
		# apply_epoch_offset(star_name, galaxy, channel, wavelength, offsets)

		# # Location and pixel phase correction
		# loc_px_corrections(star_name, galaxy, channel, wavelength)


 #    #####################################################################################################
 #    # 							CALIBRATION TO STANDARD IRAC VEGA SYSTEM
 #    ##################################################################################################### 


	# # Carry out calibration procedure to get magnitudes onto the standard IRAC Vega system of Reach et al. (2005)
	# # Steps:
	# # 1. Aperture correction 
	# # 2. Standard aperture correction
	# # 3. Zero point correction
	# # 4. Location correction
	# # 5. Pixel phase correction 
	# # Files created are .alf_cal files

	# print "Calibration to standard system"

	# for channel in ['1','2']:

	# 	if channel == '1':
	# 		wavelength = '3p6um'
	# 	else: wavelength = '4p5um'


	# 	aperture_correction(star_name, galaxy, channel, wavelength)


	# 	for epoch in range(1, num_epochs+1):

	#     	# Convert epoch into string version
	# 		if epoch < 10:
	# 			epoch_number = '0' + str(epoch)
	# 		else: epoch_number = str(epoch)

	# 		# Change directory to where this epoch is
	# 		cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/'
	# 		os.chdir(cwd)  

	# 		# Run ALLFRAME on both fields consisting of 5 dithers each
	# 		for field in [1,2]:

	# 			field = str(field)
	# 			#calibration_procedure(star_name, galaxy, channel, wavelength, epoch_number, field)
	# 			#calibration_procedure_testing(star_name, galaxy, channel, wavelength, epoch_number, field)
	# 			calibration_procedure_loc_px(star_name, galaxy, channel, wavelength, epoch_number, field)


    #####################################################################################################
    # 							FORMAT .ALF_CAL FILES TO .ALF_ALL FILES
    ##################################################################################################### 

	# Format the magnitude files so that they are in the correct format for DAOMATCH
	# Files created are .alf_all
	print "Formatting .alf_cal files to .alf_all files"

	for channel in ['1','2']:

		if channel == '1':
			wavelength = '3p6um'
			field = '2'
		else: 
			wavelength = '4p5um'
			field = '1'	

		for epoch in range(1, num_epochs+1):

	    	# Convert epoch into string version
			if epoch < 10:
				epoch_number = '0' + str(epoch)
			else: epoch_number = str(epoch)

			# Change directory to where this epoch is
			cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/'
			os.chdir(cwd)  

			# # Do on both fields
			# for field in [1,2]:

			# 	field = str(field)
			format_cal_files(star_name, galaxy, channel, wavelength, epoch_number, field)


	##################################################################################################
	# 									MATCH ALL ON TARGET FRAMES	
	##################################################################################################

	# Want to match all on target aperture magnitude files to create one giant file containing the mags
	# for all stars across all dithers, all epochs and both channels

	print "Making giant file of magnitudes"
	match_all_frames(star_name, galaxy, num_epochs)

	#################################################################################################
	#									GET AVERAGE MAGNITUDES
	#################################################################################################

	# In the giant magnitudes file is now all dithers from all epochs from both channels
	# Want to create a new dataframe which only includes the following
	# ID, X, Y, ave_mag at each epoch, ave_err at each epoch

	print "Averaging magnitudes"

	# Move up one level to star folder
	os.chdir(os.path.expanduser('../'))

	# This dataframe has the following columns
	# ID X Y 3p6_e01_ave_mag 3p6_e01_ave_err ... 4p5_e24_ave_mag 4p5_e24_ave_err 
	# Only stars with ave_mags at ALL epochs are kept
	# This dataframe will be used to search for the Cepheid
	#ave_df = average_mags(star_name, galaxy, num_epochs)
	#ave_df = average_mags_jackknife(star_name, galaxy, num_epochs)
	ave_df = average_mags_weighted(star_name, galaxy, num_epochs)


	##################################################################################################
	# 										FIND CEPHEID
	##################################################################################################

	# Then look for star by comparing coordinates between pixel and NED
	# Also if num epochs is not 24 (for LMC) then calculate variability index

	print "Finding Cepheid"
	find_target(star_name, galaxy, ra, dec, num_epochs, ave_df)


	##################################################################################################
	# 								FORMAT TO GLOESS FILE
	##################################################################################################

	print "Formatting to GLOESS format"
	gloess_files(star_name, galaxy, num_epochs, period)

	##################################################################################################
	# 									FIT GLOESS CURVE
	##################################################################################################

	# Then fit GLOESS as usual
	# Will update for it to do both [3.6], [4.5] and colour curve

	print "Fitting GLOESS light curves and colour curve"
	gloess_multiband(star_name, galaxy)


	##################################################################################################
	# 									TIDYING UP DIRECTORIES
	##################################################################################################

	print "Tidying up files"
	delete_intermediate_files(star_name, galaxy)

	shutil.rmtree('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/temp/')
	# shutil.rmtree('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch1/apc_temp/')
	# shutil.rmtree('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch2/apc_temp/')

	print '------------------------------------------------------'

	#####################################################################################################
	# 								REPEAT ON ALL STARS IN FILE 
	#####################################################################################################