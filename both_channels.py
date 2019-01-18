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

from functions import aper_phot, master_on_target, master_off_target, psf_phot, allframe, calibration_procedure, combine_dithers, ave_mag, get_mag, format_gloess, gloess_single_band, format_cal_files, match_all_frames, average_mags, find_target, gloess_files, gloess_multiband

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

	for channel in [1,2]:

		if channel == 1:
			channel = '1'
			wavelength = '3p6um'
		elif channel == 2:
			channel = '2'
			wavelength = '4p5um'
		else: wavelength = 'channel not defined'

	# 	#####################################################################################################
	# 	# 								REMOVE FILES FROM PREVIOUS RUNS 
	# 	#####################################################################################################

	# 	# Basically, remove anything that does not end in _cbcd_dn.fits or _cbcd.fits
	# 	# But double check this

	#     #####################################################################################################
	#     # 									PERFORM PHOTOMETRY
	#     #####################################################################################################

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
	# 		aper_phot(star_name, galaxy, channel, wavelength, epoch_number)


	#     #####################################################################################################
	#     # 						CREATE MEDIANED IMAGES FOR ON AND OFF TARGET FIELDS
	#     ##################################################################################################### 

	# 	# Then need to come out of the epoch folders and just inside the channel folder so that medianed
	# 	# image can be made and then copied to the epoch folders
	# 	# Make medianed image, master star list and master PSF model for on-target field i.e. field 2

	# 	print "On-target medianed image"

	# 	os.chdir(os.path.expanduser('../')) # move up a level into channel folder
	# 	master_on_target(star_name, galaxy, channel, wavelength, epoch_number, num_epochs)

	# 	print "Off-target medianed image"

	# 	# Then move back into the epoch folders again
	# 	for epoch in range(1, num_epochs+1):

	# 	    # Convert epoch into string version
	# 		if epoch < 10:
	# 			epoch_number = '0' + str(epoch)
	# 		else: epoch_number = str(epoch)

	# 		# Change directory to where this epoch is
	# 		cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/'
	# 		os.chdir(cwd) 

	# 		# Make medianed image, master star list and master PSF model for off-target field 
	# 		master_off_target(star_name, galaxy, channel, wavelength, epoch_number, num_epochs)


	#     #####################################################################################################
	#     # 								PSF PHOTOMETRY USING ALLSTAR
	#     ##################################################################################################### 

	# 	# PSF photometry on all 10 dithers for this epoch
	# 	print "PSF photometry"

	#     # Iterate over every epoch for the current star
	# 	for epoch in range(1, num_epochs+1):

	#     	# Convert epoch into string version
	# 		if epoch < 10:
	# 			epoch_number = '0' + str(epoch)
	# 		else: epoch_number = str(epoch)

	# 		# Change directory to where this epoch is
	# 		cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/'
	# 		os.chdir(cwd)  
			
	# 		# PSF photometry on all 10 dithers for this epoch
	# 		psf_phot(star_name, galaxy, channel, wavelength, epoch_number)


	#     #####################################################################################################
	#     # 										ALLFRAME
	#     ##################################################################################################### 

	# 	# Run ALLFRAME
	# 	# Files created are .alf files
	# 	print "ALLFRAME"

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
	# 			allframe(star_name, galaxy, channel, wavelength, epoch_number, field)


	#     #####################################################################################################
	#     # 							CALIBRATION TO STANDARD IRAC VEGA SYSTEM
	#     ##################################################################################################### 


	# 	# Carry out calibration procedure to get magnitudes onto the standard IRAC Vega system of Reach et al. (2005)
	# 	# Steps:
	# 	# 1. Aperture correction 
	# 	# 2. Standard aperture correction
	# 	# 3. Zero point correction
	# 	# 4. Location correction
	# 	# 5. Pixel phase correction 
	# 	# Files created are .alf_cal files

	# 	print "Calibration to standard system"


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
	# 			calibration_procedure(star_name, galaxy, channel, wavelength, epoch_number, field)


	#     #####################################################################################################
	#     # 							FORMAT .ALF_CAL FILES TO .ALF_ALL FILES
	#     ##################################################################################################### 

	# 	# Format the magnitude files so that they are in the correct format for DAOMATCH
	# 	# Files created are .alf_all
	# 	print "Formatting .alf_cal files to .alf_all files"

	# 	for epoch in range(1, num_epochs+1):

	#     	# Convert epoch into string version
	# 		if epoch < 10:
	# 			epoch_number = '0' + str(epoch)
	# 		else: epoch_number = str(epoch)

	# 		# Change directory to where this epoch is
	# 		cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number +'/'
	# 		os.chdir(cwd)  

	# 		# Do on both fields
	# 		for field in [1,2]:

	# 			field = str(field)
	# 			format_cal_files(star_name, galaxy, channel, wavelength, epoch_number, field)


	# ##################################################################################################
	# # 									MATCH ALL ON TARGET FRAMES	
	# ##################################################################################################

	# # Want to match all on target aperture magnitude files to create one giant file containing the mags
	# # for all stars across all dithers, all epochs and both channels

	# print "Making giant file of magnitudes"
	# match_all_frames(star_name, galaxy, num_epochs)

	# #################################################################################################
	# #									GET AVERAGE MAGNITUDES
	# #################################################################################################

	# # In the giant magnitudes file is now all dithers from all epochs from both channels
	# # Want to create a new dataframe which only includes the following
	# # ID, X, Y, ave_mag at each epoch, ave_err at each epoch

	# print "Averaging magnitudes"

	# # Move up one level to star folder
	# os.chdir(os.path.expanduser('../'))

	# # This dataframe has the following columns
	# # ID X Y 3p6_e01_ave_mag 3p6_e01_ave_err ... 4p5_e24_ave_mag 4p5_e24_ave_err 
	# # Only stars with ave_mags at ALL epochs are kept
	# # This dataframe will be used to search for the Cepheid
	# df = average_mags(star_name, galaxy, num_epochs)

	# ##################################################################################################
	# # 										FIND CEPHEID
	# ##################################################################################################

	# # Then look for star by comparing coordinates between pixel and NED
	# # Also if num epochs is not 24 (for LMC) then calculate variability index

	# print "Finding Cepheid"
	# find_target(star_name, galaxy, ra, dec, num_epochs, df)


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


	print '------------------------------------------------------'

	#####################################################################################################
	# 								REPEAT ON ALL STARS IN FILE 
	#####################################################################################################

end = time.time()
print(end - start)