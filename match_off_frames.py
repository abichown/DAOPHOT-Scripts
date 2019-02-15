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

from functions import aper_phot, master_on_target, master_off_target, psf_phot, allframe, calibration_procedure, combine_dithers, ave_mag, get_mag, format_gloess, gloess_single_band
from functions import format_cal_files, match_all_frames, average_mags, find_target, gloess_files, gloess_multiband, master_off_target_no_psf, master_off_target_update

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

			# Make medianed image, master star list and master PSF model for off-target field 
			# master_off_target(star_name, galaxy, channel, wavelength, epoch_number, num_epochs)
			master_off_target_update(star_name, galaxy, channel, wavelength, epoch_number, num_epochs)


	print '------------------------------------------------------'

	#####################################################################################################
	# 								REPEAT ON ALL STARS IN FILE 
	#####################################################################################################

end = time.time()
print(end - start)