'''
Purpose: Correct the instrumental magnitudes via zero point, standard aperture, aperture and location corrections
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

import numpy as np
import pandas as pd 
from astropy.io import fits
from math import log10 
import os
import sys

# Zero point correction - applied to both ap and alf mags
def zero_point(row_of_df, start_dither):

	# Image name - only need to do it once as they are all using same telescope etc
	image = stem + '_d1_cbcd_dn.fits'

	# Calculate zmag
	fits_image = fits.open(image)
	hdr = fits_image[0].header
	fluxconv = hdr['fluxconv'] * (10 ** 6) 
	px_ste = 3.3845 * (10 ** -11) # px size in steradians using online converters

	if channel == 1:
		F0 = 280.9
	elif channel == 2:
		F0 = 179.7
	else: F0 = 'Invalid'

	zmag = round(2.5 * log10(F0/(fluxconv*px_ste)), 2)

	# Load files


	return(0)

# Get onto standard aperture system - applied to ap mags only
def std_aper():
	return(0)

# Aperture correction - uses ap and alf mags to calculate correction value then applies to psf mags only
def ap_corr():
	return(0)

# Location correction - only need to apply to alf mags as won't be using ap mags after this
def loc_corr():
	return(0)

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
				MAIN PART
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Read in list of stars that need correcting
df = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star', 'Channel', 'Epoch'])

for i in range(0, len(df)):
	for j in [1,6]: 

		# Initial setup
		galaxy = df['Galaxy'][i]
		target_name = df['Star'][i]

		if df['Channel'][i] == 1:
			channel = 1
			wavelength = '3p6um'
		elif df['Channel'][i] == 2:
			channel = 2
			wavelength = '4p5um'
		else: channel = 'Invalid channel'

		if df['Epoch'][i] < 10:
			epoch_number = '0' + str(df['Epoch'][i])
		else: epoch_number = str(df['Epoch'][i])

		if j == 1:
			field = 1
		elif j == 6:
			field = 2
		else: field = "Invalid field"

		# Find absolute path of where images are
		home = '/home/ac833/Data/'
		cwd = home + str(galaxy) + '/BCD/' + target_name + '/ch' + str(df['Channel'][i]) + '/e' + str(epoch_number) + '/'

		stem = target_name + '_' + wavelength + '_e' + epoch_number

		# Change working directory to where the data is
		os.chdir(cwd)

		# Set up names of ap and alf files - have 5 alf and 5 ap files for each epoch
		# Put them in two lists to make it easier to loop over later
		ap_files = []
		alf_files = []

		for dither in range(j, j+5):
			ap_files.append(stem + '_d' + str(dither) + '_cbcd_dn.ap')
			alf_files.append(stem + '_d' + str(dither) + '_cbcd_dn.alf')

		# Zero point
		zero_point(i,j)

		# Standard aperture
		#std_aper(i,j)

		# Aperture correction
		#ap_corr(i,j)

		# Location correction
		#loc_corr(i,j)
