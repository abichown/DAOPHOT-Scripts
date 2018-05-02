'''
Purpose: Correct the instrumental magnitudes via zero point, standard aperture, aperture and location corrections
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

import numpy as np
import pandas as pd 
from astropy.io import fits
from math import log10 
import os

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
				MAIN PART
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Read in list of stars that need correcting
df = pd.read_csv('/home/ac833/DAOPHOT-Scripts/star_list.txt', header=None, delim_whitespace=True, names=['Galaxy', 'Star', 'Channel', 'Epoch'])

for i in range(0, len(df)):
	for j in [1,6]: 

		# Initial setup
		galaxy = df['Galaxy'][i]
		target_name = df['Star'][i]

		if df['Channel'] == 1:
			channel = 1
			wavelength = '3p6um'
		elif df['Channel'] == 2:
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
		else field = "Invalid field"

		# Find absolute path of where images are
		home = '/home/ac833/Data/'
		cwd = home + str(galaxy) + '/BCD/' + target_name + '/ch' + str(df['Channel'][i]) + '/e' + str(epoch_number) + '/'

		stem = target_name + '_' + wavelength + '_e' + epoch_number

		# Change working directory to where the data is
		os.chdir(cwd)

		# Zero point
		zero_point(i,j)

		# Standard aperture
		std_aper(i,j)

		# Aperture correction
		ap_corr(i,j)

		# Location correction
		loc_corr(i,j)

# Zero point correction - applied to both ap and alf mags
def zero_point():
	return(0)

# Get onto standard aperture system - applied to ap mags only
def std_aper():
	return(0)

# Aperture correction - uses ap and psf mags to calculate correction value then applies to psf mags only
def ap_corr():
	return(0)

# Location correction - only need to apply to psf mags as won't be using ap mags after this
def loc_corr():
	return(0)
