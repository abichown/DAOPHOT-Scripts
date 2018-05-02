'''
Purpose: Correct the instrumental magnitudes via zero point, standard aperture, aperture and location corrections
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

import numpy as np
import pandas as pd 
from astropy.io import fits
from math import log10 

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
				MAIN PART
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Read in list of stars that need correcting
df = pd.read_csv('/home/ac833/DAOPHOT-Scripts/star_list.txt', header=None, delim_whitespace=True, names=['Galaxy', 'Star', 'Channel', 'Epoch'])

for i in range(0, len(df)):
	for j in [1,6]: 

		# Initial setup

		# Zero point

		# Standard aperture

		# Aperture correction

		# Location correction

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
