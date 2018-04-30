'''
Purpose: Extract magnitude and error of the V*. At the moment, V* always in centre of frame
so program looks for star near (133,122). Also, V* only present in [3.6] field 2 and 
[4.5] field 1 so only looks there.
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

# Import modules
import os
import sys
import pandas as pd
import numpy as np

# Read star list
df = pd.read_csv('/home/ac833/DAOPHOT-Scripts/star_list.txt', header=None, delim_whitespace=True, names=['Galaxy', 'Star', 'Channel', 'Epoch'])

# Create list of unique star entries in df
stars = pd.unique(df['Star'])

for star in stars:

	# Find rows that match this star
	matches = df.index[(df['Star']==star)==True].tolist()

	# Just need the first occurrence to find Galaxy
	index = matches[0]

	# Then get galxy name
	galaxy = df['Galaxy'][index]

	# Wavelengths list
	wavelengths = ['3p6um', '4p5um']

	for wavelength in wavelengths:

		# Get channel number
		if wavelength == '3p6um':
			channel = '1'
			field = '2' # we want the field with the V* in
		if wavelength == '4p5um':
			channel = '2'
			field = '1' # we want the field with the V* in

		

	print galaxy