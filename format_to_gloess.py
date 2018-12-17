'''
Purpose: The magnitude file in the Magnitudes directory isn't in the appropriate format for GLOESS. This
script formats it correctly for a single channel
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

import os
import sys
import pandas as pd 
from astropy.io import fits

# Read star info
df = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star', 'Period', 'RA', 'Dec', 'Channel'])

for i in range(0, len(df)):

	# Get all the relevant star info
	galaxy = df['Galaxy'][i]
	target_star = df['Star'][i]
	channel = df['Channel'][i]
	period = df['Period'][i]

	# Change channel number to wavelength
	if channel == 1:
		wavelength = '3p6um'
	elif channel == 2:
		wavelength = '4p5um'
	else: wavelength = 'not_defined'

	# Get number of observations based on whether it is LMC or SMC
	if galaxy == 'LMC':
		n_obs = 24
	else: n_obs = 12

	# Get HMJD list for the observations
	hmjd = []

	for j in range(1,n_obs+1):

		if j < 10:
			filename = '/home/ac833/Data/'+str(galaxy)+'/BCD/'+str(target_star)+'/ch'+str(channel)+'/e0'+str(j)+'/'+str(target_star)+'_'+str(wavelength)+'_e0'+str(j)+"_d1_cbcd_dn.fits"
		else: filename = '/home/ac833/Data/'+str(galaxy)+'/BCD/'+str(target_star)+'/ch'+str(channel)+'/e'+str(j)+'/'+str(target_star)+'_'+str(wavelength)+'_e'+str(j)+"_d1_cbcd_dn.fits"

		image = fits.open(filename)
		header = image[0].header

		hmjd.append(header['HMJD_OBS'])

	# Get magnitudes and errors
	mag_file = '/home/ac833/Magnitudes/'+str(galaxy)+'/'+str(target_star)+'_'+wavelength+'.txt'

	mags = pd.read_csv(mag_file, header=0, delim_whitespace=True)

	magnitudes = list(mags['Mag'])
	errors = list(mags['Std_err'])

	# Write out file in correct format
	output_filename = '/home/ac833/GLOESS_files/'+str(target_star)+'_'+wavelength+'_gloess.txt'

	f = open(output_filename, 'w')
	f.write(str(target_star) + '\n')
	f.write(str(period) + '\n')
	f.write(str(1) + '\n') # doesn't matter what goes here as it's not used by GLOESS code

	if period <= 12:
		f.write(str(0.25) + '\n') 
	else: f.write(str(0.2) + '\n')

	for x in range(0, n_obs):
		f.write(str(hmjd[x]) + ' ' + str(magnitudes[x]) + ' ' + str(errors[x]) + '\n')

	f.close()


