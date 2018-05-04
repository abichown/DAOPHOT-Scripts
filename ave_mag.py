'''
Purpose: Determines the average magnitude of each star across the 5 dithers for each 
field (dependent on whether it was actually in the frame). Also determines average
error. Writes out these new values to a file with extension .ave
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

# Import modules
import sys
import os
import pandas as pd
from math import sqrt, log10

# Set up data frame from txt file of stars (sys.argv[1]) to do it on
df = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star','Channel','Epoch'])

# Loop over each row in txt file
# Also loop over the two dither combinations
for i in range(0, len(df)):
	for j in [1,6]:

		# INITIAL SETUP 

		galaxy = df['Galaxy'][i]
		target_name = df['Star'][i]

		if df['Channel'][i] == 1:
			channel = 1
			wavelength = '3p6um'
		elif df['Channel'][i] == 2:
			channel = 2
			wavelength = '4p5um'
		else: wavelength = 'channel not defined'

		if df['Epoch'][i] < 10:
			epoch_number = '0' + str(df['Epoch'][i])
		else: epoch_number = str(df['Epoch'][i])

		if j == 1:
			field = '1'
		elif j == 6:
			field = '2'
		else: field = 'Invalid start dither'

    	# Find absolute path of where images are
		home = '/home/ac833/Data/'

		cwd = home + str(galaxy) + '/BCD/' + target_name + '/ch' + str(df['Channel'][i]) + '/e' + str(epoch_number) + '/'
		stem = target_name + '_' + wavelength + '_e' + epoch_number

		# Change directory to where image is
		os.chdir(cwd)

		# Remove any previous runs of this particular script
		if (os.path.isfile(stem+'_f'+field+'.ave')):
			os.remove(stem+'_f'+field+'.ave')

		# Import data
		filename = stem + '_f' + field + '.alf_all'
		data = pd.read_csv(filename, delim_whitespace=True, header=0)	

		print "Working on: " + target_name + '    epoch: ' + epoch_number + '    field: ' + field 

		# Get right f0 for the channel
		if wavelength == '3p6um':
			f0 = 280.9
		if wavelength == '4p5um':
			f0 = 179.7


		# Loop over all rows of df to convert any non 99.9999 mags to a flux
		for q in range(0, len(data)):
			if data['M1'][q] != 99.9999:
				f = f0 * (10 ** (-data['M1'][q]/2.5))
				data.loc[q, 'M1'] = f
			if data['M2'][q] != 99.9999:
				f = f0 * (10 ** (-data['M2'][q]/2.5))
				data.loc[q, 'M2'] = f
			if data['M3'][q] != 99.9999:
				f = f0 * (10 ** (-data['M3'][q]/2.5))
				data.loc[q, 'M3'] = f
			if data['M4'][q] != 99.9999:
				f = f0 * (10 ** (-data['M4'][q]/2.5))
				data.loc[q, 'M4'] = f
			if data['M5'][q] != 99.9999:
				f = f0 * (10 ** (-data['M5'][q]/2.5))
				data.loc[q, 'M5'] = f

		# Create new columns for average magnitude and errors
		data['m_ave'] = 0
		data['e_ave'] = 0

		# Determine average flux and error
		for k in range(0, len(data)):
			total = 0
			error = 0
			count = 0

			if data['M1'][k] != 99.9999:
				total += data['M1'][k]
				error += data['E1'][k]
				count += 1
			if data['M2'][k] != 99.9999:
				total += data['M2'][k]
				error += data['E2'][k]
				count += 1
			if data['M3'][k] != 99.9999:
				total += data['M3'][k]
				error += data['E3'][k]
				count += 1
			if data['M4'][k] != 99.9999:
				total += data['M4'][k]
				error += data['E4'][k]
				count += 1
			if data['M5'][k] != 99.9999:
				total += data['M5'][k]
				error += data['E5'][k]
				count += 1

			data.loc[k, 'm_ave'] = total/count
			data.loc[k, 'e_ave'] = error/count

		# Convert back to magnitudes
		for p in range(0, len(data)):
			m = -2.5 * log10(data['m_ave'][p]/f0)
			data.loc[p, 'm_ave'] = m

		# Finally, write new values to a file
		filename = filename.replace('.alf_all', '.ave')
		data.to_csv(filename, columns=['ID', 'X', 'Y', 'm_ave', 'e_ave'], header=True, sep=" ", index=None)


        print "Complete"

