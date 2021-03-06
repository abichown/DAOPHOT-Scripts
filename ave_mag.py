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
import numpy as np
from math import sqrt, log10

# Set up data frame from txt file of stars (sys.argv[1]) to do it on
df = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star', 'Period', 'RA', 'Dec', 'Channel'])

# Loop over each row in txt file
# Also loop over the two dither combinations
for i in range(0, len(df)):
	for j in [6]: #[1,6] when field 1 working 

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

		if j == 1:
			field = '1'
		elif j == 6:
			field = '2'
		else: field = 'Invalid start dither'

		# Get number of epochs
		if galaxy == 'LMC':
				num_epochs = 24
		elif galaxy == 'SMC':
				num_epochs = 12
		else: num_epochs = 0

			# Run over each epoch
		for epoch in range(1, num_epochs+1):

			if epoch < 10:
				epoch_number = '0' + str(epoch)
			else: epoch_number = str(epoch)

	    	# Find absolute path of where images are
			home = '/home/ac833/Data/'

			cwd = home + str(galaxy) + '/BCD/' + target_name + '/ch' + str(df['Channel'][i]) + '/e' + str(epoch_number) + '/'
			stem = target_name + '_' + wavelength 

			# Change directory to where image is
			os.chdir(cwd)

			# Remove any previous runs of this particular script
			if (os.path.isfile(stem+'_f'+field+'.ave')):
				os.remove(stem+'_f'+field+'.ave')

			# Import data
			filename = target_name + '_' + wavelength + '_f' + field + '.cal'

			data = pd.read_csv(filename, header=0, delim_whitespace=True)	

			print "Working on: " + target_name + '    epoch: ' + epoch_number + '    field: ' + field 


			# Get right f0 for the channel
			if wavelength == '3p6um':
				f0 = 280.9
			if wavelength == '4p5um':
				f0 = 179.7		

			# New columns for fluxes
			data['F1'] = 0
			data['F2'] = 0
			data['F3'] = 0
			data['F4'] = 0
			data['F5'] = 0

			# Create new columns in dataframe for the average magnitude, average error and standard error on the mean for each star
			data['f_ave'] = 0
			data['m_ave'] = 0
			data['e_ave'] = 0
			data['std_err'] = 0	

			# Loop over all rows of dataframe to convert any non 99.9999 mags to a flux
			for q in range(0, len(data)):
				if data['M1'][q] != 99.9999:
					f = f0 * (10 ** (-data['M1'][q]/2.5))
					data.loc[q, 'F1'] = f
				if data['M2'][q] != 99.9999:
					f = f0 * (10 ** (-data['M2'][q]/2.5))
					data.loc[q, 'F2'] = f
				if data['M3'][q] != 99.9999:
					f = f0 * (10 ** (-data['M3'][q]/2.5))
					data.loc[q, 'F3'] = f
				if data['M4'][q] != 99.9999:
					f = f0 * (10 ** (-data['M4'][q]/2.5))
					data.loc[q, 'F4'] = f
				if data['M5'][q] != 99.9999:
					f = f0 * (10 ** (-data['M5'][q]/2.5))
					data.loc[q, 'F5'] = f
		

			# Determine average flux and error
			for k in range(0, len(data)):
				total = 0
				error = 0
				count = 0

				if data['M1'][k] != 99.9999:
					total += data['F1'][k]
					error += data['E1'][k]
					count += 1
				if data['M2'][k] != 99.9999:
					total += data['F2'][k]
					error += data['E2'][k]
					count += 1
				if data['M3'][k] != 99.9999:
					total += data['F3'][k]
					error += data['E3'][k]
					count += 1
				if data['M4'][k] != 99.9999:
					total += data['F4'][k]
					error += data['E4'][k]
					count += 1
				if data['M5'][k] != 99.9999:
					total += data['F5'][k]
					error += data['E5'][k]
					count += 1

				if count != 0:

					data.loc[k, 'f_ave'] = total/count
					data.loc[k, 'e_ave'] = round(error/count, 4)

			# Convert average flux to an average magnitude and put in m_ave column
			for q in range(0, len(data)):

				if data['f_ave'][q] != 0:

					m = -2.5 * log10(data['f_ave'][q]/f0)
					data.loc[q, 'm_ave'] = round(m, 4)


			# Calculate standard error on the mean for each star
			for k in range(0, len(data)):

				mags = []

				if data['M1'][k] != 99.9999:
					mags.append(data['M1'][k])
				if data['M2'][k] != 99.9999:
					mags.append(data['M2'][k])
				if data['M3'][k] != 99.9999:
					mags.append(data['M3'][k])
				if data['M4'][k] != 99.9999:
					mags.append(data['M4'][k])
				if data['M5'][k] != 99.9999:
					mags.append(data['M5'][k])

				data.loc[k, 'std_err'] = round(np.std(mags)/sqrt(len(mags)),4)

			# Write average values to a file
			filename = filename.replace('.cal', '.ave')
			data.to_csv(filename, columns=['ID', 'X', 'Y', 'm_ave', 'e_ave', 'std_err'], header=True, sep=" ", index=False)



