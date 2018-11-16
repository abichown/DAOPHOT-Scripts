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
df = pd.read_csv('/home/ac833/DAOPHOT-Scripts/star_list.txt', header=None, delim_whitespace=True, names=['Galaxy', 'Star', 'Period', 'Channel'])

# Create list of unique star entries in df
#stars = pd.unique(df['Star'])

for i in range(0, len(df)):

	# Find rows that match this star
	#matches = df.index[(df['Star']==star)==True].tolist()

	# Just need the first occurrence to find Galaxy
	#index = matches[0]

	# Then get galxy name
	galaxy = df['Galaxy'][i]
	#galaxy = df['Galaxy'][index]

	# Get star name
	star = df['Star'][i]

	# Wavelengths list
	#wavelengths = ['3p6um', '4p5um']

    # Get channel and convert to wavelength
	if df['Channel'][i] == 1:
		channel = '1'
		wavelength = '3p6um'
		field = '2' # we want the field with the V* in it
	elif df['Channel'][i] == 2:
		channel = '2'
		wavelength = '4p5um'
		field = '1' # we want the field with the V* in it
	else: wavelength = 'channel not defined'

	# Number of epochs
	if galaxy == 'LMC':
		num_epochs = 24
	elif galaxy == 'SMC':
		num_epochs = 12
	else: num_epochs = 0

	#for wavelength in wavelengths:

		# # Get channel number
		# if wavelength == '3p6um':
		# 	channel = '1'
		# 	field = '2' # we want the field with the V* in
		# if wavelength == '4p5um':
		# 	channel = '2'
		# 	field = '1' # we want the field with the V* in

	# Open file to write mag and error to 
	filename = '/home/ac833/Magnitudes/' + galaxy + '/' + star + '_' + wavelength +'.txt'
	f = open(filename, 'w')
	f.write("Epoch Mag Error Std_err \n")

	for epoch in range(1,num_epochs+1):

		# Get epoch in correct form
		if epoch < 10:
			epoch = '0' + str(epoch)
		else: epoch = str(epoch)

		# Change cwd to folder with data in
		cwd = '/home/ac833/Data/' + galaxy + '/BCD/' + star + '/ch' + channel + '/e' + epoch + '/'
		os.chdir(cwd)

		# Open file containing average magnitudes
		ave_file = star + '_' + wavelength + '_f' + field + '.ave'

		if (os.path.isfile(ave_file)):
			ave = pd.read_csv(ave_file, header=0, delim_whitespace=True)

			# Search for star near (133,122)
			# I've played with these ranges to ensure HV00872 ch1 and ch2 data get included
			# Might need to be changed
			for i in range(0, len(ave)):
				if (ave['X'][i] < 135.00 and ave['X'][i] > 128.00):
					if (ave['Y'][i] < 125.00 and ave['Y'][i] > 119.00):
						# write to file
						f.writelines("%s %f %f %f \n" % (epoch, float(ave['m_ave'][i]), float(ave['e_ave'][i]), float(ave['std_err'][i])))

	f.close()
