'''
Purpose: Take the calibrated magnitudes (.alf_cal) from calibration_steps.py and put them into a file suitable for DAOMASTER.
Then run them through DAOMASTER to get file needed for averaging magnitudes in next step
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

import pandas as pd 
import os
import sys
import pexpect
import shutil

def format_alf(row_of_df, start_dither):

	# Format all 5 dithers for this star, epoch and field
	for i in range(start_dither, start_dither+5):

		# Get filenames
		alf = stem + '_d' + str(i) + '_cbcd_dn.alf' # has everything but the correct magnitudes
		alf_apc = stem + '_d' + str(i) + '_cbcd_dn.alf_cal' # has the magnitudes we want

		f = open(alf)

		# Write the 3 line header to a new file
		new_filename = stem + '_d' + str(i) + '_cbcd_dn.alf_all'
		header = f.read().splitlines()[0:3]
		g = open(new_filename, 'w')
		g.writelines(header[0] + '\n' + header[1] + '\n' + header[2] + '\n')

		# Open alf and alf_apc
		df = pd.read_csv(alf, delim_whitespace=True, header=None, skiprows=3, names=['ID', 'X', 'Y', 'Uncorr_Mag', 'Error', 'Sky', 'Iters', 'Chi', 'Sharp'])
		df2 = pd.read_csv(alf_apc, skiprows=1, delim_whitespace=True, header=None, names=['ID_nn', 'X_nn', 'Y_nn', 'mag', 'error'])

		# Concat the two df's and drop unnecessary columns
		data = pd.concat((df,df2), axis=1)
		data.drop(['Uncorr_Mag', 'ID_nn','X_nn', 'Y_nn', 'error'], axis=1, inplace=True)
		data = data[['ID', 'X', 'Y', 'mag', 'Error', 'Sky', 'Iters', 'Chi', 'Sharp']]

		g.close()
		f.close()

		# Append to the new file with the header at the top
		data.to_csv(new_filename, mode='a', header=None, sep=' ', index=False)

	return(0)


def daomaster(row_of_df, start_dither):

	# File with the coordinate transformations in
	match_file = target_name + '_' + wavelength + '_f' + str(field) + '_master.mch'

	# Change the file extensions in the .mch file from .alf to .alf_all
	f = open(match_file, 'r')
	filedata = f.read()
	f.close()

	newdata = filedata.replace(".alf",".alf_all")

	f = open(match_file,'w')
	f.write(newdata)
	f.close()

	# Copy epoch 1 dither 1/6 .alf_all file
	if epoch_number != '01':
		shutil.copy(home + str(galaxy) + '/BCD/' + target_name + '/ch' + str(df['Channel'][i]) + '/e01/' + target_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.alf_all', target_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.alf_all')

	# Run DAOMASTER - use coordinate transformations from previous steps
	daomaster = pexpect.spawn('daomaster')

	daomaster.expect("File with list of input files:")
	daomaster.sendline(match_file)
	daomaster.expect("Minimum number, minimum fraction, enough frames:")
	daomaster.sendline("1, 0.5, 5") # play around with these values
	daomaster.expect("Maximum sigma:")
	daomaster.sendline("0.5") # play around with this value
	daomaster.expect("Your choice:")
	daomaster.sendline("6") # solve for 6 degrees of freedom
	daomaster.expect("Critical match-up radius:")
	daomaster.sendline("7") # play around with this

	if epoch_number == '01':

		for dither in range(start_dither+1,start_dither+5):
			#daomaster.expect(star_name + '_' + wavelength + '_e' + epoch + '_d' + str(dither) + '_cbcd_dn.alf')
			daomaster.sendline("")

	else:

		#daomaster.expect(star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.alf')
		#daomaster.sendline("")

		for dither in range(start_dither,start_dither+5):
			#daomaster.expect('.alf')
			#daomaster.expect(star_name + '_' + wavelength + '_e' + epoch + '_d' + str(dither) + '_cbcd_dn.alf')
			daomaster.sendline("")

	for dither in range(start_dither+1,start_dither+5):
		daomaster.expect(stem+'_d'+str(dither)+'_cbcd_dn.alf_all')
		daomaster.sendline("")

	# Repeat with decreasing match up size
	for match_up in range(7,-1,-1):
		daomaster.expect("New match-up radius")
		daomaster.sendline(str(match_up))

	# Options for different output files - only want transformations according to cookbook
	daomaster.expect("Assign new star IDs?")
	daomaster.sendline("y") # assign new ids so all frames have same ids - already done previously
	daomaster.expect("A file with mean magnitudes and scatter?")
	daomaster.sendline("n")
	daomaster.expect("A file with corrected magnitudes and errors?")
	daomaster.sendline("n")
	daomaster.expect("A file with raw magnitudes and errors?")
	daomaster.sendline("y")
	daomaster.expect("Output file name")
	daomaster.sendline(target_name + '_' + wavelength + '_f' + str(field) + '.cal')
	daomaster.expect("A file with the new transformations?")
	daomaster.sendline("e") # exits rest of options

	# If not epoch 1, then want to remove the 2 extra columns in the '.cal file'
	if epoch_number != '01':

		cal_data = pd.read_csv(target_name + '_' + wavelength + '_f' + str(field) + '.cal', skiprows=3, header=None, delim_whitespace=True, names=['ID', 'X', 'Y', 'E01_mag', 'E01_err', 'M1', 'E1', 'M2', 'E2', 'M3', 'E3', 'M4', 'E4', 'M5', 'E5'])
		cal_data.dropna(axis=0, subset=['M3'], inplace=True)
		cal_data.reset_index(drop=True, inplace=True)
		
		# Now drop 'EO1' columns
		cal_data.drop(['E01_mag', 'E01_err'], axis=1, inplace=True)

		# Write out to new file
		cal_data.to_csv(target_name + '_' + wavelength + '_f' + str(field) + '.cal', sep=' ', index=False)

	else:

		# Now just want to get rid of header so e01 and not e01 have same format of file
		cal_data = pd.read_csv(target_name + '_' + wavelength + '_f' + str(field) + '.cal', skiprows=3, header=None, delim_whitespace=True, names=['ID', 'X', 'Y', 'M1', 'E1', 'M2', 'E2', 'M3', 'E3', 'M4', 'E4', 'M5', 'E5'], usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12])
		cal_data.to_csv(target_name + '_' + wavelength + '_f' + str(field) + '.cal', sep=' ', index=False)

	return(0)


'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
				MAIN PART
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Read in list of stars that need converting
df = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star', 'Period', 'Channel'])

for i in range(0, len(df)):
	for j in [6]: # [1,6] when both fields work 

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

		if j == 1:
			field = 1
		elif j == 6:
			field = 2
		else: field = "Invalid field"

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

			stem = target_name + '_' + wavelength + '_e' + epoch_number

			# Change working directory to where the data is
			os.chdir(cwd)

			print "Working on: " + str(target_name) + "    Epoch: " + str(epoch_number) + "     Field: " + str(field)

			# Delete any previous runs
			if (os.path.isfile(stem+'_f'+str(field)+'_corrected.raw')):
				os.remove(stem+'_f'+str(field)+'_corrected.raw')		

			# Convert alf_apc file to standard format for DAOMATCH
			format_alf(i,j)

			# Do DAOMATCH and DAOMASTER
			daomaster(i,j)
