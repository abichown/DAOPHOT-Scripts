'''
Purpose: Take the .alf_apc magnitudes and put them into a file suitable for DAOMATCH/DAOMASTER.
Then run them through daomatch and daomaster to get file needed for averaging magnitudes
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

import pandas as pd 
import os
import sys
import pexpect

def format_alf(row_of_df, start_dither):

	# Format all 5 dithers for this star, epoch and field
	for i in range(start_dither, start_dither+5):

		# Get filenames
		alf = stem + '_d' + str(i) + '_cbcd_dn.alf' # has everything but the correct magnitudes
		alf_apc = stem + '_d' + str(i) + '_cbcd_dn.alf_apc' # has the magnitudes we want

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


def dao(row_of_df, start_dither):

	# Spawn DAOMATCH
	daomatch = pexpect.spawn('daomatch')

	# Give DAOMATCH all the 5 dithers to be used in making the medianed image
	daomatch.expect("Master input file:")
	daomatch.sendline(stem+'_d'+str(start_dither)+'_cbcd_dn.alf_all') # Give it the first BCD phot file
	daomatch.expect("Output file name")
	daomatch.sendline(stem+'_f'+str(field)+'_corrected.mch')
	daomatch.expect("Next input file")
	daomatch.sendline(stem+'_d'+str(start_dither + 1)+'_cbcd_dn.alf_all') # Give it the second BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline(stem+'_d'+str(start_dither + 2)+'_cbcd_dn.alf_all') # Give it the third BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline(stem+'_d'+str(start_dither + 3)+'_cbcd_dn.alf_all') # Give it the fourth BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline(stem+'_d'+str(start_dither + 4)+'_cbcd_dn.alf_all') # Give it the fifth BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline("") # exit

	# Run DAOMASTER - refine the coordinate transformations from DAOMATCH
	daomaster = pexpect.spawn('daomaster')

	daomaster.expect("File with list of input files:")
	daomaster.sendline(stem+'_f'+str(field)+'_corrected.mch')
	daomaster.expect("Minimum number, minimum fraction, enough frames:")
	daomaster.sendline("2, 0.5, 5") # play around with these values
	daomaster.expect("Maximum sigma:")
	daomaster.sendline("0.5") # play around with this value
	daomaster.expect("Your choice:")
	daomaster.sendline("6") # solve for 6 degrees of freedom
	daomaster.expect("Critical match-up radius:")
	daomaster.sendline("7") # play around with this

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
	daomaster.sendline(stem + '_f' + str(field) + '_corrected.raw')
	daomaster.expect("A file with the new transformations?")
	daomaster.sendline("e") # exits rest of options

	return(0)


'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
				MAIN PART
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Read in list of stars that need converting
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

		print "Working on: " + str(target_name) + "    Epoch: " + str(epoch_number) + "     Field: " + str(field)

		# Delete any previous runs
		if (os.path.isfile(stem+'_f'+str(field)+'_corrected.raw')):
			os.remove(stem+'_f'+str(field)+'_corrected.raw')		

		# Convert alf_apc file to standard format for DAOMATCH
		format_alf(i,j)

		# Do DAOMATCH and DAOMASTER
		dao(i,j)
