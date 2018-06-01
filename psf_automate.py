'''
Purpose: Select candidate PSF stars from star list. Then run PSF once to find the residuals.
Discard stars whose residual > 0.040. Update PSF star list. Then rerun PSF with new star list.
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

# Import modules
import pandas as pd
import numpy as np
import pexpect
import sys
import fnmatch
import os
import shutil

import time
start = time.time()

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

# find files 
stars = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star','Channel','Epoch'])

# Iterate over every line in text file
for i in range(0, len(stars)):

	# Grab galaxy, star name, channel and epoch
	galaxy = stars['Galaxy'][i]
	target_name = stars['Star'][i]

	if stars['Channel'][i] == 1:
		wavelength = '3p6um'
	elif stars['Channel'][i] == 2:
		wavelength = '4p5um'
	else: wavelength = 'channel not defined'

	if stars['Epoch'][i] < 10:
		epoch_number = '0' + str(stars['Epoch'][i])
	else: epoch_number = str(stars['Epoch'][i])

    # Find absolute path of image	
	home = '/home/ac833/Data/'
    
	cwd = home + str(galaxy) + '/BCD/' + target_name + '/ch' + str(stars['Channel'][i]) + '/e' + str(epoch_number) + '/'
	stem = target_name + '_' + wavelength + '_e' + epoch_number

	# Change directory to where image is
	os.chdir(cwd)

	print "Working on: " + stem

	# Loop over all 10 dithers
	for i in range(1,11):

		filename = stem + '_d' + str(i) + '_cbcd_dn' # get filename to use for this dither

		# Remove any previous runs of this particular script
		extensions = ['.lst', '.nei', '.psf', '_psf_log.txt']
		for ext in extensions:
			if (os.path.isfile(filename+ext)):
				os.remove(filename+ext)

		# Open DAOPHOT
		daophot = pexpect.spawn('daophot')

		# Set up PSF log file to be used to find residuals
		fout = file(filename+'_psf_log.txt','w')
		daophot.logfile = fout

		# Attach image
		daophot.expect("Command:")
		daophot.sendline("at " + filename)

		# Select candidate stars
		daophot.expect("Command:")
		daophot.sendline("pi")
		daophot.expect("Input file name")
		daophot.sendline("")
		daophot.expect("Desired number of stars, faintest magnitude:")
		daophot.sendline("40,20")
		daophot.expect("Output file name")
		daophot.sendline("")

		# Run PSF to get residuals outputted to logfile
		daophot.expect("Command:")
		daophot.sendline("psf")
		daophot.expect("File with aperture results")
		daophot.sendline("")
		daophot.expect("File with PSF stars")
		daophot.sendline("")
		daophot.expect("File for the PSF")
		daophot.sendline("")

		# Close DAOPHOT
		daophot.expect("Command:")
		daophot.sendline("exit")
		daophot.close(force=True)	

		# Delete stars whose residuals are higher than 0.040

		# Obtain list of good PSF stars by going through lst file and finding ones with residuals < 0.040
		f = open(filename+'_psf_log.txt', 'r')
		lines = f.read().splitlines()
		f.close()

		residuals = lines[-14:-6] # obtain the residuals part of the file

		new_residuals = residuals[0].split() + residuals[1].split() + residuals[2].split() + residuals[3].split() + residuals[4].split() + residuals[5].split() + residuals[6].split() + residuals[7].split()

		# Remove occurences of '?' and '*' symbols
		new_residuals = list(filter(lambda x: x != '?', new_residuals))
		new_residuals = list(filter(lambda x: x != '*', new_residuals))

		# NEED TO WORK ON FROM HERE...

		# convert list to np array and then reshape to basically get table of ID with residuals
		res_df = pd.DataFrame(np.array(new_residuals).reshape(len(new_residuals)/2,2), columns = ['ID', 'Residual'])

		# Drop columns whose residual is saturated or defective
		res_df.drop(res_df[res_df.Residual == 'saturated'].index, inplace=True)
		res_df.drop(res_df[res_df.Residual == 'defective'].index, inplace=True)

		res_df.reset_index(inplace=True, drop=True) # then reset index

		# List of star IDs to keep
		keep_stars = []

		for j in range(0, len(res_df)):
			if float(res_df['Residual'][j]) < 0.040:
				keep_stars.append(res_df['ID'][j])


		# KEEP STARS NOW HAS ALL THE IDS OF THE STARS TO BE USED IN THE NEW PSF MODEL

		# NOW MODIFY THE LST FILE TO ONLY KEEP THESE GOOD PSF STARS

		df = pd.read_csv(filename+".lst", skiprows=3, header=None, delim_whitespace=True, names=['ID', 'X', 'Y', 'Mag', 'Error'])

		# Get the index values of the stars to be used
		index_list = []

		for ID in keep_stars:
			for i in range(0, len(df)):
				if df.iloc[i]['ID'] == int(ID):
					index_list.append(i)

		# Now create new df with only these stars in
		new_df = df.filter(index_list, axis=0)

		# WRITE OUT TO NEW LST FILE

		# Get header of lst file
		f = open(filename+".lst", 'r')
		header = f.read().splitlines()[0:3]
		f.close()

		# Now overwrite lst file
		f = open(filename+".lst", 'w')
		f.writelines(header[0] + '\n' + header[1] + '\n' + header[2] + '\n')

		new_df.to_csv(f, sep=' ', mode='a', header=None, index=False)

		f.close()

		# Delete psf log - don't need anymore
		os.remove(filename+'_psf_log.txt')

	# THEN NEED TO REDO PSF MODEL WITH THESE STARS ONLY

	# Choose the one with the most PSF stars as the final psf model to be made

	curr_max = 0

	for i in range(1,11):

		filename = stem + '_d' + str(i) + '_cbcd_dn.lst'

		if file_len(filename) > curr_max:
			winning_dither = i

	winning_file = stem + '_d' + str(winning_dither) + '_cbcd_dn'

	# Open DAOPHOT
	daophot = pexpect.spawn('daophot')

	# Attach image
	daophot.expect("Command:")
	daophot.sendline("at " + winning_file+'.fits')

	# Get final psf mode
	daophot.expect("Command:")
	daophot.sendline("psf")
	daophot.expect("File with aperture results")
	daophot.sendline("")
	daophot.expect("File with PSF stars")
	daophot.sendline("")
	daophot.expect("File for the PSF")
	daophot.sendline("")
	daophot.expect("New output file name")
	daophot.sendline("")
	daophot.expect("New output file name")
	daophot.sendline("")

	# Close DAOPHOT
	daophot.expect("Command:")
	daophot.sendline("exit")
	daophot.close(force=True)	

	# Copy the winning dither psf model to the other 9 dithers
	for j in range(1,11):

		psf_model = winning_file+'.psf'
		filename = stem + '_d' + str(j) + '_cbcd_dn.psf'

		os.remove(stem + '_d' + str(j) + '_cbcd_dn.nei')

		if psf_model != filename:
			shutil.copy(psf_model, filename)

end = time.time()
print(end - start)