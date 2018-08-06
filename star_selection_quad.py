'''
Purpose: Select candidate PSF stars from star list. Split stars into quadrants and then choose the two 
brightest stars from each quadrant to be the PSF stars... so 8 stars in total per dither.
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
	for index in range(1,11):

		filename = stem + '_d' + str(index) + '_cbcd_dn' # get filename to use for this dither

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

		# Close DAOPHOT
		daophot.expect("Command:")
		daophot.sendline("exit")
		daophot.close(force=True)	

		print "Working on " + filename

		# Open LST file and extract stars
		df = pd.read_csv(filename+".lst", header=None, names=['ID', 'X', 'Y', 'Mag', 'Error'], skiprows=3, delim_whitespace=True)

		# Bin stars into quadrants 

		# 4 quadrants of image
		q1 = [] # upper-left i.e. x in [0,128] and y in [128,256]
		q2 = [] # upper-right i.e. x in [128,256] and y in [128,256]
		q3 = [] # lower-left i.e. x in [0,128] and y in [0, 128]
		q4 = [] # lower-right i.e. x in [128,256] and y in [0,128]

		for j in range(0, len(df)):

		 	if df['X'][j] < 128:
		 		if df['Y'][j] < 128:
		 			q3.append(df['ID'][j])
		 		else: 
		 			q1.append(df['ID'][j])
		 	else:
		 		if df['Y'][j] < 128:
		 			q4.append(df['ID'][j])
		 		else:
		 			q2.append(df['ID'][j])


		# Create dataframe for each of the quadrants
		q1_mag = []
		q1_x = []
		q1_y = []
		q1_err = []
		q2_mag = []
		q2_x = []
		q2_y = []
		q2_err = []
		q3_mag = []
		q3_x = []
		q3_y = []
		q3_err = []
		q4_mag = []
		q4_x = []
		q4_y = []
		q4_err = []

		for star_id in q1:
			for k in range(0, len(df)):
				if df.iloc[k]['ID'] == star_id:
					q1_mag.append(round(df['Mag'][k], 3))
					q1_x.append(round(df['X'][k], 2))
					q1_y.append(round(df['Y'][k], 2))
					q1_err.append(round(df['Error'][k], 3))

		for star_id in q2:
			for k in range(0, len(df)):
				if df.iloc[k]['ID'] == star_id:
					q2_mag.append(round(df['Mag'][k], 3))
					q2_x.append(round(df['X'][k], 2))
					q2_y.append(round(df['Y'][k], 2))
					q2_err.append(round(df['Error'][k], 3))

		for star_id in q3:
			for k in range(0, len(df)):
				if df.iloc[k]['ID'] == star_id:
					q3_mag.append(round(df['Mag'][k], 3))
					q3_x.append(round(df['X'][k], 2))
					q3_y.append(round(df['Y'][k], 2))
					q3_err.append(round(df['Error'][k], 3))


		for star_id in q4:
			for k in range(0, len(df)):
				if df.iloc[k]['ID'] == star_id:
					q4_mag.append(round(df['Mag'][k], 3))
					q4_x.append(round(df['X'][k], 2))
					q4_y.append(round(df['Y'][k], 2))
					q4_err.append(round(df['Error'][k], 3))


		q1_df = pd.DataFrame({'ID': q1, 'X': q1_x, 'Y': q1_y, 'Mag': q1_mag, 'Error': q1_err})
		q2_df = pd.DataFrame({'ID': q2, 'X': q2_x, 'Y': q2_y, 'Mag': q2_mag, 'Error': q2_err})
		q3_df = pd.DataFrame({'ID': q3, 'X': q3_x, 'Y': q3_y, 'Mag': q3_mag, 'Error': q3_err})
		q4_df = pd.DataFrame({'ID': q4, 'X': q4_x, 'Y': q4_y, 'Mag': q4_mag, 'Error': q4_err})

		# Choose the two brighest stars from each quadrant to be PSF stars

		psf_list = pd.concat([q1_df[0:4], q2_df[0:4], q3_df[0:4], q4_df[0:4]], axis=0)

		# Change order of columns
		cols = ['ID', 'X', 'Y', 'Mag', 'Error']
		psf_list = psf_list[cols]

		# Write out in LST file format

		# Get header of LST file
		f = open(filename+'.lst', 'r')
		header = f.read().splitlines()[0:3]
		f.close()

		# Now overwrite LST file
		f = open(filename+'.lst', 'w')
		f.writelines(header[0] + '\n' + header[1] + '\n' + header[2] + '\n')

		psf_list.to_csv(f, sep=' ', mode='a', header=None, index=False)

		f.close()

		# Delete log file
		os.remove(filename+'_psf_log.txt')


end = time.time()
print(end - start)