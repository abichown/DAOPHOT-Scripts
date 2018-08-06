
'''
Purpose: Make the PSF model from the stars selected in the previous script (star_selection.py)
Then run ALLSTAR to get PSF magnitudes for all the stars
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

		# Open DAOPHOT
		daophot = pexpect.spawn('daophot')

		# Attach image
		daophot.expect("Command:")
		daophot.sendline("at " + filename)

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

		# Open ALLSTAR
		allstar = pexpect.spawn('allstar')

		allstar.expect("OPT>")
		allstar.sendline("")

		allstar.expect("Input image name:")
		allstar.sendline(filename)

		allstar.expect("File with the PSF")
		allstar.sendline("")
		allstar.expect("Input file")
		allstar.sendline("")
		allstar.expect("File for results")
		allstar.sendline("")
		allstar.expect("Name for subtracted image")
		allstar.sendline("")

		# Close ALLSTAR
		allstar.expect("Good bye")
		allstar.close(force=True)




