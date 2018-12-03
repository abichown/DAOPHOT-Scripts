'''
Purpose: Perform PSF photometry using ALLSTAR on the individual BCD images.
The PSF used is created in master_image.py using PSF stars chosen from the master image.
Outputs are individual BCD .als files
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

# Import modules
import pandas as pd
import os
import pexpect
import shutil
import sys

# Load stars and epochs to perform photometry on
stars = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star','Period', 'RA', 'Dec', 'Channel'])


# Iterate over every line in text file
for i in range(0, len(stars)):

    # Get target star and galaxy
    galaxy = str(stars['Galaxy'][i])
    star_name = str(stars['Star'][i])
    
    # Get channel and convert to wavelength
    if stars['Channel'][i] == 1:
    	channel = '1'
        wavelength = '3p6um'
    elif stars['Channel'][i] == 2:
    	channel = '2'
        wavelength = '4p5um'
    else: wavelength = 'channel not defined'

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

		# Where all my data is kept
		home = '/home/ac833/Data/'

		# Current working directory is the home of this image 
		cwd = home + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch_number + '/'
		os.chdir(cwd)

		# Copy allstar.opt file to cwd 
		shutil.copy('/home/ac833/daophot-options-files/allstar.opt', 'allstar.opt')

		# For each field (currently only doing for field 2)
		for field in [2]: #[1,2]

			if field == 1:
				start_dither = 1
			else: start_dither = 6

			for dither in range(start_dither, start_dither+5):

				# Spawn ALLSTAR and perform PSF photometry using PSF model made from master image
				allstar = pexpect.spawn('allstar')

				fout = file('allstar_log.txt', 'w')
				allstar.logfile = fout

				allstar.expect("OPT>")
				allstar.sendline("")
				allstar.expect("Input image name:")
				allstar.sendline(star_name + '_' + wavelength + '_e'+ epoch_number + '_d' + str(dither) + '_cbcd_dn.fits')
				allstar.expect("File with the PSF")
				allstar.sendline(star_name + '_' + wavelength + '_f' + str(field) + '_master.psf')

				allstar.expect("Input file")
				allstar.sendline(star_name + '_' + wavelength + '_e'+ epoch_number + '_d' + str(dither) + '_cbcd_dn.ap')
				allstar.expect("File for results")
				allstar.sendline(star_name + '_' + wavelength + '_e'+ epoch_number + '_d' + str(dither) + '_cbcd_dn.als')
				allstar.expect("Name for subtracted image")
				allstar.sendline(star_name + '_' + wavelength + '_e'+ epoch_number + '_d' + str(dither) + '_cbcd_dns.fits')

				allstar.expect("Good bye")
				allstar.close(force=True)			

