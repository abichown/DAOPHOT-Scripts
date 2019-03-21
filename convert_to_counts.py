'''
Purpose: Convert all BCD images (ending in _cbcd.fits) in a given directory to counts.
It uses the keywords in the FITS image to do the conversion.
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

import re
import shutil
import os
import fnmatch
from astropy.io import fits

# Directory to look in
home = '/home/ac833/Data/SMC/BCD/'
os.chdir(home)

# Find all directories in the home directory
directories = os.listdir(home)

# Walk over all directories to find files to convert 
for directory in directories:
	for root,dirs,files in os.walk(directory):
		for filename in files:

			# If filename matches with _cbcd.fits then convert
			if fnmatch.fnmatch(filename, '*_cbcd.fits'):

				# This is the full path to the file
				full_filename = os.path.join(home,root,filename)

				print full_filename

				# Proceed to convert image
				new_name = re.sub('.fits', '_dn.fits', full_filename)

				if os.path.isfile(new_name) == False:
					shutil.copy(full_filename, new_name)
					hdulist = fits.open(new_name, mode='update')
					prihdr = hdulist[0].header
					scidata = hdulist[0].data 
					exptime = prihdr['exptime']
					fluxconv = prihdr['fluxconv']
					scidata *= exptime/fluxconv # this part converts to counts