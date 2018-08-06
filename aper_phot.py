'''
Purpose: Reads in list of star names, channel and epoch. For a given star and epoch, FIND stars 
and does aperture PHOT on them. Then choose the candidate PSF stars
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

import sys
import pexpect
import shutil
import os
import fnmatch
import re
import pandas as pd

import time
start = time.time()


# Open list of star names
# Format of the file is: GALAXY STAR_NAME CHANNEL EPOCH
stars = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star','Channel','Epoch'])

# Iterate over every line in text file
for i in range(0, len(stars)):


    '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		    INITIAL SETUP
    '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

    # Get target star and galaxy
    galaxy = stars['Galaxy'][i]
    target_name = stars['Star'][i]
    
    # Get channel and convert to wavelength
    if stars['Channel'][i] == 1:
    	channel = '1'
        wavelength = '3p6um'
    elif stars['Channel'][i] == 2:
    	channel = '2'
        wavelength = '4p5um'
    else: wavelength = 'channel not defined'
    
    # Get current epoch to be working with 
    if stars['Epoch'][i] < 10:
        epoch_number = '0' + str(stars['Epoch'][i])
    else: epoch_number = str(stars['Epoch'][i])
        
    # Iterate over each of the 10 dithers that make up one epoch
    for j in range(1,11):
        
	    # Dither number
	    dither_number = str(j)
	       
	    # File to work on
	    image = target_name + '_' + wavelength + '_e' + epoch_number + '_d' + dither_number + '_cbcd_dn.fits'
	    image_nf =  image.replace('.fits','') # name without fits ext

	    #print 'Currently working on file:' + image

	    # Where all my data is kept
	    home = '/home/ac833/Data/'

	    # Current working directory is the home of this image 
	    cwd = home + galaxy + '/BCD/' + target_name + '/ch' + channel + '/e' + epoch_number + '/'
	    os.chdir(cwd)

		# Remove any previous runs of this script
	    extensions = ['.coo', '.ap']
	    for ext in extensions:
		    if (os.path.isfile(image_nf+ext)):
			    os.remove(image_nf+ext)

		# Copy daophot options file to current working directory
	    if wavelength == '3p6um':
		    shutil.copy('/home/ac833/daophot-options-files/daophot.opt', 'daophot.opt')
	    if wavelength == '4p5um':
		    shutil.copy('/home/ac833/daophot-options-files/daophot.opt', 'daophot.opt')

		# Copy aperture photometry options files to current working directory
	    shutil.copy('/home/ac833/daophot-options-files/photo.opt', 'photo.opt')

	    '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
	     		RUN DAOPHOT
	    '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

	    daophot = pexpect.spawn("daophot")

		# Attach the image
	    daophot.expect("Command:")
	    daophot.sendline("at " + image_nf)

		# Find the stars
	    daophot.expect("Command:")
	    daophot.sendline("fi")
	    daophot.expect("Number of frames averaged, summed:")
	    daophot.sendline("1,1")
	    daophot.expect("File for positions")
	    daophot.sendline("")
	    daophot.expect("Are you happy with this?")
	    daophot.sendline("y")

		# Perform aperture photometry
	    daophot.expect("Command:")
	    daophot.sendline("ph")
	    daophot.expect("File with aperture radii")
	    daophot.sendline("")
	    daophot.expect("PHO>")
	    daophot.sendline("")
	    daophot.expect("Input position file")
	    daophot.sendline(image_nf + ".coo")
	    daophot.expect("Output file")
	    daophot.sendline(image_nf + ".ap")

	    # Choose candidate PSF stars
	    daophot.expect("Command:")
	    daophot.sendline("pi")
	    daophot.expect("Input file name")
	    daophot.sendline(image_nf + '.ap')
	    daophot.expect("Desired number of stars, faintest magnitude:")
	    daophot.sendline("40,19")
	    daophot.expect("Output file name")
	    daophot.sendline(image_nf + '.lst')   

		# Exit daophot
	    daophot.expect("Command:")
	    daophot.sendline("exit")
	    daophot.close(force=True)

end = time.time()
print(end - start)