'''
Purpose: Reads in list of star names, channel and epoch. For a given star and epoch, FIND stars 
and does aperture PHOT on them
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
stars = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star','Channel','Epoch'])

# Iterate over every line in text file
for i in range(0, len(stars)):


    '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		    INITIAL SETUP
    '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

    # Grab name of target star and galaxy
    galaxy = stars['Galaxy'][i]
    target_name = stars['Star'][i]
    
    # Grab channel
    if stars['Channel'][i] == 1:
        wavelength = '3p6um'
    elif stars['Channel'][i] == 2:
        wavelength = '4p5um'
    else: wavelength = 'channel not defined'
    
    # Grab epoch 
    if stars['Epoch'][i] < 10:
        epoch_number = '0' + str(stars['Epoch'][i])
    else: epoch_number = str(stars['Epoch'][i])
        
    for j in range(1,11):
        
	    # Dither number
	    dither_number = str(j)
	       
	    # File to work on
	    image = target_name + '_' + wavelength + '_e' + epoch_number + '_d' + dither_number + '_cbcd_dn.fits'
	    image_nf =  image.replace('.fits','')

	    #print 'Currently working on file:' + image

	    home = '/home/ac833/Data/'

	    for root,dirs,files in os.walk(home+galaxy):
		    for filename in files:
			    if fnmatch.fnmatch(filename,image):
				    path_to_image = os.path.join(home,root) + '/'

	    # print 'Changed directory to where the image is: ' + path_to_image
	    os.chdir(path_to_image)

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

	    # print 'Opening DAOPHOT...'

	    daophot = pexpect.spawn("daophot")

		# Attach the image
	    daophot.expect("Command:")
	    daophot.sendline("at " + image_nf)

	    # print 'Attached file: ' + image_nf

		# Find the stars
	    daophot.expect("Command:")
	    daophot.sendline("fi")
	    daophot.expect("Number of frames averaged, summed:")
	    daophot.sendline("1,1")
	    daophot.expect("File for positions")
	    daophot.sendline("")
	    daophot.expect("Are you happy with this?")
	    daophot.sendline("y")

	    # print "FIND complete"

		# Aperture photometry

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

	    # print "PHOT complete"

		# Exit daophot

	    daophot.expect("Command:")
	    daophot.sendline("exit")
	    daophot.close(force=True)


print 'FIND and PHOT done for all dithers at these epochs for these stars'

end = time.time()
print(end - start)