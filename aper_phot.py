'''
Purpose: Really a practice of how to control DAOPHOT with Python but also script to
do initial aperture photometry
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

import sys
import pexpect
import shutil
import os
import fnmatch
import re

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		INITIAL SETUP
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Info entered on command line to identify the image to work on
target_name = sys.argv[1]
channel = sys.argv[2]
epoch_number = sys.argv[3]
dither_number = sys.argv[4]

# Convert channel number to wavelength
if channel == '1':
	wavelength = '3p6um'
elif channel == '2':
	wavelength = '4p5um'
else: wavelength = 'channel not defined'

# Get epoch number into right format
if len(epoch_number) == 1:
	epoch_number = '0' + epoch_number

# File to work on
image = target_name + '_' + wavelength + '_e' + epoch_number + '_d' + dither_number + '_cbcd_dn.fits'
image_nf =  target_name + '_' + wavelength + '_e' + epoch_number + '_d' + dither_number + '_cbcd_dn'

print 'Currently working on file:' + image

## Get absolute path of image as need to move to this directory

home = '/home/ac833/Data/'
    
directories = os.listdir(home)

for directory in directories:
	for root,dirs,files in os.walk(home+directory):
		for filename in files:
			if fnmatch.fnmatch(filename,image): 
				path_to_image = os.path.join(home,root) +'/'

print 'Changed directory to where the image is: ' + path_to_image
os.chdir(path_to_image)

# Remove any previous runs
extensions = ['.coo', '.ap', '.lst', '.nei', '.psf', '.als', 's.fits']
for ext in extensions:
	if (os.path.isfile(image_nf+ext)):
		os.remove(image_nf+ext)

# Copy daophot options file to current working directory
if channel == '1':
	shutil.copy('/home/ac833/daophot-options-files/daophot.opt', 'daophot.opt')
if channel == '2':
	shutil.copy('/home/ac833/daophot-options-files/daophot.opt', 'daophot.opt')


'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		RUN DAOPHOT
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

print 'Opening DAOPHOT...'

daophot = pexpect.spawn("daophot")

# Attach the image
daophot.expect("Command:")
daophot.sendline("AT " + image_nf)

print 'Attached file: ' + image_nf

# Find the stars
daophot.expect("Command:")
daophot.sendline("FI")
daophot.expect("Number of frames averaged, summed:")
daophot.sendline("1,1")
daophot.expect("File for positions")
daophot.sendline("")
daophot.expect("Are you happy with this?")
daophot.sendline("y")

print "FIND complete"

# Aperture photometry

