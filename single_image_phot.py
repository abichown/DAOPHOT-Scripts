'''
Purpose: Carry out photometry on a single image
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
if channel == '1' or channel == '3p6' or channel == '3p6um' or channel == 'ir1':
	wavelength = '3p6um'
elif channel == '2' or channel == '4p5' or channel == '4p5um' or channel == 'ir2':
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

home = '/home/ac833/HV00872_comparison/' # for current project
#home = '/home/ac833/Data/' # for project in general
    
directories = os.listdir(home)

for directory in directories:
	for root,dirs,files in os.walk(home+directory):
		for filename in files:
			if fnmatch.fnmatch(filename,image): 
				path_to_image = os.path.join(home,root) +'/'

print 'Changed directory to where the image is: ' + path_to_image
os.chdir(path_to_image)

# Remove any previous runs
extensions = ['.coo', '.ap', '.lst', '.nei', '.psf', '.als', 's.fits', '_allstar_log.txt', '_daophot_log.txt']
for ext in extensions:
	if (os.path.isfile(image_nf+ext)):
		os.remove(image_nf+ext)

# Copy daophot options file to current working directory
if channel == '1':
	shutil.copy('/home/ac833/daophot-options-files/daophot.opt', 'daophot.opt')
if channel == '2':
	shutil.copy('/home/ac833/daophot-options-files/daophot.opt', 'daophot.opt')

# Copy aperture photometry and allstar options files to current working directory
shutil.copy('/home/ac833/daophot-options-files/photo.opt', 'photo.opt')
shutil.copy('/home/ac833/daophot-options-files/allstar.opt', 'allstar.opt')

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		RUN DAOPHOT TO CREATE PSF MODEL
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

print 'Opening DAOPHOT...'

daophot = pexpect.spawn("daophot")

# Set up log file
fout = file(image_nf+'_daophot_log.txt','w')
daophot.logfile = fout

# Attach the image
daophot.expect("Command:")
daophot.sendline("at " + image_nf)

print 'Attached file: ' + image_nf

# Find the stars

daophot.expect("Command:")
daophot.sendline("fi")
daophot.expect("Number of frames averaged, summed:")
daophot.sendline("1,1")
daophot.expect("File for positions")
daophot.sendline("")
daophot.expect("Are you happy with this?")
daophot.sendline("y")

print "FIND complete"

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

print "PHOT complete"

# Pick candidate PSF stars

daophot.expect("Command:")
daophot.sendline("pi")
daophot.expect("Input file name")
daophot.sendline("")
daophot.expect("Desired number of stars, faintest magnitude:")
daophot.sendline("25,19")
daophot.expect("Output file name")
daophot.sendline("")

print "Candidate PSF stars chosen"

# Exit daophot

daophot.expect("Command:")
daophot.sendline("exit")
daophot.close(force=True)

# PSF MODEL CREATION
# This part is still to be done by hand for now 

print "Go to DAOPHOT and make PSF model..."

psf_done = raw_input("Type 'done' when PSF model created: ")

while psf_done != 'done':
	print "PSF not created... go make the PSF!"
	psf_done = raw_input("Type 'done' when PSF model created: ")

print "PSF model created"


'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		RUN ALLSTAR ON ALL STARS IN IMAGE
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

allstar = pexpect.spawn("allstar")

# Set up log file
fout = file(image_nf+'_allstar_log.txt','w')
allstar.logfile = fout

allstar.expect('OPT>')
allstar.sendline("")
allstar.expect('Input image name:')
allstar.sendline(image_nf)
allstar.expect('File with the PSF')
allstar.sendline("")

# Fit PSF model to all stars
allstar.expect('Input file')
allstar.sendline("")

# Creates .als file
allstar.expect('File for results')
allstar.sendline("")
# allstar.expect('New output file name')
# allstar.sendline("")

# Creates _dns file
allstar.expect('Name for subtracted image')
allstar.sendline("")

# Exit allstar
allstar.expect("Good bye.")
allstar.close(force=True)

print "All stars now have psf magnitudes and have been subtracted"


'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		DELETE FILES NO LONGER NEEDED
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# # Remove any files not needed
# extensions = ['.coo', '.nei', '.psf', '_allstar_log.txt', '_daophot_log.txt']
# for ext in extensions:
# 	if (os.path.isfile(image_nf+ext)):
# 		os.remove(image_nf+ext)	

print "Photometry on " + image_nf + " complete"