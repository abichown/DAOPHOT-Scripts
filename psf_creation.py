'''
Purpose: Select candidate PSF stars from star list. Then let user manually create PSF model for 
one of the dithers. Then copy this PSF model file to the remaining 9 dithers.
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

# Import modules
import pandas as pd
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

	# File to work on - this is the dither 1 image for this epoch
	image = target_name + '_' + wavelength + '_e' + epoch_number + '_d1_cbcd_dn.fits'
	image_nf =  target_name + '_' + wavelength + '_e' + epoch_number + '_d1_cbcd_dn'

	print "Working on " + image_nf

    # Find absolute path of image
	home = '/home/ac833/Data/'
    
	for root,dirs,files in os.walk(home+galaxy):
		for filename in files:
			if fnmatch.fnmatch(filename,image):
				path_to_image = os.path.join(home,root) + '/'


	# Change directory to where image is
	os.chdir(path_to_image)

	# Remove any previous runs of this particular script
	extensions = ['.lst', '.nei', '.psf']
	for ext in extensions:
		if (os.path.isfile(image_nf+ext)):
			os.remove(image_nf+ext)

	# Open DAOPHOT
	daophot = pexpect.spawn('daophot')

	# Set up logfile
	fout = file(image_nf+'_daophot_log.txt','w')
	daophot.logfile = fout

	# Attach image
	daophot.expect("Command:")
	daophot.sendline("at " + image_nf)

	# Select candidate stars
	daophot.expect("Command:")
	daophot.sendline("pi")
	daophot.expect("Input file name")
	daophot.sendline("")
	daophot.expect("Desired number of stars, faintest magnitude:")
	daophot.sendline("25,19")
	daophot.expect("Output file name")
	daophot.sendline("")

	print "Candidate PSF stars chosen"

	# Close DAOPHOT
	daophot.expect("Command:")
	daophot.sendline("exit")
	daophot.close(force=True)

	# Create PSF model for Dither 1

	print "Now go and make PSF model for dither 1..."

	psf_done = raw_input("Type 'done' when PSF model is made: ")

	while psf_done != 'done':
		print "PSF not created.. go make it!"
		psf_done = raw_input("Type 'done' when PSF model is made: ")

	# Copy PSF model to other 9 dithers
	for j in range(2,11):
		new_name = image_nf.replace('d1','d'+str(j)) + '.psf'
		if (os.path.isfile(new_name)):
			os.remove(new_name)
		shutil.copy(image_nf+'.psf', new_name)

	print "PSF model made for " + image_nf + " and copied to all other dithers"

print "Complete"