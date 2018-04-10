'''
Purpose: Create a medianed image from a set of BCDs. Will make a total of 5 medianed images.
1. Ch1 d1-5 
2. Ch2 d6-10 
3. Ch1 d6-10
4. Ch2 d1-5
Images 4 & 5 have target V* star present. Images 1 and 2 do not.
Input is the same txt file 
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

# Import modules
import pandas as pd 
import sys
import pexpect
import os
import fnmatch

# Table of channels, start and end dithers for each image to be made
# images_to_make = [[1,1,5],[2,6,10],[1,6,10],[2,1,5]]

# Make medianed image function to call
def median(row_of_df, start_dither):

	# MAKE HOME OF IMAGES CURRENT WORKING DIRECTORY

	# Grab galaxy, star name, channel and epoch

	galaxy = df['Galaxy'][i]
	target_name = df['Star'][i]

	if df['Channel'][i] == 1:
		channel = 1
		wavelength = '3p6um'
	elif df['Channel'][i] == 2:
		channel = 2
		wavelength = '4p5um'
	else: wavelength = 'channel not defined'

	if df['Epoch'][i] < 10:
		epoch_number = '0' + str(df['Epoch'][i])
	else: epoch_number = str(df['Epoch'][i])

	# File to work on - this is the dither 1 image for this epoch
	# image = target_name + '_' + wavelength + '_e' + epoch_number + '_d1_cbcd_dn.fits'
	# image_nf =  target_name + '_' + wavelength + '_e' + epoch_number + '_d1_cbcd_dn'

    # Find absolute path of where images are
	home = '/home/ac833/Data/'

	cwd = home + str(galaxy) + '/BCD/' + target_name + '/ch' + str(df['Channel'][i]) + '/e' + str(epoch_number) + '/'
	stem = target_name + '_' + wavelength + '_e' + epoch_number

	# Change directory to where image is
	os.chdir(cwd)

	print "Working on star: " + target_name + "     ch: " + wavelength + "     epoch: " + epoch_number

	# Remove any previous runs of this particular script - needs completing
	#extensions = ['']
	#for ext in extensions:
	 #	for i in range(0,11):
	# 		if (os.path.isfile(target_name + '_' + wavelength + '_e' + epoch_number + '_d' + ))
	# 	if (os.path.isfile(image_nf+ext)):
	# 		os.remove(image_nf+ext)

	# Spawn DAOMATCH - put through the each BCD for the correct dithers
	daomatch = pexpect.spawn('daomatch')

	# Set up log file
	fout = file(stem+'_daomatch_log.txt','w')
	daomatch.logfile = fout

	# Give DAOMATCH all the 5 dithers to be used in making the medianed image
	daomatch.expect("Master input file:")
	daomatch.sendline(stem+'_d'+str(start_dither)+'_cbcd_dn.ap') # Give it the first BCD phot file
	daomatch.expect("Output file name")
	daomatch.sendline(stem+'_f1.mch')
	daomatch.expect("Next input file:")
	daomatch.sendline(stem+'_d'+str(start_dither + 1)+'_cbcd_dn.ap') # Give it the second BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline(stem+'_d'+str(start_dither + 2)+'_cbcd_dn.ap') # Give it the third BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline(stem+'_d'+str(start_dither + 3)+'_cbcd_dn.ap') # Give it the fourth BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline(stem+'_d'+str(start_dither + 4)+'_cbcd_dn.ap') # Give it the fifth BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline("") # exit

	print "DAOMATCH has made preliminary coordinate transformations"
	print "Checking how good they are..."

	# Open .mch file to check coefficients
	coeffs = pd.read_csv(stem+'_f1.mch', header=None, delim_whitespace=True, usecols=[2,3,4,5,6,7], names=['A', 'B', 'C', 'D', 'E', 'F'])

	if len(coeffs[(coeffs['C'] < 1.01) & (coeffs['C'] > 0.99)]) == 5:
		if len(coeffs[(coeffs['F'] < 1.01) & (coeffs['F'] > 0.99)]) == 5:
			if len(coeffs[(coeffs['D'] < 0.01) & (coeffs['D'] > -0.01)]) == 5:
				if len(coeffs[(coeffs['E'] < 0.01) & (coeffs['E'] > -0.01)]) == 5:
					print "All coefficients are good"
				else: 
					print "Coeff E is bad"
			else: 
				print "Coeff D is bad"
		else: 
			print "Coeff F is bad"
	else: 
		print "Coeff C is bad"

	# Run DAOMASTER - refine the coordinate transformations from DAOMATCH
	daomaster = pexpect.spawn('daomaster')

	# Set up log file
	fout = file(stem+'_daomaster_log.txt','w')
	daomaster.logfile = fout

	daomaster.expect("File with list of input files:")
	daomaster.sendline(stem+'_f1.mch')
	daomaster.expect("Minimum number, minimum fraction, enough frames:")
	daomaster.sendline("2, 0.5, 5") # play around with these values
	daomaster.expect("Maximum sigma:")
	daomaster.sendline("0.5") # play around with this value
	daomaster.expect("Your choice:")
	daomaster.sendline("6") # solve for 6 degrees of freedom
	daomaster.expect("Critical match-up radius:")
	daomaster.sendline("7") # play around with this

	for dither in range(start_dither+1,start_dither+5):
		daomaster.expect(stem+'_d'+str(dither)+'_cbcd_dn.ap')
		daomaster.sendline("")

	# Repeat with decreasing match up size
	for match_up in range(7,-1,-1):
		daomaster.expect("New match-up radius")
		daomaster.sendline(str(match_up))

	# Options for different output files - only want transformations according to cookbook
	daomaster.expect("Assign new star IDs?")
	daomaster.sendline("y") # assign new ids so all frames have same ids
	daomaster.expect("A file with mean magnitudes and scatter?")
	daomaster.sendline("n")
	daomaster.expect("A file with corrected magnitudes and errors?")
	daomaster.sendline("n")
	daomaster.expect("A file with raw magnitudes and errors?")
	daomaster.sendline("n")
	daomaster.expect("A file with the new transformations?")
	daomaster.sendline("y")
	daomaster.expect("Output file name")
	daomaster.sendline(stem+'_f1.mch_mast')
	daomaster.expect("A file with the transfer table?")
	daomaster.sendline("e") # exits rest of options

	# Run MONTAGE2 to actually make medianed image

	montage2 = pexpect.spawn('montage2')

	# Set up log file
	fout = file(stem+'_montage_log.txt','w')
	montage2.logfile = fout

	montage2.expect("File with transformations:")
	montage2.sendline(stem+'_f1.mch_mast')
	montage2.expect("Image-name suffix:")
	montage2.sendline("")
	montage2.expect("Minimum number of frames, percentile:")
	montage2.sendline("2,0.5") # play around with minimum number of frames
	montage2.expect("X limits of output image:")
	montage2.sendline("e")
	montage2.expect("Y limits of output image:")
	montage2.sendline("e")
	montage2.expect("Expansion factor:")
	montage2.sendline("1") # creates image with same scale as bcd images
	montage2.expect("Determine sky from overlap region?")
	montage2.sendline("n")
	montage2.expect("Name for output image")
	montage2.sendline(stem+'_f1')

	# Write down X and Y offsets
	log = open(stem+'_montage_log.txt', 'r')
	split = log.read().split()
	offsets = [split[-31], split[-30]]

	# Add back in sky value 
	# THIS IS GETTING THE ERROR OF FLOATING POINT INVALID OPERATION WHEN TRYING TO DO THIS THROUGH IRAF

	return(0)

# Set up data frame from txt file of stars (sys.argv[1]) to do it on
df = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star','Channel','Epoch'])

# # Loop over each row in txt file

for i in range(0, 2): # len(df)
	median(i, 1)

# for i in range(0, len(df)):

# 	# Make the 5 images described above
# 	for j in range(0, len(images_to_make)):

# 		if # element of images_to_make is only one so only one filter

# 		median(i, images_to_make[j][0], images_to_make[j][1], images_to_make[j][2])

# 		else # it is combining two filters




