'''
Purpose: Create a medianed image from a set of BCDs. Will make a total of 5 medianed images.
1. Ch1 d1-5 
2. Ch2 d6-10 
3. Ch1 d6-10
4. Ch2 d1-5
5. Ch1 d6-10 and Ch2 d1-5 combined
Images 3-5 have target V* star present. Images 1 and 2 do not.
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
# images_to_make = [[1,1,5],[2,6,10],[1,6,10],[2,1,5],[[1,6,10],[2,1,5]]]


# Make medianed image function to call
def median(row_of_df):

	# MAKE HOME OF IMAGES CURRENT WORKING DIRECTORY

	# Grab galaxy, star name, channel and epoch

	galaxy = df['Galaxy'][i]
	target_name = df['Star'][i]

	if df['Channel'][i] == 1:
		wavelength = '3p6um'
	elif df['Channel'][i] == 2:
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

	# Change directory to where image is
	os.chdir(cwd)

	print "Working on star: " + target_name + "     ch: " + wavelength + "     epoch: " + epoch_number

	# Run DAOMATCH - put through the each BCD for the correct dithers
	# pexpect.spawn('')

	# Run DAOMASTER - refine the coordinate transformations from DAOMATCH

	# Run MONTAGE2 to actually make medianed image

	# Write down X and Y offsets

	# Add back in sky value 

	return(0)

# Set up data frame from txt file of stars (sys.argv[1]) to do it on
df = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star','Channel','Epoch'])

# # Loop over each row in txt file

for i in range(0, len(df)):
	median(i)

# for i in range(0, len(df)):

# 	# Make the 5 images described above
# 	for j in range(0, len(images_to_make)):

# 		if # element of images_to_make is only one so only one filter

# 		median(i, images_to_make[j][0], images_to_make[j][1], images_to_make[j][2])

# 		else # it is combining two filters




