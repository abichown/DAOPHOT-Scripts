'''
Purpose: Create a medianed image from a set of BCDs. Will make 2 medianed images for each row.
of the star list: one for dithers 1-5 and one for dithers 6-10. Then performs allframe 
photometry. Final file that you want is the .alf file for the field which has the psf mags
for each star in each frame.
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

# Import modules
import pandas as pd 
import sys
import pexpect
import os
import fnmatch
import shutil

import matplotlib.pyplot as plt
from kneed import DataGenerator, KneeLocator

# Make medianed image function to call
def median(row_of_df, start_dither):

	print "Making median image of star: " + target_name + "     ch: " + wavelength + "     epoch: " + epoch_number

	# Spawn DAOMATCH - put through the each BCD for the correct dithers
	daomatch = pexpect.spawn('daomatch')

	# Set up log file
	fout = file(stem+'_daomatch_log.txt','w')
	daomatch.logfile = fout

	# Give DAOMATCH all the 5 dithers to be used in making the medianed image
	# The '/' after the file names lets DAOMATCH know that the scale is the same and it is only the rotation and shifts that might have changed
	daomatch.expect("Master input file:")
	daomatch.sendline(stem+'_d'+str(start_dither)+'_cbcd_dn.als') # Give it the first BCD phot file
	daomatch.expect("Output file name")
	daomatch.sendline(stem+'_f'+field+'.mch')
	daomatch.expect("Next input file:")
	daomatch.sendline(stem+'_d'+str(start_dither + 1)+'_cbcd_dn.als/') # Give it the second BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline(stem+'_d'+str(start_dither + 2)+'_cbcd_dn.als/') # Give it the third BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline(stem+'_d'+str(start_dither + 3)+'_cbcd_dn.als/') # Give it the fourth BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline(stem+'_d'+str(start_dither + 4)+'_cbcd_dn.als/') # Give it the fifth BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline("") # exit

	# Run DAOMASTER - refine the coordinate transformations from DAOMATCH
	daomaster = pexpect.spawn('daomaster')

	# Set up log file
	#fout = file(stem+'_daomaster_log.txt','w')
	#daomaster.logfile = fout

	daomaster.expect("File with list of input files:")
	daomaster.sendline(stem+'_f'+field+'.mch')
	daomaster.expect("Minimum number, minimum fraction, enough frames:")
	daomaster.sendline("2, 0.5, 5") # play around with these values
	daomaster.expect("Maximum sigma:")
	daomaster.sendline("0.5") # play around with this value
	daomaster.expect("Your choice:")
	daomaster.sendline("6") # solve for 6 degrees of freedom
	daomaster.expect("Critical match-up radius:")
	daomaster.sendline("7") 
	
	for dither in range(start_dither+1,start_dither+5):
		daomaster.expect(stem+'_d'+str(dither)+'_cbcd_dn.als')
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
	daomaster.sendline("y")
	daomaster.expect("Output file name")
	daomaster.sendline(stem + '_f' + field + '.raw')
	daomaster.expect("A file with the new transformations?")
	daomaster.sendline("y")
	daomaster.expect("Output file name")
	daomaster.sendline(stem+'_f'+field+'.mch')
	daomaster.expect("New output file name")
	daomaster.sendline("")
	daomaster.expect("A file with the transfer table?")
	daomaster.sendline("e") # exits rest of options

	# Run MONTAGE2 to actually make medianed image

	montage2 = pexpect.spawn('montage2')

	# Set up log file - need this for the offsets
	fout = file(stem+'_montage_log.txt','w')
	montage2.logfile = fout

	montage2.expect("File with transformations:")
	montage2.sendline(stem+'_f'+field+'.mch')
	montage2.expect("Image-name suffix:")
	montage2.sendline("")
	montage2.expect("Minimum number of frames, percentile:")
	montage2.sendline("1,0.5") # play around with minimum number of frames
	montage2.expect("X limits of output image:")
	montage2.sendline("e")
	montage2.expect("Y limits of output image:")
	montage2.sendline("e")
	montage2.expect("Expansion factor:")
	montage2.sendline("1") # creates image with same scale as bcd images
	montage2.expect("Determine sky from overlap region?")
	montage2.sendline("y")
	montage2.expect("Name for output image")
	montage2.sendline(stem+'_f'+field)

	# Write down X and Y offsets
	log = open(stem+'_montage_log.txt', 'r')
	split = log.read().split()
	offsets = [split[-31], split[-30]]

	# Add back in sky value 
	# THIS IS GETTING THE ERROR OF FLOATING POINT INVALID OPERATION WHEN TRYING TO DO THIS THROUGH IRAF

	return(offsets)

# Find stars of medianed image
def find_stars(row_of_df, start_dither):

	print "Finding stars of " + target_name + "     ch: " + wavelength + "     epoch: " + epoch_number

	# Run DAOPHOT
	daophot = pexpect.spawn('daophot') 

	# Set up logfile
	fout = file(stem+'_daophot_log.txt','w')
	daophot.logfile = fout

	# Attach medianed image
	daophot.expect("Command:")
	daophot.sendline("at " + stem + '_f' + field + '.fits')

	# Experiment with threshold to find the value near the 'elbow' of the curve
	daophot.expect("Command:")
	daophot.sendline("opt")
	daophot.expect("File with parameters")
	daophot.sendline("")
	daophot.expect("OPT>")
	daophot.sendline("th=2") # set initial threshold to 2
	daophot.expect("OPT>")
	daophot.sendline("")

	daophot.expect("Command:")
	daophot.sendline("fi")
	daophot.expect("Number of frames averaged, summed:")
	daophot.sendline("5,1") # for this work with 5 dithers making the median image
	daophot.expect("File for positions")
	daophot.sendline("")

	# Needs to store the number of detections for the threshold
	for threshold in range(3,20):

		daophot.expect("Are you happy with this?")
		daophot.sendline("n")
		daophot.expect("New threshold")
		daophot.sendline(str(threshold))
		daophot.expect("Output file name")
		daophot.sendline("")
		daophot.expect("New output file name")
		daophot.sendline("")

	# Last threshold = 20
	daophot.expect("Are you happy with this?")
	daophot.sendline("y")

	# Open log file, extract number of detections at each threshold
	# Then work out the ratio between two consecutive detections
	# Once ratio exceeds 0.6, this is the desired threshold
	log = open(stem+'_daophot_log.txt', 'r')
	split = log.read().split()

	num_det = []

	for i in range(0,18):
		index = -8 -i*40
		num_det.append(split[index])

	# Reverse order so threshold increases and detections decrease
	num_det = num_det[::-1]
	thresholds = range(2,20)

	# plt.clf()

	# plt.scatter(thresholds, num_det)
	# plt.show()

	# Package that finds the elbow/knee of a curve - this will be the optimal detection threshold
	kn = KneeLocator(thresholds, num_det, direction='decreasing')
	#print "Knee = " + str(kn.knee)


	# # Now find optimal threshold near elbow of curve
	# ratio = 0
	# i = 0

	# while ratio < 0.6:
	# 	i += 1
	# 	ratio = float(num_det[i])/float(num_det[i-1])
	# 	threshold = i+2 # this is the desired threshold

	# print "threshold = " + str(threshold)

	# Then run FIND with this desired threshold
	daophot.expect("Command:")
	daophot.sendline("opt")
	daophot.expect("File with parameters")
	daophot.sendline("")
	daophot.expect("OPT>")
	daophot.sendline("th="+str(kn.knee)) # now set to best threshold
	daophot.expect("OPT>")
	daophot.sendline("")

	daophot.expect("Command:")
	daophot.sendline("fi")
	daophot.expect("Number of frames averaged, summed:")
	daophot.sendline("5,1") # for this work with 5 dithers making the median image
	daophot.expect("File for positions")
	daophot.sendline("")
	daophot.expect("New output file name")
	daophot.sendline("")
	daophot.expect("Are you happy with this?")
	daophot.sendline("y")

	# Run PHOT
	daophot.expect("Command:")
	daophot.sendline("ph")
	daophot.expect("File with aperture radii")
	daophot.sendline("")
	daophot.expect("PHO>")
	daophot.sendline("")
	daophot.expect("Input position file")
	daophot.sendline("")
	daophot.expect("Output file")
	daophot.sendline("")

	# Close DAOPHOT
	daophot.expect("Command:")
	daophot.sendline("ex")
	daophot.close(force=True)

	# Copy allstar option file if it doesn't already exist in cwd
	shutil.copy('/home/ac833/daophot-options-files/allstar.opt', 'allstar.opt')

	# Open ALLSTAR
	allstar = pexpect.spawn('allstar')

	fout = file(stem+'_allstar_log.txt', 'w')
	allstar.logfile = fout

	allstar.expect("OPT>")
	allstar.sendline("")
	allstar.expect("Input image name:")
	allstar.sendline(stem + '_f' + field + '.fits')
	allstar.expect("File with the PSF")
	allstar.sendline(stem+'_d1_cbcd_dn.psf')
	allstar.expect("Input file")
	allstar.sendline(stem + '_f' + field + '.ap')
	allstar.expect("File for results")
	allstar.sendline(stem + '_f' + field + '.als')
	allstar.expect("Name for subtracted image")
	allstar.sendline(stem + '_f' + field + 's.fits')

	allstar.expect("Good bye")
	allstar.close(force=True)

	# Run DAOPHOT ONE LAST TIME
	daophot = pexpect.spawn('daophot')

	# Set up logfile
	fout = file(stem+'_daophot_log.txt','w')
	daophot.logfile = fout

	daophot.expect("Command:")
	daophot.sendline("off") # offsets to put x and y back in
	daophot.expect("Input file name:")
	daophot.sendline(stem + '_f' + field + '.als') # output file from ALLSTAR
	daophot.expect("Additive offsets ID, DX, DY, DMAG:")
	daophot.sendline("0," + xyoff[0] + "," + xyoff[1] + ",0")
	daophot.expect("Output file name")
	daophot.sendline(stem + '_f' + field + '.mag')

	daophot.expect("Command:")
	daophot.sendline("ex")
	daophot.close(force=True)

	return(0)

def run_allframe(row_of_df, start_dither):

	print "Running ALLFRAME on " + target_name + "     ch: " + wavelength + "     epoch: " + epoch_number

	# Run allstar on every aper phot file for each bcd image
	# for i in range(start_dither, start_dither+5):

	# 	allstar = pexpect.spawn('allstar')

	# 	allstar.expect("OPT>")
	# 	allstar.sendline("")
	# 	allstar.expect("Input image name:")
	# 	allstar.sendline(stem + '_d' + str(i) + '_cbcd_dn.fits')
	# 	allstar.expect("File with the PSF")
	# 	allstar.sendline(stem + '_d' + str(i) + '_cbcd_dn.psf')
	# 	allstar.expect("Input file")
	# 	allstar.sendline(stem + '_d' + str(i) + '_cbcd_dn.ap')
	# 	allstar.expect("File for results")
	# 	allstar.sendline(stem + '_d' + str(i) + '_cbcd_dn.als')
	# 	allstar.expect("Name for subtracted image")
	# 	allstar.sendline(stem + '_d' + str(i) + '_cbcd_dns')

	# 	allstar.expect("Good bye")
	# 	allstar.close(force=True)

	# Change .ap file names in .mch_mast file to .als

	# # Read in the file
	# with open(target_name +'_' + wavelength + '_e'+ epoch_number + '_f'+field+'.mch', 'r') as file:
 #  		filedata = file.read()

	# # Replace the target string
	# filedata = filedata.replace('.ap', '.als')

	# # Write the file out again
	# with open(target_name+'_' + wavelength + '_e' + epoch_number + '_f'+field+'.mch', 'w') as file:
 #  		file.write(filedata)

	# Copy allframe option file if it doesn't already exist in cwd
	shutil.copy('/home/ac833/daophot-options-files/allframe.opt', 'allframe.opt')

	# Open ALLFRAME
	allframe = pexpect.spawn('allframe')

	allframe.expect("OPT>")
	allframe.sendline("")
	allframe.expect("File with list of images:")
	allframe.sendline(stem + '_f' + field + '.mch')
	allframe.expect("File with list of stars")
	allframe.sendline(stem + '_f' + field + '.mag')

	allframe.expect("Good bye")
	allframe.close(force=True)

	return(0)

# Match alf files from allframe for each of the BCDs
def match_alf(row_of_df, start_dither):

	print "Running DAOMATCH on " + target_name + "     ch: " + wavelength + "     epoch: " + epoch_number

	# Change suffix of mch_mast file from als to alf
	# Read in the file
	with open(target_name+'_' + wavelength + '_e'+ epoch_number + '_f'+field+'.mch', 'r') as file:
 		filedata = file.read()

	# Replace the als with alf
	filedata = filedata.replace('.als', '.alf')

	# Write the file out again
	with open(target_name+'_' + wavelength + '_e' + epoch_number + '_f'+field+'.mch', 'w') as file:
  		file.write(filedata)


  	# Now open DAOMATCH to match up the alf files as these are the files you want to obtain magnitudes from!
  	daomatch = pexpect.spawn('daomatch')

  	#fout = file(stem+'_daomatch_log.txt','w')
	#daomatch.logfile = fout

	daomatch.expect("Master input file:")
	daomatch.sendline(stem+'_d'+str(start_dither)+'_cbcd_dn.alf') # Give it the first BCD phot file
	daomatch.expect("Output file name")
	daomatch.sendline(stem+'_f'+field+'.mch')
	daomatch.expect("New output file name")
	daomatch.sendline("")
	daomatch.expect("Next input file:")
	daomatch.sendline(stem+'_d'+str(start_dither + 1)+'_cbcd_dn.alf/') # Give it the second BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline(stem+'_d'+str(start_dither + 2)+'_cbcd_dn.alf/') # Give it the third BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline(stem+'_d'+str(start_dither + 3)+'_cbcd_dn.alf/') # Give it the fourth BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline(stem+'_d'+str(start_dither + 4)+'_cbcd_dn.alf/') # Give it the fifth BCD phot file
	daomatch.expect("Next input file")
	daomatch.sendline("") # exit
	daomatch.expect("Good bye")
	daomatch.close(force=True)	

	# Run DAOMASTER - refine the coordinate transformations from DAOMATCH
	daomaster = pexpect.spawn('daomaster')

	daomaster.expect("File with list of input files:")
	daomaster.sendline(stem+'_f'+field+'.mch')
	daomaster.expect("Minimum number, minimum fraction, enough frames:")
	daomaster.sendline("2, 0.5, 5") # play around with these values
	daomaster.expect("Maximum sigma:")
	daomaster.sendline("0.5") # play around with this value
	daomaster.expect("Your choice:")
	daomaster.sendline("6") # solve for 6 degrees of freedom
	daomaster.expect("Critical match-up radius:")
	daomaster.sendline("7") # play around with this

	for dither in range(start_dither+1,start_dither+5):
		daomaster.expect(stem+'_d'+str(dither)+'_cbcd_dn.alf')
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
	daomaster.sendline("y")
	daomaster.expect("Output file name")
	daomaster.sendline(stem + '_f' + field + '.alf')
	daomaster.expect("A file with the new transformations?")
	daomaster.sendline("y")
	daomaster.expect("Output file name")
	daomaster.sendline(stem+'_f'+field+'.mch')
	daomaster.expect("New output file name")
	daomaster.sendline("")
	daomaster.expect("A file with the transfer table?")
	daomaster.sendline("e") # exits rest of options		

	return(0)



'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
			MAIN PART OF PROGRAM TO LOOP OVER ALL STARS
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Set up data frame from txt file of stars (sys.argv[1]) to do it on
df = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star','Channel','Epoch'])

# Loop over each row in txt file
# Also loop over the two dither combinations
for i in range(0, len(df)):
	for j in [1,6]: 
		
		# INITIAL SETUP 

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

		if j == 1:
			field = '1'
		elif j == 6:
			field = '2'
		else: field = 'Invalid start dither'

    	# Find absolute path of where images are
		home = '/home/ac833/Data/'

		cwd = home + str(galaxy) + '/BCD/' + target_name + '/ch' + str(df['Channel'][i]) + '/e' + str(epoch_number) + '/'
		stem = target_name + '_' + wavelength + '_e' + epoch_number

		# Change directory to where image is
		os.chdir(cwd)

		# Remove any previous runs of this particular script
		if (os.path.isfile(stem+'_montage_log.txt')):
			os.remove(stem+'_montage_log.txt')
		if (os.path.isfile(stem+'_daomatch_log.txt')):
			os.remove(stem+'_daomatch_log.txt')
		if (os.path.isfile(stem+'_daomaster_log.txt')):
			os.remove(stem+'_daomaster_log.txt')
		if (os.path.isfile(stem+'_f'+field+'.mch')):
			os.remove(stem+'_f'+field+'.mch')
		if (os.path.isfile(stem+'_f'+field+'.mch_mast')):
			os.remove(stem+'_f'+field+'.mch_mast')
		if (os.path.isfile(stem+'_f'+field+'.fits')):
			os.remove(stem+'_f'+field+'.fits')

		print "Working on " + stem + "field " + field

		# MAKE MEDIANED IMAGE - OUTPUT OF MEDIAN FUNCTION IS THE X AND Y OFFSETS THAT NEED 
		# TO BE APPLIED IN FIND STARS
		# DEFINITELY WORKS
		xyoff = median(i,j) 

		# FIND STARS ON MEDIANED IMAGE
		find_stars(i,j)

		# RUN ALLFRAME ON MEDIANED IMAGE
		run_allframe(i,j)

		# MATCH ALF FILES
		#match_alf(i,j)

