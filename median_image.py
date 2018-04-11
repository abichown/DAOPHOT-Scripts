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

# Make medianed image function to call
def median(row_of_df, start_dither):

	print "Making median image of star: " + target_name + "     ch: " + wavelength + "     epoch: " + epoch_number

	# Spawn DAOMATCH - put through the each BCD for the correct dithers
	daomatch = pexpect.spawn('daomatch')

	# Set up log file
	fout = file(stem+'_daomatch_log.txt','w')
	daomatch.logfile = fout

	# Give DAOMATCH all the 5 dithers to be used in making the medianed image
	daomatch.expect("Master input file:")
	daomatch.sendline(stem+'_d'+str(start_dither)+'_cbcd_dn.ap') # Give it the first BCD phot file
	daomatch.expect("Output file name")
	daomatch.sendline(stem+'_f'+field+'.mch')
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

	#print "DAOMATCH has made preliminary coordinate transformations"
	#print "Checking how good they are..."

	# # Open .mch file to check coefficients
	# coeffs = pd.read_csv(stem+'_f'+field+'.mch', header=None, delim_whitespace=True, usecols=[2,3,4,5,6,7], names=['A', 'B', 'C', 'D', 'E', 'F'])

	# if len(coeffs[(coeffs['C'] < 1.01) & (coeffs['C'] > 0.99)]) == 5:
	# 	if len(coeffs[(coeffs['F'] < 1.01) & (coeffs['F'] > 0.99)]) == 5:
	# 		if len(coeffs[(coeffs['D'] < 0.01) & (coeffs['D'] > -0.01)]) == 5:
	# 			if len(coeffs[(coeffs['E'] < 0.01) & (coeffs['E'] > -0.01)]) == 5:
	# 				#print "All coefficients are good"
	# 		    else: print "Coeff E is bad"
	# 	    else: print "Coeff D is bad"
	#     else: print "Coeff F is bad"
 #    else: print "Coeff C is bad"

	# Run DAOMASTER - refine the coordinate transformations from DAOMATCH
	daomaster = pexpect.spawn('daomaster')

	# Set up log file
	fout = file(stem+'_daomaster_log.txt','w')
	daomaster.logfile = fout

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
	daomaster.sendline(stem+'_f'+field+'.mch_mast')
	daomaster.expect("A file with the transfer table?")
	daomaster.sendline("e") # exits rest of options

	# Run MONTAGE2 to actually make medianed image

	montage2 = pexpect.spawn('montage2')

	# Set up log file
	fout = file(stem+'_montage_log.txt','w')
	montage2.logfile = fout

	montage2.expect("File with transformations:")
	montage2.sendline(stem+'_f'+field+'.mch_mast')
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
	montage2.sendline(stem+'_f'+field)

	# Write down X and Y offsets
	log = open(stem+'_montage_log.txt', 'r')
	split = log.read().split()
	offsets = [split[-31], split[-30]]

	# Add back in sky value 
	# THIS IS GETTING THE ERROR OF FLOATING POINT INVALID OPERATION WHEN TRYING TO DO THIS THROUGH IRAF

	return("Complete")

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

	# Then run FIND with this desired threshold

	# Run PHOT

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


		# MAKE MEDIANED IMAGE
		median(i,j)

		# FIND STARS ON MEDIANED IMAGE
		find_stars(i,j)

		# RUN ALLFRAME ON MEDIANED IMAGE

