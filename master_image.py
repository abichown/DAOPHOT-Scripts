'''
Purpose: Create a medianed image from all the 120 (24 x 5) or 60 (12 x 5) frames for a field for a channel.
This image will then be used to obtain the master star list that gets put through ALLFRAME.
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

# Import modules
import pandas as pd
import os
import pexpect
import shutil


# Get list of fields to make images for i.e. count unique fields in test_star.txt
df = pd.read_csv('/home/ac833/DAOPHOT-Scripts/test_star.txt', header=None, delim_whitespace=True, usecols=[0,1,2], names=['Galaxy', 'Star', 'Channel'])
df.drop_duplicates(inplace=True)
df.reset_index(drop=True, inplace=True)

print df

for i in range(0,len(df)):

	galaxy = str(df['Galaxy'][i]) # generalise the 0
	star_name = str(df['Star'][i]) # generalise the 0
	channel = str(df['Channel'][i]) # generalise the 0

	if channel == '1':
		wavelength = '3p6um'
	else: wavelength = '4p5um'

	if galaxy == 'LMC':
		num_epochs = 24
	elif galaxy == 'SMC':
		num_epochs = 12
	else: num_epochs = 0 # will error later on in script

	num_images = num_epochs * 5 # for a 5 dither pattern

	print galaxy, star_name, wavelength, num_epochs

	# Create list of all files to be used to make the image
	field1_files = []
	field2_files = []

	# Populate lists

	# Field 1
	for epoch in range(1,num_epochs+1):
		for dither in range(1,6):

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			dither = str(dither)

			filename = star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.ap'
			field1_files.append(filename)


	# Field 2
	for epoch in range(1,num_epochs+1):
		for dither in range(6,11):

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			dither = str(dither)

			filename = star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.ap'
			field2_files.append(filename)


	temp = '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/temp/'

	# Delete temp folder if it already exists - NEED TO ADD IF IT EXISTS TO DELETE
	#shutil.rmtree(temp)

	# Make temp folder
	os.mkdir(temp)

	for start_dither in [1,6]:

		# Copy all FITS images, psfs and phot files to temp folder
		for epoch in range(1,num_epochs+1):

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			cwd = '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'

			for dither in range(start_dither,start_dither+5):

				dither = str(dither)

				shutil.copyfile(cwd + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.ap', temp + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.ap')
				shutil.copyfile(cwd + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.fits', temp + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.fits')
				shutil.copyfile(cwd + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.psf', temp + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.psf')

		# Change to temp folder 
		os.chdir(temp)

		# Use DAOMATCH to get initial transformations
		daomatch = pexpect.spawn('daomatch', timeout=30)

		fout = file('log.txt','w')
		daomatch.logfile = fout

		# Input epoch 1 first dither of field as master input file
		daomatch.expect("Master input file:")

		if start_dither == 1:
			daomatch.sendline(field1_files[0]) # this is epoch 1 dither 1 file
			daomatch.expect("Output file name")
			daomatch.sendline(star_name + '_' + wavelength + '_f1.mch')
		else: 
			daomatch.sendline(field2_files[0]) # this is epoch 1 dither 1 file
			daomatch.expect("Output file name")
			daomatch.sendline(star_name + '_' + wavelength + '_f2.mch')		


		for j in range(1,120):
			print j

			try:
				daomatch.expect("Next input file")

				if start_dither == 1:
					daomatch.sendline(field1_files[j]+'/')
				else:
					daomatch.sendline(field2_files[j]+'/')
				pass

			except: 
				print "Bad"
				daomatch.expect("Write this transformation?")
				daomatch.sendline("No")


		daomatch.expect("Next input file")
		daomatch.sendline("") # exit


		# Use DAOMASTER to refine the transformations
		daomaster = pexpect.spawn('daomaster')

		fout = file('daomaster_log.txt','w')
		daomaster.logfile = fout

		print "Running DAOMASTER"

		daomaster.expect("File with list of input files:")

		if start_dither == 1:
			daomaster.sendline(star_name + '_' + wavelength + '_f1.mch')
		else: 
			daomaster.sendline(star_name + '_' + wavelength + '_f2.mch')

		daomaster.expect("Minimum number, minimum fraction, enough frames:")
		daomaster.sendline("1, 0.5, 120") # play around with these values
		daomaster.expect("Maximum sigma:")
		daomaster.sendline("99") # play around with this value
		daomaster.expect("Your choice:")
		daomaster.sendline("6") # solve for 6 degrees of freedom
		daomaster.expect("Critical match-up radius:")
		daomaster.sendline("7") 

		for j in range(1,120):

			if start_dither == 1:
				daomaster.expect(field1_files[j])
				daomaster.sendline("")
			else:
				daomaster.expect(field2_files[j])
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

		if start_dither == 1:
			daomaster.sendline(star_name + '_' + wavelength + '_f1.raw')
		else: 
			daomaster.sendline(star_name + '_' + wavelength + '_f2.raw')		

		daomaster.expect("A file with the new transformations?")
		daomaster.sendline("y")
		daomaster.expect("Output file name")

		if start_dither == 1:
			daomaster.sendline(star_name + '_' + wavelength + '_f1.mch')
		else:
			daomaster.sendline(star_name + '_' + wavelength + '_f2.mch')

		daomaster.expect("New output file name")
		daomaster.sendline("")
		daomaster.expect("A file with the transfer table?")
		daomaster.sendline("e") # exits rest of options

		print "Running MONTAGE2"

		# Use MONTAGE2 to actually make image
		montage2 = pexpect.spawn('montage2')

		# Set up log file
		fout = file('montage_log.txt','w')
		montage2.logfile = fout

		montage2.expect("File with transformations:")

		if start_dither == 1:
			montage2.sendline(star_name + '_' + wavelength + '_f1.mch')
		else:
			montage2.sendline(star_name + '_' + wavelength + '_f2.mch')

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

		if start_dither == 1:
			montage2.sendline(star_name + '_' + wavelength + '_f1.fits')
		else:
			montage2.sendline(star_name + '_' + wavelength + '_f2.fits')

		montage2.expect("Good bye")
		montage2.close(force=True)

		# Use DAOPHOT to find optimal threshold (elbow of curve) and then create star list
		shutil.copy('/home/ac833/daophot-options-files/daophot.opt', 'daophot.opt')
		shutil.copy('/home/ac833/daophot-options-files/photo.opt', 'photo.opt')
		shutil.copy('/home/ac833/daophot-options-files/allstar.opt', 'allstar.opt')

		daophot = pexpect.spawn('daophot') 

		print "Running DAOPHOT"

		# Set up logfile
		fout = file('daophot_log.txt','w')
		daophot.logfile = fout

		# Attach medianed image
		daophot.expect("Command:")

		if start_dither == 1:
			daophot.sendline("at " + star_name + '_' + wavelength + '_f1.fits')
		else:
			daophot.sendline("at " + star_name + '_' + wavelength + '_f2.fits')

		# Experiment with threshold to find the value near the 'elbow' of the curve
		daophot.expect("Command:")
		daophot.sendline("opt")
		daophot.expect("File with parameters")
		daophot.sendline("")
		daophot.expect("OPT>")
		daophot.sendline("th=40") # set initial threshold to 2
		daophot.expect("OPT>")
		daophot.sendline("")

		daophot.expect("Command:")
		daophot.sendline("fi")
		daophot.expect("Number of frames averaged, summed:")
		daophot.sendline(str(num_images)+",1")
		daophot.expect("File for positions")
		daophot.sendline("")
		daophot.expect("Are you happy with this?")
		daophot.sendline("y")

		daophot.expect("Command")
		daophot.sendline("ex")

		daophot.close(force=True)

		# Now copy the master image and the master star list to each of the epochs 
		for epoch in range(1,num_epochs+1):

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			if start_dither == 1:
				shutil.copy(star_name + '_' + wavelength + '_f1.fits', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'+star_name + '_' + wavelength + '_f1_master_image.fits')
				shutil.copy(star_name + '_' + wavelength + '_f1.coo', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'+star_name + '_' + wavelength + '_f1_master_star_list.coo')
			else:
				shutil.copy(star_name + '_' + wavelength + '_f2.fits', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'+star_name + '_' + wavelength + '_f2_master_image.fits')
				shutil.copy(star_name + '_' + wavelength + '_f2.coo', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'+star_name + '_' + wavelength + '_f2_master_star_list.coo')


	# Delete temp folder
	#shutil.rmtree(temp)