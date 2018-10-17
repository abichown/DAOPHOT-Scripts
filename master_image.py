'''
Purpose: Create a medianed image from all the 120 (24 x 5) or 60 (12 x 5) frames for a field for a channel.
Then obtain the master star list which gets put through ALLFRAME in the next script.
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

	num_images = num_epochs * 5 # for a 5 dither pattern like we have for CHP data

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

		if start_dither == 1:
			field = '1'
		else: field = '2'

		# Copy all FITS images, psfs and phot files to temp folder
		for epoch in range(1,num_epochs+1):

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			cwd = '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'

			for dither in range(start_dither,start_dither+5):

				dither = str(dither)

				shutil.copyfile(cwd + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.ap', temp + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.ap') #.ap
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

		if field == '1':
			daomatch.sendline(field1_files[0]) # this is epoch 1 dither 1 file
			daomatch.expect("Output file name")
			daomatch.sendline(star_name + '_' + wavelength + '_f1.mch')
		else: 
			daomatch.sendline(field2_files[0]) # this is epoch 1 dither 1 file
			daomatch.expect("Output file name")
			daomatch.sendline(star_name + '_' + wavelength + '_f2.mch')		


		for j in range(1,num_images):

			if field == '1':
				print field1_files[j]
			else:
				print field2_files[j]


			daomatch.expect("Next input file")

			if field == '1':
				daomatch.sendline(field1_files[j]+'/')
			else:
				daomatch.sendline(field2_files[j]+'/')


		daomatch.expect("Next input file")
		daomatch.sendline("") # exit


		# Use DAOMASTER to refine the transformations
		daomaster = pexpect.spawn('daomaster')

		fout = file('daomaster_log.txt','w')
		daomaster.logfile = fout

		print "Running DAOMASTER"

		daomaster.expect("File with list of input files:")
		daomaster.sendline(star_name + '_' + wavelength + '_f' + field + '.mch')
		daomaster.expect("Minimum number, minimum fraction, enough frames:")
		daomaster.sendline("1, 0.5, " + str(num_images)) # play around with these values
		daomaster.expect("Maximum sigma:")
		daomaster.sendline("99") # play around with this value
		daomaster.expect("Your choice:")
		daomaster.sendline("6") # solve for 6 degrees of freedom
		daomaster.expect("Critical match-up radius:")
		daomaster.sendline("7") 

		for j in range(1,num_images):

			if field == '1':
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
		daomaster.sendline("n")
		daomaster.expect("A file with the new transformations?")
		daomaster.sendline("y")
		daomaster.expect("Output file name")
		daomaster.sendline(star_name + '_' + wavelength + '_f' + field + '.mch')
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
		montage2.sendline(star_name + '_' + wavelength + '_f' + field + '.mch')
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
		montage2.sendline(star_name + '_' + wavelength + '_f' + field + '.fits')
		montage2.expect("Good bye")
		montage2.close(force=True)

		# Write down X and Y offsets
		log = open('montage_log.txt', 'r')
		lines = log.readlines()

		offsets = []
		
		for line in lines:
			if "Offsets" in line:

				offsets.append(line.split(' ')[-3])
				offsets.append(line.split(' ')[-2])


		print offsets

		# Use DAOPHOT to create star list
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
		daophot.sendline("at " + star_name + '_' + wavelength + '_f' + field +'.fits')
		daophot.expect("Command:")
		daophot.sendline("opt")
		daophot.expect("File with parameters")
		daophot.sendline("")
		daophot.expect("OPT>")
		daophot.sendline("th=20") 
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

		daophot.expect("Command:")
		daophot.sendline("ph")
		daophot.expect("File with aperture radii")
		daophot.sendline("")
		daophot.expect("PHO>")
		daophot.sendline("")
		daophot.expect("Input position file")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.coo')
		daophot.expect("Output file")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.ap')

		daophot.expect("Command")
		daophot.sendline("ex")

		daophot.close(force=True)

		# Open ALLSTAR
		allstar = pexpect.spawn('allstar')

		fout = file('allstar_log.txt', 'w')
		allstar.logfile = fout

		allstar.expect("OPT>")
		allstar.sendline("")
		allstar.expect("Input image name:")
		allstar.sendline(star_name + '_' + wavelength + '_f'+ field + '.fits')
		allstar.expect("File with the PSF")

		if field == '1':
			allstar.sendline(star_name + '_' + wavelength + '_e01_d1_cbcd_dn.psf')
		else:
			allstar.sendline(star_name + '_' + wavelength + '_e01_d6_cbcd_dn.psf')

		allstar.expect("Input file")
		allstar.sendline(star_name + '_' + wavelength + '_f' + field + '.ap')
		allstar.expect("File for results")
		allstar.sendline(star_name + '_' + wavelength + '_f' + field + '.als')
		allstar.expect("Name for subtracted image")
		allstar.sendline(star_name + '_' + wavelength + '_f' + field + '_dns.fits')

		allstar.expect("Good bye")
		allstar.close(force=True)

		# Run DAOPHOT to add offsets back in
		daophot = pexpect.spawn('daophot')

		daophot.expect("Command:")
		daophot.sendline("off") # offsets to put x and y back in
		daophot.expect("Input file name:")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.als')
		daophot.expect("Additive offsets ID, DX, DY, DMAG:")
		daophot.sendline("0," + offsets[0] + "," + offsets[1] + ",0")
		daophot.expect("Output file name")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.mag')

		daophot.expect("Command:")
		daophot.sendline("ex")
		daophot.close(force=True)	

		# Change .ap file names in .mch file to .als
		f = open(star_name +'_' + wavelength + '_f'+field+'.mch', 'r')
		filedata = f.read()
		f.close()

		newdata = filedata.replace(".ap",".als")

		f = open(star_name +'_' + wavelength + '_f'+field+'.mch','w')
		f.write(newdata)
		f.close()	

		# Now copy the master image and the master star list to each of the epochs 
		for epoch in range(1,num_epochs+1):

			# Get epoch in correct format
			if epoch < 10:
				epoch = '0' + str(epoch)
			else: epoch = str(epoch)

			shutil.copy(star_name + '_' + wavelength + '_f' + field + '.fits', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'+star_name + '_' + wavelength + '_f' + field + '_master.fits')
			shutil.copy(star_name + '_' + wavelength + '_f' + field + '.mch', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'+star_name + '_' + wavelength + '_f' + field + '_master.mch')
			shutil.copy(star_name + '_' + wavelength + '_f' + field + '.mag', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'+star_name + '_' + wavelength + '_f' + field + '_master.mag')

	# Delete temp folder
	shutil.rmtree(temp)