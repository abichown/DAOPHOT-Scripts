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
import astropy.io.fits as fits


# Get list of fields to make master images for 
df = pd.read_csv('/home/ac833/DAOPHOT-Scripts/star_list.txt', header=None, delim_whitespace=True, usecols=[0,1,2,3], names=['Galaxy', 'Star', 'Period', 'Channel'])

# No longer need these as won't have duplicates
#df.drop_duplicates(inplace=True)
#df.reset_index(drop=True, inplace=True)

for i in range(0,len(df)):

	galaxy = str(df['Galaxy'][i]) 
	star_name = str(df['Star'][i]) 
	channel = str(df['Channel'][i]) 

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
	if (os.path.isdir(temp)):
		shutil.rmtree(temp)

	# Make temp folder
	os.mkdir(temp)

	for start_dither in [6]: # [1,6] just checking field 2 works right now 

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
				#shutil.copyfile(cwd + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.psf', temp + star_name + '_' + wavelength + '_e' + epoch + '_d' + dither + '_cbcd_dn.psf')

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

		# Now choose brightest 20 stars in image to be candidate PSF stars
		daophot.expect("Command:")
		daophot.sendline("pi")
		daophot.expect("Input file name")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.ap')
		daophot.expect("Desired number of stars, faintest magnitude:")
		daophot.sendline("20,99")
		daophot.expect("Output file name")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.lst')
		daophot.expect("Command:")
		daophot.sendline("ex")

		daophot.close(force=True)

		# Now run these candidate PSF stars through the series of tests to get rid of bad stars

		# Read in FITS image
		hdulist = fits.open(star_name + '_' + wavelength + '_f' + field +'.fits')

		# Access the primary header-data unit (HDU)
		hdu = hdulist[0]
		data = hdu.data

		df2 = pd.read_csv(star_name + '_' + wavelength + '_f' + field + '.lst', delim_whitespace=True, skiprows=3, header=None, names=['ID', 'X', 'Y', 'Mag', 'Error'], index_col=0)

		print df2

		# Carry out all the tests on each star in the df
		for index, row in df2.iterrows():

			execute = 1

			# TEST 1 : TOO CLOSE TO EDGE OF FRAME

			# If X < 64 or X >330, drop row
			if row['X'] < 64 or row['X'] > 330:
				df2.drop(index, inplace=True)
				print "Deleting star %d because it is too close to edge of frame" % index
				execute = 0 # don't need to carry out rest of tests

			if execute == 1:

				# If Y < 56 or Y > 369, drop row
				if row['Y'] < 56 or row['Y'] > 369:
					df2.drop(index, inplace=True)
					print "Deleting star %d because it is too close to edge of frame" % index
					execute = 0 # don't need to carry out rest of tests

			# TEST 2 : NOT BRIGHT ENOUGH

			# Get x and y coords of the star in question
			x_coord = int(round(row['X'] - 1)) # zero-indexed in data and must be rounded to nearest integer
			y_coord = int(round(row['Y'] - 1)) # zero-indexed in data and must be rounded to nearest integer

			if execute == 1:
				if data[y_coord, x_coord] < 150:
					df2.drop(index, inplace=True)
					print "Deleting star %d because it is not bright enough" % index
					execute = 0 # don't need to carry out rest of tests

			# TEST 3 : NOT STAR-LIKE

			# Value of pixel at centre
			centre = data[y_coord, x_coord]

			# Set up of rings for the nearest neighbours
			ring1 = [data[(y_coord, x_coord - 1)], data[(y_coord, x_coord + 1)], data[(y_coord - 1, x_coord)], data[(y_coord + 1, x_coord)]]

			ring2 = [data[(y_coord, x_coord - 2)], data[(y_coord + 1, x_coord - 1)], data[(y_coord + 2, x_coord)], data[(y_coord + 1, x_coord + 1)], data[(y_coord, x_coord + 2)], data[(y_coord - 1, x_coord + 1)], data[(y_coord - 2, x_coord)], data[(y_coord - 1, x_coord - 1)]]

			ring3 = [data[(y_coord, x_coord - 3)], data[(y_coord + 1, x_coord - 2)], data[(y_coord + 2, x_coord - 1)], data[(y_coord + 3, x_coord)], data[(y_coord + 2, x_coord + 1)], data[(y_coord + 1, x_coord + 2)], data[(y_coord, x_coord + 3)], data[(y_coord - 1, x_coord + 2)], data[(y_coord - 2, x_coord + 1)], data[(y_coord - 3, x_coord)], data[(y_coord - 2, x_coord - 1)], data[(y_coord - 1, x_coord - 2)]]

			# Get average values of the rings
			ave_1 = sum(ring1)/len(ring1)
			ave_2 = sum(ring2)/len(ring2)
			ave_3 = sum(ring3)/len(ring3)

			if execute == 1:

				# Is centre > ave1 > ave2 > ave3?
				if ave_1 > centre or ave_2 > centre or ave_3 > centre or ave_2 > ave_1 or ave_3 > ave_1 or ave_3 > ave_2: 
					df2.drop(index, inplace=True)
					print "Deleting star %d because it isn't star-like" % index
					execute = 0 # don't need to do rest of tests

			# TEST 4 : ELIMINATE HOT PIXELS

			if execute == 1:

				# Check ratio between centre and ring1 is > 10%
				if ave_1 / centre < 0.1:
					df2.drop(index, inplace=True)
					print "Deleting star %d because it is just a hot pixel" % index
					execute = 0


			# # TEST 5 : CLOSE NEIGHBOURS

			# if execute == 1:

			# 	# Initially set neighbour = 0 i.e. assume no neighbours
			# 	neighbour = 0

			# 	# Write out coordinates of star
			# 	coords = []

			# 	for i in range(-6,7):
			# 		for j in range(-6,7):
			# 			coords.append((y_coord+j, x_coord+i))

			# 	# Check 20 x 20 grid for stars nearby
			# 	for i in range(-10,11):
			# 		for j in range(-10,11):

			# 			value = 100

			# 			# Check the value isn't in the star's area
			# 			if (y_coord+j, x_coord+i) not in coords:
			# 				if y_coord+j < 256 and y_coord >= 0:
			# 					if x_coord+i < 256 and x_coord+i >= 0:

			# 						# Check for neighbours
			# 						if data[(y_coord+j, x_coord+i)] > value:
			# 							detection = data[(y_coord+j, x_coord+i)]
			# 							neighbour = 1

			# 	if neighbour == 1:
			# 		print "Deleting star %d because there is a pixel value of %f nearby" % (index, detection)
			# 		df.drop(index, inplace=True)  

			# Write out final list of stars to the lst file in the correct format

			# Get header of lst file
			f = open(star_name + '_' + wavelength + '_f' + field + '.lst', 'r')
			header = f.read().splitlines()[0:3]
			f.close()

			# Now overwrite this file
			f = open(star_name + '_' + wavelength + '_f' + field + '.lst', 'w')
			f.writelines(header[0] + '\n' + header[1] + '\n' + header[2] + '\n')

			# Send stars to lst file
			df2.to_csv(f, sep=' ', mode='a', header=None) #.iloc[0:6]

			f.close()

		daophot = pexpect.spawn('daophot')

		daophot.expect("Command:")
		daophot.sendline("at " + star_name + '_' + wavelength + '_f' + field +'.fits')

		# Now create PSF model to be used in allstar in next step of this script
		daophot.expect("Command:")
		daophot.sendline("psf")
		daophot.expect("File with aperture results")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.ap')
		daophot.expect("File with PSF stars")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.lst')
		daophot.expect("File for the PSF")
		daophot.sendline(star_name + '_' + wavelength + '_f' + field + '.psf')

		daophot.expect("Command:")
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
		allstar.sendline(star_name + '_' + wavelength + '_f' + field + '.psf')

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
			shutil.copy(star_name + '_' + wavelength + '_f' + field + '.psf', '/home/ac833/Data/'+galaxy+'/BCD/'+star_name+'/ch'+channel+'/e'+epoch+'/'+star_name + '_' + wavelength + '_f' + field + '_master.psf')

	# Delete temp folder
	#shutil.rmtree(temp)