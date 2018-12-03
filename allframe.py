'''
Purpose: Run master star list through ALLFRAME with individual images.
Obtain PSF magnitudes for all stars in master list.
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

# Import modules
import pandas as pd
import os
import pexpect
import shutil


# Get list of fields to make images for i.e. count unique fields in test_star.txt
df = pd.read_csv('/home/ac833/DAOPHOT-Scripts/star_list.txt', header=None, delim_whitespace=True, names=['Galaxy', 'Star', 'Period', 'RA', 'Dec', 'Channel'])

# Don't need this now since won't have any duplicate lines
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

	print galaxy, star_name, wavelength, num_epochs

	for epoch in range(1,num_epochs+1):

		# Get epoch in correct format
		if epoch < 10:
			epoch = '0' + str(epoch)
		else: epoch = str(epoch)

		for start_dither in [6]: #[1,6]

			if start_dither == 1:
				field = '1'
			else: field = '2'

			# Go to folder with data
			directory = '/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e' + epoch +'/'
			os.chdir(directory)

			if epoch != '01':
				shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e01/' + star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.ap', star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.ap')
				shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e01/' + star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.als', star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.als')
				shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e01/' + star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.psf', star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.psf')
				shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e01/' + star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.fits', star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.fits')

			# Open mch file and delete rows not relevant to current epoch
			mch_file = star_name + '_' + wavelength + '_f' + field + '_master.mch'
			mch_df = pd.read_csv(mch_file, delim_whitespace=True, header=None, names=['Filename', 'Apostrophe', 'A', 'B', 'C', 'D', 'E', 'F', 'Mag_offset', 'Scatter'])

			indices_to_remove = []

			for index, row in mch_df.iterrows():

				if index != 0:

					# this string is what we want to look for in all the filenames of the mch file
					epoch_string = 'e' + epoch

					# get list of indices that contain files we want to keep
					if epoch_string not in row['Filename']:
						indices_to_remove.append(index)


			# Drop rows that are not relevant
			mch_df.drop(indices_to_remove, inplace=True)

			# Write out to new mch file
			mch_df.to_csv(mch_file, header=None, sep=' ', index=False)

			# Run ALLFRAME

			# Copy allframe option file if it doesn't already exist in cwd
			shutil.copy('/home/ac833/daophot-options-files/allframe.opt', 'allframe.opt')

			# Copy the master PSF to individual dither names so that ALLFRAME can use them
			for dither in range(start_dither, start_dither+5):
				shutil.copy(star_name + '_' + wavelength + '_f' + field + '_master.psf', star_name + '_' + wavelength + '_' + epoch_string + '_d' + str(dither) + '_cbcd_dn.psf')

			# Open ALLFRAME
			allframe = pexpect.spawn('allframe')

			fout = file('allframe_log.txt', 'w')
			allframe.logfile = fout

			allframe.expect("OPT>")
			allframe.sendline("")
			allframe.expect("File with list of images:")
			allframe.sendline(mch_file)
			allframe.expect("File with list of stars")
			allframe.sendline(star_name + '_' + wavelength + '_f' + field + '_master.mag')

			allframe.expect("Good bye")
			allframe.close(force=True)

			# Match .alf files with DAOMATCH

			# First want to replace .als with .alf in mch file
			f = open(star_name +'_' + wavelength + '_f'+field+'_master.mch', 'r')
			filedata = f.read()
			f.close()

			newdata = filedata.replace(".als",".alf")

			f = open(star_name +'_' + wavelength + '_f'+field+'_master.mch','w')
			f.write(newdata)
			f.close()

			if epoch != '01':
				shutil.copy('/home/ac833/Data/' + galaxy + '/BCD/' + star_name + '/ch' + channel + '/e01/' + star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.alf', star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.alf')
	
			

			"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
							IF YOU WANT A FILE WITH ALL THE INSTRUMENTAL MAGNITUDES IN THEN UNCOMMENT THIS
			"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

			# daomaster = pexpect.spawn('daomaster')

			# fout = file('daomaster_log.txt', 'w')
			# daomaster.logfile = fout

			# daomaster.expect("File with list of input files:")
			# daomaster.sendline(star_name + '_' + wavelength + '_f' + field + '_master.mch')
			# daomaster.expect("Minimum number, minimum fraction, enough frames:")
			# daomaster.sendline("1, 0.5, 5") # play around with these values
			# daomaster.expect("Maximum sigma:")
			# daomaster.sendline("0.5") # play around with this value
			# daomaster.expect("Your choice:")
			# daomaster.sendline("6") # solve for 6 degrees of freedom
			# daomaster.expect("Critical match-up radius:")
			# daomaster.sendline("7") # play around with this

			# if epoch == '01':

			# 	for dither in range(start_dither+1,start_dither+5):
			# 		#daomaster.expect(star_name + '_' + wavelength + '_e' + epoch + '_d' + str(dither) + '_cbcd_dn.alf')
			# 		daomaster.sendline("")

			# else:

			# 	#daomaster.expect(star_name + '_' + wavelength + '_e01_d' + str(start_dither) + '_cbcd_dn.alf')
			# 	#daomaster.sendline("")

			# 	for dither in range(start_dither,start_dither+5):
			# 		#daomaster.expect('.alf')
			# 		#daomaster.expect(star_name + '_' + wavelength + '_e' + epoch + '_d' + str(dither) + '_cbcd_dn.alf')
			# 		daomaster.sendline("")

			# # Repeat with decreasing match up size
			# for match_up in range(7,-1,-1):
			# 	daomaster.expect("New match-up radius")
			# 	daomaster.sendline(str(match_up))

			# # Options for different output files - only want transformations according to cookbook
			# daomaster.expect("Assign new star IDs?")
			# daomaster.sendline("y") # assign new ids so all frames have same ids
			# daomaster.expect("A file with mean magnitudes and scatter?")
			# daomaster.sendline("n")
			# daomaster.expect("A file with corrected magnitudes and errors?")
			# daomaster.sendline("n")
			# daomaster.expect("A file with raw magnitudes and errors?")
			# daomaster.sendline("y")

			# # Obtain final magnitude file for epoch across all 5 dithers

			# daomaster.expect("Output file name")
			# daomaster.sendline(star_name + '_' + wavelength + '_f' + field + '.final_mag')
			# daomaster.expect("A file with the new transformations?")
			# daomaster.sendline("y")
			# daomaster.expect("Output file name")
			# daomaster.sendline(star_name + '_' + wavelength + '_f' + field + '_master.mch')
			# daomaster.expect("New output file name")
			# daomaster.sendline("")
			# daomaster.expect("A file with the transfer table?")
			# daomaster.sendline("e") # exits rest of options	

			# daomaster.close(force=True)	



