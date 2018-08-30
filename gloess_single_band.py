
#!/usr/bin/env/python

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import sys
import gloess_fits as gf
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
import os
os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
matplotlib.rc('text',usetex=True)
from matplotlib import rcParams
from math import sqrt
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Garamond']


df = pd.read_csv('/home/ac833/DAOPHOT-Scripts/star_info.txt', header=None, names=['Galaxy', 'Star', 'Channel', 'Period'], delim_whitespace=True)

for i in range(0, len(df)):

	print i 

	# Get all important info
	galaxy = str(df['Galaxy'][i])
	target_star = str(df['Star'][i])
	channel = str(df['Channel'][i])
	period = str(df['Period'][i])

	# Change channel number to wavelength
	if channel == '1':
		wavelength = '3p6um'
	elif channel == '2':
		wavelength = '4p5um'
	else: wavelength = 'not_defined'

	input_filename = '/home/ac833/GLOESS_mags/' + target_star + '_' + wavelength + '_gloess.txt'	

	print input_filename	


	dir1 = [] # magnitudes
	deir1 = [] # errors
	dmjd = [] # hmjd of observations

	counter = 0

	## Want to know whether the IRAC data is phased or not. 
	## If it is phased, must reduce the uncertainty by another factor of sqrt(N)
	## if phased == 1 then true. if phased == 0, false

	for line in open(input_filename):
		data = line.split()
		# Counter always starts equal to 0 so this will always be the first line
		if counter == 0:
			# cepname (Cepheid name?) is the first word in the first line of the input file	
			cepname = data[0]
		if counter == 1:
			# if counter = 1 then the first word in this line of the input file is the period of the cepheid
			period = float(data[0])
			if period > 0:
				# positive number means the period has been phased
				phased = 1
			else:
				# otherwise its not been phased
				phased = 0
		if counter == 2:
			# if counter = 2, the first word in this line of the input file is the number of lines
			nlines = float(data[0])
		if counter == 3:
			# if counter = 3, this line of data contains the x values in all the 12 bands
			xir1 = float(data[0])
		if counter > 3:
			# if counter>3, then this line of the input consists of many values that get appended to the lists at the beginning of the code
			# I think the d in front of du say is just so that when convert to numpy arrays the u isn't already taken
			dmjd.append(float(data[0]))
			dir1.append(float(data[1]))
			deir1.append(float(data[2]))

		# Counter incrememnts so that each line from the input file is read 
		counter  = counter + 1	
			
	## Read in all the data from the file and filled the arrays. Need to convert these to numpy arrays.

	# Number of lines of data as opposed to the Name, period etc
	number = counter - 4 # Number data lines in the file

	# Convert to numpy arrays
	ir1 = np.array(dir1)

	# Associated errors in the bands
	eir1 = np.array(deir1)


	mjd = np.array(dmjd)


	# Sum the values of each band - without the <50 it would just sum the values but the <50 means tell me how many elements have values less than 50
	nir1 = sum(ir1<50)

	# Phases don't need to be done individually by band - only depends on P - check my method of phasing still gives the same results 
	# This might just be a more general case
	phase = (mjd / period) - np.floor(mjd / period)
	#phase = (mjd - mjd[0])/period
	multiple_phase = np.concatenate((phase,(phase+1.0),(phase+2.0),(phase+3.0),(phase+4.0)))


	# Usage:  fit_one_band(data,err,phases,n,smooth):
	maxvals = []
	minvals = []

	# amax is the max value and amin is the min value 
	# so np.amax(u[u<50]) + 3.0 is the largest value of u such that the value is less than 50 and then add 3.0 (I assume this is an offset for ease in plotting)
	# Append the max and min values to the maxvals and minvals lists
	if nir1 > 0:
		maxvals.append(np.amax(ir1[ir1<50]))
		minvals.append(np.amin(ir1[ir1<50]))


	# Convert the max and min values lists to np arrays
	maxvals = np.array(maxvals)
	minvals = np.array(minvals)

	# Max and min across all bands are the largest values of the max and min arrays
	max = np.max(maxvals)
	min = np.min(minvals)

	print cepname, ' ---- Period =', period, 'days'
	print '------------------------------------------------------'

	# Set up names for output files
	#avname = cepname + '.glo_avs'
	#avsout = open(avname,'w')

	# Set max and min limits
	maxlim = max + 0.5
	minlim = min - 0.5

	# Clear figure
	plt.clf()

	#fig = plt.figure()
	#ax1 = fig.add_subplot(111)
	#mp.figure(figsize=(16.0,10.0))


	# gridspec is a module that specifies the location of the subplot in the figure
	# GridSpec specifies the geometry of the grid that a subplot will be placed. 
	# The number of rows and number of columns of the grid need to be set. Optionally, the subplot layout parameters (e.g., left, right, etc.) can be tuned.
	# So we have 3 rows and 4 columns
	gs = gridspec.GridSpec(3, 4)

	# Surely these won't work as you haven't import matplotlib.pyplot as plt, you've imported it as mp
	ax1 = plt.subplot(gs[:,:])


	# Values in the [] are xmin,xmax,ymin,ymax
	ax1.axis([1,3.5,(maxlim),(minlim)])

	titlestring = cepname + ', P = ' + str(period) + ' days'
	plt.suptitle(titlestring, fontsize=20)

	ax1.set_ylabel('Magnitude')
	ax1.set_xlabel('Phase $\phi$')


	## Fitting and plotting for each band
	if nir1 > 0:
		ir11, ir1x, yir1, yeir1, xphaseir1 = gf.fit_one_band(ir1,eir1,phase,nir1,xir1)
		ax1.plot(ir1x,ir11,'k-')
	 	ax1.plot(xphaseir1,yir1,color='MediumVioletRed',marker='o',ls='None', label='mag')
	## for RRLyrae WISE plots:
	#	ax1.plot(ir1x,ir11+1.,'k-')
	# 	ax1.plot(xphaseir1,yir1+1.,color='Turquoise',marker='o',ls='None', label='W1+1.0')
		aveir1, adevir1, sdevir1, varir1, skewir1, kurtosisir1, ampir1 = gf.moment(ir11[200:300],100)
		if phased == 1:
			factor = sqrt(nir1)
		if phased == 0:
			factor = 1 
		if nir1 > 1:
			#print >> avsout, '<mag> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f} N I1 = {3}'.format(aveir1, sdevir1/factor, ampir1,nir1)
			print  '<mag> = {0:.3f}    std dev = {1:.3f}     amplitude = {2:.3f}' .format(aveir1, sdevir1/factor, ampir1)
		if nir1 == 1:
			#print >> avsout, 'mag = {0:.3f} --- single point'.format(aveir1)
			print  'mag = {0:.3f} --- single point'.format(aveir1)


	handles, labels = ax1.get_legend_handles_labels() 
	# ax1.legend(handles[::-1],labels[::-1],loc=4, numpoints=1,prop={'size':10})

	lc_name = '/home/ac833/Light_Curves/' + galaxy +'/' + target_star + '_' + wavelength +'.png'
	plt.savefig(lc_name)

	# Output all relevant information to a file
	# Galaxy, star, channel, period, mean mag, std_dev, amplitude

	master_mags_file = '/home/ac833/master_stars.txt'

	f = open(master_mags_file, 'a')

	f.write(str(galaxy) + ' ' + str(target_star) + ' ' + str(channel) + ' ' + str(period) + ' ' + str(round(aveir1,2)) + ' ' + str(round(sdevir1/factor,2)) + ' ' + str(round(ampir1,2)) ) 

	f.close()

