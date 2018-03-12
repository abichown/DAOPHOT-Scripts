'''
Purpose: Really a practice of how to control DAOPHOT with Python but also looking to set up
optional parameters ready to do photometry
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

import sys
import pexpect
import shutil

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		INITIAL SETUP
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Info entered on command line to identify the image to work on
target_name = sys.argv[1]
channel = sys.argv[2]
epoch_number = sys.argv[3]
dither_number = sys.argv[4]

# Convert channel number to wavelength
if channel == '1':
	wavelength = '3p6um'
elif channel == '2':
	wavelength = '4p5um'
else: wavelength = 'channel not defined'

# Get epoch number into right format
if len(epoch_number) == 1:
	epoch_number = '0' + epoch_number

# File to work on
filename = target_name + '_' + wavelength + '_e' + epoch_number + '_d' + dither_number + '_cbcd_dn'

print 'Currently working on file:' + filename

# Copy daophot options file to current working directory
if channel == '1':
	shutil.copy('/home/ac833/daophot-options-files/daophot.opt', 'daophot.opt')
if channel == '2':
	shutil.copy('/home/ac833/daophot-options-files/daophot.opt', 'daophot.opt')



'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		RUN DAOPHOT
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

daophot = pexpect.spawn("daophot")

# Attach the image
daophot.expect("Command:")
daophot.sendline("AT " + filename)


# Setup log file

