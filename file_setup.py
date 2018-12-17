'''
Purpose: Sort all 'r...' folders into a better file structure in my home directory
New file structure: Data > Galaxy > BCD/PBCD > Star > Channel > Epoch > dither
Must create Data folder with subdirectories for each Galaxy and BCD and PBCD folders prior to running
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

from astropy.io import fits
import os
import shutil
import glob
import fnmatch

# This is a searching script that looks for all files in the file system called certain things
os.chdir('/home/ac833/Python-Scripts/')
from search_files import unsorted_files

def find_bcd_a_home():

	# Home to all the unorganised 'raw' (for want of a better word) data
	os.chdir('/home/ac833/Downloads/P61004/')
	#os.chdir('/mnt/ac833-XDrive/Physics/ResearchProjects/VScowcroft/EA-PH1166/Raw_Data/LMC/')

	# Search folders for files and write to text file '~/unsorted_files.txt'
	unsorted_files()

	# Open this text file 
	f = open('/home/ac833/unsorted_files.txt', 'r')
	frames = []

	# Write filenames to 'frames' list
	for file in f:
	    frames.append(file.rstrip()) # writes all filenames to the list and also removes the '\n' 

	# Root data directory of where the files will be stored
	root = '/home/ac833/Data/'

	# Places each file in the list in its correct place
	for filename in frames:

	    # Open FITS image
	    image = fits.open(filename)
	    
	    # Get header
	    header = image[0].header
	    
	    # Grab AORLabel and split into galaxy [0], cepheid [1] and epoch [2]
	    galaxy = header['AORLABEL'].split()[0]
	    cepheid = header['AORLABEL'].split()[1]
	    epoch = header['AORLABEL'].split()[2]
	    
	    # Grab Channel number and whether image is BCD or PBCD
	    path = os.path.dirname(filename)
	    channel = path.split('/')[-2]
	    bcd_or_pbcd = path.split('/')[-1]

	    # Grab dither position
	    dither = header['DITHPOS']
	    
	    # BCD or PBCD
	    if bcd_or_pbcd.lower() == 'bcd':
	        bcd_or_pbcd = 'bcd'
	    elif bcd_or_pbcd.lower() == 'pbcd':
	        bcd_or_pbcd = 'pbcd'
	    else:
	        bcd_or_pbcd = 'error'
	        print "Not a BCD or PBCD image"
	    
	    # Create directory for star if it doesn't already exist
	    pathname = root + str(galaxy) + '/' + str(bcd_or_pbcd.upper()) + '/' + str(cepheid)

	    if os.path.exists(pathname) == False:
	        os.mkdir(pathname)
	        
	    # Find out if it is channel 1 or channel 2
	    if channel == 'ch1':
	        channel = 'ch1'
	    elif channel == 'ch2':
	        channel = 'ch2'
	    else:
	        channel = 'error'
	        print "Not Channel 1 or Channel 2"
	    
	    # Make channel directory if it doesn't already exist
	    pathname = root + str(galaxy) + '/' + str(bcd_or_pbcd.upper()) + '/' + str(cepheid) + '/' + str(channel)

	    if os.path.exists(pathname) == False:
	        os.mkdir(pathname)
	        
	    # Find out what epoch
	    if len(epoch) == 7: #hence it is a double-digit number from 10-24
	        epoch = str(epoch[-2:])
	    else: epoch = '0' + str(epoch[-1:]) #so it is between 1-9
	    
	    # Make epoch directory if it doesn't already exist
	    pathname = root + str(galaxy) + '/' + str(bcd_or_pbcd.upper()) + '/' + str(cepheid) + '/' + str(channel) + '/e' + epoch

	    if os.path.exists(pathname) == False:
	        os.mkdir(pathname)
	        
	    # Make dither position run from 1 - 10 not 1-5 then 1-5 again
	    checks = ['_0005_0000', '_0006_0000', '_0007_0000', '_0008_0000', '_0009_0000']
	    counter = 0
	    
	    for word in checks:
	        if filename.find(word) != -1:
	            counter += 1
	            
	    if counter != 0:
	        dither = dither + 5
	    
	    # # Make dither directory if it doesn't already exist
	    # pathname = root + str(galaxy) + '/' + str(bcd_or_pbcd.upper()) + '/' + str(cepheid) + '/' + str(channel) + '/e' + epoch + '/d' + str(dither)

	    # if os.path.exists(pathname) == False:
	    #     os.mkdir(pathname)
	    
	    if channel == 'ch1':
	        ch = '3p6um'
	    else:
	        ch = '4p5um'
	    
	    # Place FITS file here with better name
	    new_name = pathname + '/'+ str(cepheid) + '_' + ch + '_e' + str(epoch) + '_d' + str(dither) + '_cbcd.fits'

	    # Copy file
	    shutil.copyfile(filename, new_name)
	    
	    image.close() 

	return("All BCD images found a home")

def find_pbcd_a_home():

	os.chdir('/mnt/ac833-XDrive/Physics/ResearchProjects/VScowcroft/Scowcroft_Astro/Spitzer/LMC/AORs/')

	# Search folders for files and write to text file '~/unsorted_files.txt'
	unsorted_files()

	# Open this text file 
	f = open('/home/ac833/unsorted_pbcd_files.txt', 'r')
	frames = []

	# Write filenames to 'frames' list
	for file in f:
	    frames.append(file.rstrip()) # writes all filenames to the list and also removes the '\n' 

	# Root data directory of where the files will be stored
	root = '/home/ac833/Data/'

	# Places each file in the list in its correct place
	for filename in frames:
    
	    # Open FITS image
	    image = fits.open(filename)
	    
	    # Get header
	    header = image[0].header
	    
	    # Grab AORLabel and split into galaxy [0], cepheid [1] and epoch [2]
	    galaxy = header['AORLABEL'].split()[0]
	    cepheid = header['AORLABEL'].split()[1]
	    epoch = header['AORLABEL'].split()[2]
	    
	    # Grab Channel number and whether image is BCD or PBCD
	    path = os.path.dirname(filename)
	    channel = path.split('/')[-2]
	    bcd_or_pbcd = path.split('/')[-1]

	    # BCD or PBCD
	    if bcd_or_pbcd.lower() == 'bcd':
	        bcd_or_pbcd = 'bcd'
	    elif bcd_or_pbcd.lower() == 'pbcd':
	        bcd_or_pbcd = 'pbcd'
	    else:
	        bcd_or_pbcd = 'error'
	        print "Not a BCD or PBCD image"
	    
	    # Create directory for star if it doesn't already exist
	    pathname = root + str(galaxy) + '/' + str(bcd_or_pbcd.upper()) + '/' + str(cepheid)

	    if os.path.exists(pathname) == False:
	        os.mkdir(pathname)
	        
	    # Find out if it is channel 1 or channel 2
	    if channel == 'ch1':
	        channel = 'ch1'
	    elif channel == 'ch2':
	        channel = 'ch2'
	    else:
	        channel = 'error'
	        print "Not Channel 1 or Channel 2"
	    
	    # Make channel directory if it doesn't already exist
	    pathname = root + str(galaxy) + '/' + str(bcd_or_pbcd.upper()) + '/' + str(cepheid) + '/' + str(channel)

	    if os.path.exists(pathname) == False:
	        os.mkdir(pathname)
	        
	    # Find out what epoch
	    if len(epoch) == 7: #hence it is a double-digit number from 10-24
	        epoch = str(epoch[-2:])
	    else: epoch = '0' + str(epoch[-1:]) #so it is between 1-9
	    
	    # Make epoch directory if it doesn't already exist
	    pathname = root + str(galaxy) + '/' + str(bcd_or_pbcd.upper()) + '/' + str(cepheid) + '/' + str(channel) + '/e' + epoch + '/'

	    if os.path.exists(pathname) == False:
	        os.mkdir(pathname)
	        
	    filetype = filename[-9:-5]
	    
	    if channel == 'ch1':
	        ch = '3p6um'
	    else:
	        ch = '4p5um'
	    
	    # Place FITS file here with better name
	    new_name = pathname + '/'+ str(cepheid) + '_' + ch + '_' + filetype + '.fits'

	    # Copy file
	    shutil.copyfile(filename, new_name)
	    
	    image.close()

	return("All PBCD images found a home")