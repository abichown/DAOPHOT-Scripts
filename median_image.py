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


# Make medianed image function to call
def median(row_of_df, start_dither, end_dither):

	# Run DAOMATCH - put through the each BCD for the correct dithers

	# Run DAOMASTER - refine the coordinate transformations from DAOMATCH

	# Run MONTAGE2 to actually make medianed image

	# Write down X and Y offsets

	# Add back in sky value 

	return(0)

# Set up data frame from txt file of stars (sys.argv[1]) to do it on
df = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True, names=['Galaxy', 'Star','Channel','Epoch'])

# Loop over each row in txt file 

# Make the 5 images described above 