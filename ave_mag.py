'''
Purpose: Determines the average magnitude of each star across the 5 dithers for each 
field (dependent on whether it was actually in the frame). Also determines average
error. Writes out these new values to a file with extension .ave
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

# Import modules
import pandas as pd
from math import sqrt, log10

