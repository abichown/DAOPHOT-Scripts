'''
Purpose: Extract magnitude and error of the V*. At the moment, V* always in centre of frame
so program looks for star near (133,122). Also, V* only present in [3.6] field 2 and 
[4.5] field 1 so only looks there.
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

# Import modules
import os
import sys
import pandas as pd
import numpy as np