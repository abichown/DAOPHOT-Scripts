'''
Purpose: Select candidate PSF stars from star list. 
Then do series of tests to remove stars that wouldn't make good psf stars
The following tests are carried out:
1. Is the star too close to the edge of the frame?
2. Is the star bright enough?
3. Does the star look star-like?
4. Is it just a hot pixel?
The stars that pass all the test are outputted to the lst file to be used in the psf creation process
Written by: Abi Chown A.H.Chown@bath.ac.uk
'''

import 