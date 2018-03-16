# DAOPHOT

## This repository contains all the scripts required to perform precise photometry of Spitzer Data using DAOPHOT.

The scripts should be performed in this order:
1. file_setup.py - give all your files sensible names and put them in a more logical directory structure
2. convert_to_counts.ipynb - convert all the files in new file system to counts
2. aper_phot.py - find the stars & perform aperture photometry
3. Pick stars & creation of PSF
4. Apply PSF model & subtract stars
5. Apply zero point correction and calibrate aperture to standard system
6. Aperture correction
7. Location correction
8. Pixel phase correction
9. Extract mag and merr for each image for V* to be analysed
10. Plot light curve

These steps will be updated with the script names as and when I make them.

(See full_procedure.pdf for more information)
