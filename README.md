# DAOPHOT

## This repository contains all the scripts required to perform precise photometry of Spitzer Data using DAOPHOT.

The scripts should be performed in this order:
1. file_setup.py - give all your files sensible names and put them in a more logical directory structure
2. convert_to_counts.ipynb - convert all the files in new file system to counts
2. aper_phot.py - find the stars & perform aperture photometry.
3. psf_creation.py - pick stars & creation of PSF for one dither. Then copies this PSF model to the other nine dithers.
4. median_image.py - make the median image, find the stars in it and use ALLFRAME to perform the photometry.
5. Apply zero point correction and calibrate aperture to standard system
6. Aperture correction
7. Location correction
8. ave_mag - take alf file and calculate average mag and error for each star at each epoch
9. Extract mag and err for each image for V* to be analysed
10. Plot light curve

These steps will be updated with the script names as and when I make them.

(See full_procedure.pdf for more information)

Additional scripts:
- single_image_phot.py - performs photometry on a single image 

