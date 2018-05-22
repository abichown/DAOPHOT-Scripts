# DAOPHOT

## This repository contains all the scripts required to perform precise photometry of Spitzer Data using DAOPHOT.

The scripts should be performed in this order:
1. file_setup.py - give all your files sensible names and put them in a more logical directory structure
2. convert_to_counts.ipynb - convert all the files in new file system to counts
2. aper_phot.py - find the stars & perform aperture photometry.
3. psf_creation.py - pick stars & creation of PSF for one dither. Then copies this PSF model to the other nine dithers.
4. median_image.py - make the median image, find the stars in it and use ALLFRAME to perform the photometry.
5. corrections.py - apply std aperture correction, zero point calibration, aperture correction and location correction.
6. convert_alf.py - run the .alf_apc and .alf files through here to create files which are in the correct format to be put through DAOMATCH and DAOMASTER
6. ave_mag - take alf_all file from step 5 and calculate average mag, average error and std dev for each star at each epoch. Output is .ave file
7. get_mag.py - take .ave file created in step 6 for each epoch and find mag, error and std dev of V* and write to file. Saved in Magnitudes directory.
8. plot_lc.py - plot light curve

These steps will be updated with the script names as and when I make them.

(See full_procedure.pdf for more information)

Additional scripts:
- single_image_phot.py - performs photometry on a single image 

