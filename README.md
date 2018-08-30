# DAOPHOT

## This repository contains all the scripts required to perform precise photometry of Spitzer Data using DAOPHOT.

The scripts should be performed in this order:
1. file_setup.py - give all your files sensible names and put them in a more logical directory structure
2. convert_to_counts.ipynb - convert all the files in new file system to counts
3. aper_phot.py - find the stars & perform aperture photometry. Also picks candidate PSF stars
4. star_selection.py - looks through candidate PSF stars and runs a series of tests to determine the suitability of the star to be a PSF star. Outputs all suitable ones to file.
5. make_psf.py - makes the PSF model using the PSF stars from step 4.
6. run_allframe.py - creates the medianed image, puts medianed image through the pipeline and then runs ALLFRAME. Useful outputs are the '.alf' photometry files
7. corrections.py - apply std aperture correction, zero point calibration, aperture correction and location correction to aperture and PSF instrumental magnitudes to get calibrated magnitudes.
8. convert_alf.py - bookkeeping step as the corrected PSF photometry files (.alf_apc) are not in correct format to be put through DAOMATCH/DAOMASTER. This script just puts it in the correct format - called 'alf_all'.
9. ave_mag - take 'alf_all' file from step 8 and calculate average mag, average error and std dev for each star at each epoch. Output is '.ave' file
10. get_mag.py - take '.ave' file from step 9 for each epoch and find mag, error and std dev of V* and write to file. Saved in Magnitudes directory.
11. plot_lc.py - plot light curve


(See full_procedure.pdf for more information)

Additional scripts:
- single_image_phot.py - performs photometry on a single image 
- median_image.py - old redundant script
- psf_automate.py - old redundant script
- psf_creation.py - old redundant script
- star_selection_quad.py - splits field into quadrants and determines the brighest 2 stars from each quadrant to be the PSF stars i.e. 8 stars in total

