from astropy.io import fits
import numpy as np
import os
import glob


os.chdir('/grp/hst/wfc3v/martlin/Myles_GD_files/14393_quad_filters')

flc_files = glob.glob('*flc.fits')
flt_files = glob.glob('*_flt.fits')

for image in flc_files:
	hdu=fits.open(image)
	print(hdu[0].header['filename'],hdu[0].header['targname'],hdu[0].header['filter'],hdu[0].header['date-obs'],
	hdu[0].header['detector'],hdu[0].header['IDCTAB'],hdu[0].header['POSTARG1'],hdu[0].header['POSTARG2'])
	
	
