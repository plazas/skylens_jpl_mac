#!/usr/bin/python
import astropy.io.fits as pyfits
import aplpy
import numpy as np

#cutout_r=pyfits.open("./testimage_noisy_WFIRST_F184.fits")[0].data
#cutout_g=pyfits.open("./testimage_noisy_WFIRST_H158.fits")[0].data
#cutout_b=pyfits.open("./testimage_noisy_WFIRST_J129.fits")[0].data

cutout_r=pyfits.open("./old_fits/testimage_noisy_jwst_F115W.fits")[0].data
cutout_g=pyfits.open("./old_fits/testimage_noisy_jwst_F090W.fits")[0].data
cutout_b=pyfits.open("./old_fits/testimage_noisy_jwst_F070W.fits")[0].data

image_cube = np.zeros((3, cutout_r.shape[0], cutout_r.shape[1]), dtype=np.float32)
image_cube[0,:,:]=cutout_r
image_cube[1,:,:]=cutout_g
image_cube[2,:,:]=cutout_b

pyfits.writeto('new.fits', image_cube, clobber=True)

aplpy.make_rgb_image('new.fits',"color_JWST.png")
