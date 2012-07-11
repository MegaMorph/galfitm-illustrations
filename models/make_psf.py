import gaussian  # SPB's own function

# produce a symmetrical 2d gaussian centred in stamp
# with sigma = 1 pixel, hence FWHM = 2.355 pixels.

d = gaussian.gaussian2d((25,25), 1.0, 12.0, 12.0, 1.0, 1.0)
pyfits.writeto('psf.fits', d)
