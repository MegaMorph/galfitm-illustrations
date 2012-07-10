#!/usr/bin/env python

# Note that this should be used with original GALFIT.

from glob import glob
import pyfits
import sys, os
import numpy

shape = (100,100)

bands = ['u', 'g', 'r', 'i', 'z', 'Y', 'J', 'H', 'K']
zp = numpy.array([17.021,15.909,15.041,14.549,14.262,14.146,14.034,13.793,13.603])

def make_images(noiselevel=5,
                bandsel=['u', 'g', 'r', 'i', 'z', 'Y', 'J', 'H', 'K']):
    
    noisebands = 10**(-0.4*(zp-15.0)) * noiselevel

    noise = []
    for n in noisebands:
        noise.append(numpy.random.normal(0.0, n, shape))

    gals = glob('*.galfit')

    for g in gals:
        os.system('galfit %s'%g)
        imgname = g.replace('.galfit', '')
        img = pyfits.open(imgname+'.fits')
        for j, b in enumerate(bands):
            if b in bandsel:
                ext = img['MODEL_'+b]
                print ext.name, j, noisebands[j]
                ext.data += noise[j]
                ext.writeto(imgname+'_%s_n%i.fits'%(b, noiselevel), clobber=True)

if __name__ =='__main__':
    make_images(5)
    make_images(50, ['H'])
