#!/usr/bin/env python

from glob import glob
import pyfits
import sys, os
import numpy

noisetype = 'realistic'  # like SDSS and UKIDSS
noisetype = 'flat'  # all bands same zeropoint and sky
noisetype = 'simple'  # realistic sky noise, but no noise associated with source counts

shape = (100,100)

bands = ['u', 'g', 'r', 'i', 'z', 'Y', 'J', 'H', 'K']

# averages from SDSS CAS

#Submitted query: select avg(nMgyPerCount_u), avg(nMgyPerCount_g), avg(nMgyPerCount_r), avg(nMgyPerCount_i), avg(nMgyPerCount_z), stdev(nMgyPerCount_u), stdev(nMgyPerCount_g), stdev(nMgyPerCount_r), stdev(nMgyPerCount_i), stdev(nMgyPerCount_z), avg(gain_u), avg(gain_g), avg(gain_r), avg(gain_i), avg(gain_z), avg(sky_u), avg(sky_g), avg(sky_r), avg(sky_i), avg(sky_z) from field

#Result:
zp_sdss = numpy.array([0.0102, 0.0038, 0.0051, 0.0067, 0.0329])  # nMgy/count
gain_sdss = numpy.array([1.71, 3.85, 4.73, 4.94, 4.62]) # e-/count
sky_sdss = numpy.array([1.558, 2.038, 4.778, 8.655, 26.622])  # nMgy/arcsec2
exp_sdss = numpy.array([1.0, 1.0, 1.0, 1.0, 1.0])

sky_sdss *= 0.396**2 * gain_sdss / zp_sdss
zp_sdss = -2.5*numpy.log10(zp_sdss * 1e-9 / gain_sdss) # AB magnitude corresponding to 1 electron in final image

# averages from UKIDSS LAS WSA database

#Submitted query: select filterID, avg(nightZPcat) as nightZPcat, stdev(nightZPcat) as stdev_nightZPcat, avg(totalExpTime) as totalExpTime, avg(gain) as gain, avg(readnoise) as readnoise, avg(skylevel) as skylevel from MultiframeDetector M, lasFrameSets F where M.multiframeID = F.multiframeID and nightZPcat > 1e-5 group by filterID

#Result:
zp_las = numpy.array([22.842, 23.000, 23.237, 22.561])  # magnitude corresponding to 1 count / sec
exptime_las = numpy.array([45.732, 48.157, 45.360, 45.424])  # seconds
gain_las = numpy.array([4.5, 4.5, 4.5, 4.5])  # e- / count
readnoise_las = numpy.array([25.0, 25.0, 25.0, 25.0])
sky_las = numpy.array([466.873, 580.646, 6822.956, 5991.750])  # counts/pixel
exp_las = numpy.array([2.0, 2.0, 4.0, 4.0])
AB_Vega = numpy.array([0.634, 0.938, 1.379, 1.900])

zp_las += 2.5*numpy.log10(exptime_las*gain_las) + AB_Vega # AB magnitude corresponding to 1 electron in final image
sky_las *= gain_las # electrons/pixel

zp_realistic = numpy.concatenate((zp_sdss, zp_las))
sky_realistic = numpy.concatenate((sky_sdss, sky_las))
exp_realistic = numpy.concatenate((exp_sdss, exp_las))

zp_flat = numpy.array([29.0]*9)
sky_flat = numpy.array([1000.0]*9)
exp_flat = numpy.array([1.0]*9)

# Sizes of GAMA galaxies
# At r ~ 17, Re ~ 3 arcsec ~ 9 pixels
# At r ~ 15, Re ~ 6 arcsec ~ 18 pixels

# Keep r ~ 15, and set Re_bulge = 12 pixels, Re_disk = 24 pixels.

def make_images(model='A', brighten=0, bandsel=['u', 'g', 'r', 'i', 'z', 'Y', 'J', 'H', 'K']):
    if noisetype == 'realistic':
        zp = zp_realistic
        sky = sky_realistic
        exp = exp_realistic
    else:
        zp = zp_flat
        sky = sky_flat
        exp = exp_flat
    print 'Using zeropoints:', zp
    print 'Using sky values:', sky
    gals = glob('model%s.galfit'%model)
    for g in gals:
        os.system('nice galfit %s > %s.out'%(g,g))
        imgname = g.replace('.galfit', '')
        img = pyfits.open(imgname+'.fits')
        for j, b in enumerate(bands):
            if b in bandsel:
                ext = img['MODEL_'+b]
                print g, b, j, ext.name
                ext.data *= 10**(0.4*(zp[j]-25+brighten))
                if noisetype == 'simple':
                    sigma = numpy.sqrt(sky[j])
                else:
                    sigma = numpy.sqrt((ext.data+sky[j])/exp[j])
                ext.data += numpy.random.normal(0.0, 1.0, sigma.shape) * sigma
                pyfits.writeto(imgname+'_%i%s_%s%i_sigma.fits'%(j+1, b, noisetype[0], brighten), sigma, clobber=True)
                pyfits.writeto(imgname+'_%i%s_%s%i.fits'%(j+1, b, noisetype[0], brighten), ext.data, clobber=True)


if __name__ =='__main__':
    make_images('A', 0)
    make_images('A', 2)
    make_images('A', -2, ['H'])
    make_images('B', 0)
    make_images('B', 2)
    make_images('B', -2, ['H'])
    for x in 'abcdefghi':
        make_images('C'+x, 0, ['r'])
        make_images('C'+x, 2, ['r'])
    make_images('D', 0)
    make_images('D', 2)
    make_images('E', 0)
    make_images('E', 2)
    make_images('NA', 0)
    make_images('NA', 2)
    make_images('NB', 0)
    make_images('NB', 2)
    make_images('NC', 0)
    make_images('NC', 2)
