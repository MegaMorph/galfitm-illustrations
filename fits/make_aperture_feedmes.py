#! /usr/bin/env python
import pyfits
import os
from glob import glob

def make_aperture_feedmes(id, idbase, bandbase='z'):
    # get z-band results
    dfn = '{id}/fit{id}{band}.fits'.format(id=idbase, band=bandbase)
    d = pyfits.getdata(dfn, 'FINAL_BAND')
    # copy and adapt start fit feedmes
    feedmes = glob('{}/fit*galfit'.format(idbase))
    def idswitch(x):
        return x.replace(idbase, id)
    for fn in feedmes:
        outfn = fn.replace(idbase, id)
        copy_feedme_with_fixed_structure(fn, fn.replace(idbase, id), d,
                                         idswitch)


def copy_feedme_with_fixed_structure(fnin, fnout, data, idswitch):
# copy a feedme, fnin, to fnout, replacing structural parameters from
# data and fixing their values
    with open(fnout, 'w') as fout:
        with open(fnin) as f:
            component = 0
            for l in f:
                if l.startswith(' 0)') and 'sky' not in l:
                    component += 1
                elif l.startswith(' 1)'):
                    p = data['COMP%i_XC'%component][0]
                    q = data['COMP%i_YC'%component][0]
                    l = ' 1) {:f} {:f} 0 0\n'.format(p, q)
                elif l.startswith(' 4)'):
                    p = data['COMP%i_Re'%component][0]
                    l = ' 4) {:f} 0\n'.format(p)
                elif l.startswith(' 5)'):
                    p = data['COMP%i_n'%component][0]
                    l = ' 5) {:f} 0\n'.format(p)
                elif l.startswith(' 9)'):
                    p = data['COMP%i_AR'%component][0]
                    l = ' 9) {:f} 0\n'.format(p)
                elif l.startswith('10)'):
                    p = data['COMP%i_PA'%component][0]
                    l = '10) {:f} 0\n'.format(p)
                elif l.startswith('B)'):
                    l = idswitch(l)
                fout.write(l)

if __name__ =='__main__':
    make_aperture_feedmes('A3', 'A2')
    make_aperture_feedmes('Ah3', 'Ah2')
    make_aperture_feedmes('Bh3', 'Bh2')
    make_aperture_feedmes('A6', 'A5')
    make_aperture_feedmes('D3', 'D2')
    make_aperture_feedmes('D6', 'D5')
    make_aperture_feedmes('Dc6', 'Dc5')
    make_aperture_feedmes('Dh6', 'Dh5')
    make_aperture_feedmes('E6', 'E5')
    make_aperture_feedmes('Eh6', 'Eh5')
