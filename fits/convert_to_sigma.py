#!/usr/bin/env python

from glob import glob
import os

def convert_to_sigma():
    # One-time script
    # Used to convert feedmes to use sigma image
    ids = glob('*/')
    for id in ids:
        os.chdir(id)
        feedmes = glob('fit*galfit')
        # output starting models
        for f in feedmes:
            text = file(f).readlines()
            fout = file(f, 'w')
            for l in text:
                ls = l.split()
                if len(ls) > 0 and ls[0] == 'A)':
                    cline = l[:l.find('#')].replace('A)', 'C)').replace('.fits', '_sigma.fits')
                if len(ls) > 0 and ls[0] == 'C)':
                    fout.write(cline+'\n')
                else:
                    fout.write(l)
            fout.close()


if __name__ =='__main__':
    convert_to_sigma()
