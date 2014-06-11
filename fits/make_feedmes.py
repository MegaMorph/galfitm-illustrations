#!/usr/bin/env python

from glob import glob
import os
import re

def make_feedmes():
    # One-time script
    # Used to convert all the fit*.galfit files to fit*.diff
    ids = glob('*/')
    for id in ids:
        os.chdir(id)
        feedmes = glob('fit*diff')
        # output starting models
        for f in feedmes:
            template = r'fit.+(\d)(n|m){0,1}([ugrizYJHK]{0,1})([abcde]{0,1})'
            matchobj = re.match(template, f)
            cmd = matchobj.expand('patch -o \g<0>.galfit ../A\g<1>/'
                                  'fitA\g<1>\g<3>.galfit \g<0>.diff')
            print cmd
            os.system(cmd)
        os.chdir('..')

if __name__ =='__main__':
    make_feedmes()
