#!/usr/bin/env python

from glob import glob
import os
import re

def make_feedmes():
    # Used to convert all the fit*.diff files to fit*.galfit
    ids = glob('*/')
    for id in ids:
        feedmes = glob(id+'fit*diff')
        # output starting models
        for f in feedmes:
            template = r'.*fit(.*)(\d)(n|m){0,1}([ugrizYJHK]{0,1})([abcde]{0,1})'
            matchobj = re.match(template, f)
            if matchobj.group(1) != 'A' or matchobj.group(5) != '':
                cmd = matchobj.expand('patch -o \g<0>.galfit A\g<2>/'
                                      'fitA\g<2>\g<4>.galfit \g<0>.diff')
                os.system(cmd)

if __name__ =='__main__':
    make_feedmes()
