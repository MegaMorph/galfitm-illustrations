#!/usr/bin/env python

from glob import glob
import os
import re

def make_diffs():
    # One-time script
    # Used to convert all the model*.galfit files to model*.diff
    feedmes = glob('model*.galfit')
    # output starting models
    for f in feedmes:
        if f == 'modelA.galfit':
            continue
        diff = f.replace('.galfit', '.diff')
        cmd = 'diff -e modelA.galfit %s > %s'%(f, diff)
        os.system(cmd)

if __name__ =='__main__':
    make_diffs()
