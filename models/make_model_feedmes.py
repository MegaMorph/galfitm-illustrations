#!/usr/bin/env python

from glob import glob
import os
import re

def make_feedmes():
    # Used to convert all the model*.diff files to model*.galfit
    diffs = glob('model*.diff')
    # output starting models
    for d in diffs:
        feedme = d.replace('.diff', '.galfit')
        cmd = 'patch -o %s modelA.galfit %s'%(feedme, d)
        os.system(cmd)

if __name__ =='__main__':
    make_feedmes()
