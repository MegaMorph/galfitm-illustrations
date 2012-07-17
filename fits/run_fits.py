#!/usr/bin/env python

from glob import glob
import os

def run_all_fits():
    ids = glob('*/')
    for id in ids:
        run_fits(id)

def run_fits(id='A1'):
    os.chdir(id)
    feedmes = glob('fit*galfit')
    for f in feedmes:
        os.system('galfit %s > %s.out; if [ $? -eq 0 ]; then echo %s: success; else echo %s: failure; fi'%(f,f,f,f))
    os.chdir('..')

if __name__ =='__main__':
    run_all_fits()
