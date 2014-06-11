#!/usr/bin/env python

from glob import glob
import os, sys
import subprocess
import select
import re

def run_all_fits():
    poller = select.poll()
    subprocs = {} # map stdout pipe's file descriptor to the Popen object
    ids = glob('*/')
    for id in ids:
        os.chdir(id)
        feedmes = glob('fit*galfit')
        # output starting models
        for f in feedmes:
            s = f.replace('fit', 'start', 1)
            sfn = s.replace('.galfit', '.fits')
            cmd = 'nice galfit -o1 -f %s %s > %s.out; if [ $? -eq 0 ]; then echo %s: success; else echo %s: failure; fi'%(sfn,f,s,s,s)
            subproc = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            subprocs[subproc.stdout.fileno()] = subproc
            poller.register(subproc.stdout, select.POLLHUP)
        # do proper fits
        for f in feedmes:
            cmd = 'nice galfit %s > %s.out; if [ $? -eq 0 ]; then echo %s: success; else echo %s: failure; fi'%(f,f,f,f)
            subproc = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            subprocs[subproc.stdout.fileno()] = subproc
            poller.register(subproc.stdout, select.POLLHUP)
        os.chdir('..')
    #loop that polls until completion
    while True:
        for fd, flags in poller.poll(5000): # never more than five seconds without a UI update
            done_proc = subprocs[fd]
            del subprocs[fd]
            poller.unregister(fd)
            print '\r\033[K'+done_proc.communicate()[0],
        print '.',
        sys.stdout.flush()
        #don't forget to break when all are done
        if len(subprocs) == 0:
            print '\r\033[KAll complete'
            break


if __name__ =='__main__':
    run_all_fits()
