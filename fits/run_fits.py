#!/usr/bin/env python

from glob import glob
import os, sys
import subprocess
import select
import re

def run_all_fits(idglob='*'):
    poller = select.poll()
    subprocs = {} # map stdout pipe's file descriptor to the Popen object
    ids = glob(idglob+'/')
    commands = []
    for id in ids:
        os.chdir(id)
        feedmes = glob('fit*galfit')
        # output starting models
        for f in feedmes:
            s = f.replace('fit', 'start', 1)
            sfn = s.replace('.galfit', '.fits')
            cmd = 'cd %s; nice galfit -o1 -f %s %s > %s.out; if [ $? -eq 0 ]; then echo %s: success; else echo %s: failure; fi'%(id,sfn,f,s,s,s)
            commands.append(cmd)
        # do proper fits
        for f in feedmes:
            cmd = 'cd %s; nice galfit %s > %s.out; if [ $? -eq 0 ]; then echo %s: success; else echo %s: failure; fi'%(id,f,f,f,f)
            commands.append(cmd)
        os.chdir('..')
    #loop that polls until completion
    while True:
        if len(commands) > 0 and len(subprocs) < 10:
            cmd = commands.pop(0)
            subproc = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            subprocs[subproc.stdout.fileno()] = subproc
            poller.register(subproc.stdout, select.POLLHUP)
        else:
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
    if len(sys.argv) > 1:
        run_all_fits(sys.argv[1])
    else:
        run_all_fits()
