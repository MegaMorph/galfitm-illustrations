#!/usr/bin/env python

from glob import glob
import os
import numpy
import pyfits
from matplotlib import pyplot

bands = ['u', 'g', 'r', 'i', 'z', 'Y', 'J', 'H', 'K']
w = numpy.array([3543,4770,6231,7625,9134,10305,12483,16313,22010], numpy.float)

marker = ['o', '^', 's', 'D', 'x', '+', '*']

def plot_all():
    plot(('A1', 'A2'), '1')
    plot(('Ah1', 'Ah2'), '2')  # and fit 3
    plot(('Bh1', 'Bh2'), '3')  # and fit 3
    plot(('A1', 'A1c'), '4')  # add additional wavelength scale
    plot(('Ah1', 'Ah1c'), '4h')  # add additional wavelength scale
    plot(('A1', 'A1a', 'A1b'), '5')  # requires plotting of chebyshev functions
    plot(('Ah1', 'Ah1a', 'Ah1b'), '5h')  # requires plotting of chebyshev functions
    # plot(('A1', 'A1c', 'A1d'), '6')  # requires plotting of chebyshev functions
    # plot(('Ah1', 'Ah1c', 'Ah1d'), '6h')  # requires plotting of chebyshev functions
    # illustration 7 requires a different kind of plot
    plot(('D1',), '8')  # and D2 and D3
    
def plot(id=('A1', 'A2'), name='0'):
    nid = len(id)
    res = [fit_results(i) for i in id]
    fig = pyplot.figure(figsize=(5, 15))
    fig.subplots_adjust(bottom=0.05, top=0.94, left=0.2, right=0.95, hspace=0.05)
    ax = make_bands_plot(fig, 511, '$m$', True, False)
    simmag = numpy.array([16.935,15.964,15.0,14.562,14.267,14.183,13.992,13.672,13.547])
    pyplot.plot(w, simmag, '-k')
    plotres(res, nid, 'COMP1_MAG')
    ax = make_bands_plot(fig, 512, '$R_e$', False, False)
    plotres(res, nid, 'COMP1_Re')
    ax = make_bands_plot(fig, 513, '$n$', False, False)
    plotres(res, nid, 'COMP1_n')
    ax = make_bands_plot(fig, 514, '$b/a$', False, False)
    plotres(res, nid, 'COMP1_AR')
    ax = make_bands_plot(fig, 515, '$\\theta$', False, True)
    plotres(res, nid, 'COMP1_PA')
    #fig.show()
    fig.savefig('illustration_%s.pdf'%name)

def plotres(res, nid, field):
    mec = ['black', None] * (1+nid//2)
    mfc = ['white', 'black'] * (1+nid//2)
    color = ['grey', 'grey'] * (1+nid//2)
    ymin, ymax = (1e99, -1e99)
    for i in range(nid):
        x = w + 100 * (1+i//2) * (-1)**i
        pyplot.errorbar(x, res[i][field], res[i][field+'_ERR'], color=color[i])
        pyplot.plot(x, res[i][field], color=color[i], marker=marker[i//2],
                    mec=mec[i], markerfacecolor=mfc[i])

        ymin = min(ymin, (res[i][field]-res[i][field+'_ERR']).min())
        ymax = max(ymax, (res[i][field]+res[i][field+'_ERR']).max())
    yrange = ymax - ymin
    ymin -= 0.05 * yrange
    ymax += 0.05 * yrange
    pyplot.ylim(ymin, ymax)


def fit_results(f):
    fn = 'fits/%s/fit%s.fits'%(f,f)
    if os.path.exists(fn):
        r = pyfits.getdata('fits/%s/fit%s.fits'%(f,f), 'final_band')
    else:
        r = numpy.concatenate([pyfits.getdata('fits/%s/fit%s%s.fits'%(f,f, b), 'final_band')
                              for b in bands])
    return r

def make_bands_plot(fig, subplot=111, ylabel='', top=True, bottom=True):    
    ax1 = fig.add_subplot(subplot)
    ax2 = ax1.twiny()
    ax1.set_ylabel(ylabel)
    ax1.set_xlim(2000, 23000)
    if top:
        ax2.set_xlabel('wavelength, $\AA$')
    else:
        ax2.set_xticklabels([])
    ax2.set_xlim(2000, 23000)
    ax1.set_xticks(w)
    if bottom:
        ax1.set_xticklabels(['$'+i+'$' for i in bands])
    else:
        ax1.set_xticklabels([])
    pyplot.setp(ax1.get_xticklabels(), va='baseline')
    pyplot.setp(ax1.get_xaxis().get_major_ticks(), pad=20.)
    pyplot.setp(ax1.get_yaxis().get_major_ticks(), pad=8.)
    ax2.xaxis.labelpad = 12
    return ax1

if __name__ =='__main__':
    plot_all()
