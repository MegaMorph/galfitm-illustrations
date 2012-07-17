#!/usr/bin/env python

from glob import glob
import os
import numpy
import pyfits
import matplotlib
from matplotlib import pyplot
from numpy.polynomial.chebyshev import Chebyshev

matplotlib.rcParams.update({'font.size': 16})

bands = ['u', 'g', 'r', 'i', 'z', 'Y', 'J', 'H', 'K']

w = numpy.array([3543,4770,6231,7625,9134,10305,12483,16313,22010], numpy.float)

xlim = (2000, 23000)

simmag = numpy.array([16.935,15.964,15.0,14.562,14.267,14.183,13.992,13.672,13.547])

marker = ['o', '^', 's', 'D', 'x', '+', '*']

ylim_std = {'MAG': (17.51, 12.49), 'Re': (0.01, 11.99), 'n': (0.01, 5.99),
            'AR': (0.41, 0.79), 'PA': (30.01, 59.99)}

varlist_std = ('MAG', 'Re', 'n', 'AR', 'PA')

labels = {'MAG': '$m$', 'Re': '$R_e$', 'n': '$n$', 'AR': '$b/a$', 'PA': '$\\theta$'}

wlfuncs = {'A1c': numpy.log10, 'Ah1c': numpy.log10}

def plot_all():
    plot(('A2', 'A1'), 1, '1', 'True')
    plot(('Ah2', 'Ah1'), 1, '2', 'True')  # and fit 3
    plot(('Bh2', 'Bh1'), 1, '3', 'True')  # and fit 3
    plot(('A1c', 'A1'), 1, '4', 'True')  # add additional wavelength scale
    plot(('Ah1c', 'Ah1'), 1, '4h', 'True')  # add additional wavelength scale
    plot(('A1', 'A1a', 'A1b'), 1, '5', 'True')
    plot(('Ah1', 'Ah1a', 'Ah1b'), 1, '5h', 'True')
    # plot(('A1', 'A1c', 'A1d'), 1, '6', 'True')
    # plot(('Ah1', 'Ah1c', 'Ah1d'), 1, '6h', 'True')
    # illustration 7 requires a different kind of plot
    plot(('D1',), 1, '8', 'True')  # and D2 and D3
    plot(('A4', 'A5'), 1, '9-1', 'True') # and A6
    plot(('A4', 'A5'), 2, '9-2', 'True') # and A6
    
def plot(id=('A2', 'A1'), compno=1, name='0', show_func=False,
         varlist=varlist_std, ylim=ylim_std):
    print name, ':', id
    res = [fit_results(i) for i in id]
    if show_func:
        func = [fit_func(i) for i in id]
    else:
        func = None
    nvar = len(varlist)
    fig = pyplot.figure(figsize=(nvar, nvar*3))
    fig.subplots_adjust(bottom=0.05, top=0.94, left=0.2, right=0.95, hspace=0.05)
    for i, v in enumerate(varlist):
        ax = make_bands_plot(fig, (nvar, 1, i+1), labels[v], i==0, i==nvar-1)
        if v == 'MAG':
            pyplot.plot(w, simmag, '-xk')
        plotres(res, id, 'COMP%i_%s'%(compno, v), func)
        pyplot.ylim(ylim[v])
        if i==0:
            pyplot.legend(loc='lower right', numpoints=1, prop={'size': 16}) 
    fig.savefig('illustration_%s.pdf'%name)
    pyplot.close('all')

def plotres(res, id, field, func=None):
    nid = len(id)
    mec = ['black', None] * (1+nid//2)
    mfc = ['white', 'black'] * (1+nid//2)
    #color = ['grey', 'grey'] * (1+nid//2)
    color = ['Green', 'MediumPurple', 'Orange', 'MediumTurqoise']
    ymin, ymax = (1e99, -1e99)
    for i, iid in enumerate(id):
        if nid%2 == 0:
            x = w + 100 * (1+i//2) * (-1)**i
        else:
            x = w + 100 * (i//2) * (-1)**i            
        if func is not None and func[i] is not None:
            plotfunc(func[i][field], wlfunc=wlfuncs.get(iid), color=color[i])
        pyplot.errorbar(x, res[i][field], res[i][field+'_ERR'], color=color[i],
                        marker=marker[i//2], mec=mec[i], markerfacecolor=mfc[i], linestyle='',
                        label=iid)
        ymin = min(ymin, (res[i][field]-res[i][field+'_ERR']).min())
        ymax = max(ymax, (res[i][field]+res[i][field+'_ERR']).max())
    yrange = ymax - ymin
    ymin -= 0.05 * yrange
    ymax += 0.05 * yrange
    pyplot.ylim(ymin, ymax)

def plotfunc(func, wlfunc=None, color='red', label=''):
    dx = (xlim[1] - xlim[0]) / 1000.0
    x = numpy.arange(xlim[0], xlim[1]+dx/2.0, dx)
    if wlfunc is None:
        xfunc = x
    else:
        xfunc = wlfunc(x)
    y = func(xfunc)
    return pyplot.plot(x, y, '-', color=color, label=label)

def fit_results(f):
    fn = 'fits/%s/fit%s.fits'%(f,f)
    r = None
    if os.path.exists(fn):
        r = pyfits.getdata('fits/%s/fit%s.fits'%(f,f), 'final_band')
    else:
        r = numpy.concatenate([pyfits.getdata('fits/%s/fit%s%s.fits'%(f,f, b), 'final_band')
                              for b in bands])
    return r

def fit_func(f):
    fn = 'fits/%s/fit%s.fits'%(f,f)
    if os.path.exists(fn):
        r = {}
        d = pyfits.getdata('fits/%s/fit%s.fits'%(f,f), 'band_info')
        low = d.field('wl').min()
        high = d.field('wl').max()
        d = pyfits.getdata('fits/%s/fit%s.fits'%(f,f), 'final_cheb')
        for n in d.names:
            r[n] = Chebyshev(d.field(n), (low, high))
    else:
        r = None
    return r
    

def make_bands_plot(fig, subplot=111, ylabel='', top=True, bottom=True):    
    ax1 = fig.add_subplot(*subplot)
    ax2 = ax1.twiny()
    ax1.set_ylabel(ylabel)
    ax1.set_xlim(xlim)
    if top:
        ax2.set_xlabel('wavelength, $\AA$')
    else:
        ax2.set_xticklabels([])
    ax2.set_xlim(xlim)
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
