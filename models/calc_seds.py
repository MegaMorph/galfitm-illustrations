import pyfits, numpy
from scipy.optimize import fmin_powell as fmin
from scipy.stats import scoreatpercentile
from matplotlib import pyplot

d = pyfits.getdata('/Users/spb/Work/projects/MegaMorph/gama_cats/SersicCatAllv07.fits', 1, memmap=True)

def func(r, w, order):
    f = 0.0
    for i in range(order):
        f += r[-i-1] * w**i
    return f

def minfunc(r, w, m, order):
    return ((func(r, w, order) - m)**2).sum()

def combfunc(x, m1, m2):
    x = max(0, x)
    f1 = 10**(-0.4*m1)
    f2 = 10**(-0.4*m2)
    f = x*f1 + (1-x)*f2
    m = -2.5 * numpy.log10(f)
    return m

def mincombfunc(x, m1, m2, match):
    return ((combfunc(x, m1, m2) - match)**2).sum()

def calc_seds(order=5, extremes=False):
    GI = (d.field('GAL_MAG_G') - d.field('GAL_MAG_I'))
    GI10 = scoreatpercentile(GI[(GI < 5) & (GI > -5)], 10)
    GI90 = scoreatpercentile(GI[(GI < 5) & (GI > -5)], 90)
    print GI10
    fig = pyplot.figure()
    pyplot.hist(GI, 28, range=(-2, 5))
    sp = (d.GAL_INDEX_R > 0.5) & (d.GAL_INDEX_R < 1.5)
    if extremes:
        sp &= GI < GI10
    print len(sp.nonzero()[0])
    el = (d.GAL_INDEX_R > 3.0) & (d.GAL_INDEX_R < 6.0)
    if extremes:
        el &= GI > GI90
    print len(el.nonzero()[0])
    bd = (d.GAL_INDEX_R > 1.5) & (d.GAL_INDEX_R < 3.0)
    print len(bd.nonzero()[0])
    bands = ['u', 'g', 'r', 'i', 'z', 'Y', 'J', 'H', 'K']
    w = numpy.array([3543,4770,6231,7625,9134,10305,12483,16313,22010], numpy.float)
    x = (w-10000)/10000.
    offset_sp=[]
    offset_el=[]
    offset_bd=[]
    offset=[]
    for b in bands:
        off = d.field('GAL_MAG_'+b) - d.field('GAL_MAG_R')
        ok = d.field('GAL_MAG_'+b) > 9
        ok &= d.field('GAL_MAG_'+b) < 99
        ok &= d.field('GAL_MAG_R') > 9
        ok &= d.field('GAL_MAG_R') < 99
        ok &= d.field('GAL_INDEX_R') < 6
        ok &= d.field('GAL_INDEX_R') > 0.5
        print b, len(ok.nonzero()[0])
        offset.append(numpy.median(off[ok]))
        offset_sp.append(numpy.median(off[sp&ok]))
        offset_el.append(numpy.median(off[el&ok]))
        offset_bd.append(numpy.median(off[bd&ok]))
    print offset
    print offset_sp
    print offset_el
    print offset_bd    
    offset = numpy.array(offset)
    offset_sp = numpy.array(offset_sp)
    offset_el = numpy.array(offset_el)
    offset_bd = numpy.array(offset_bd)
    #fcomb = fmin(mincombfunc, 0.5, args=[offset_sp, offset_el, offset])
    fcomb = 0.5
    print fcomb
    offset_combined = combfunc(fcomb, offset_sp, offset_el)
    print offset_combined
    # Fitting and reconstruction (disabled)
    # r = fmin(minfunc, [0.]*order, args=[x, offset, order], xtol=1e-6, ftol=1e-6, maxiter=1000, maxfun=100000)
    # r_sp = fmin(minfunc, r, args=[x, offset_sp, order], xtol=1e-6, ftol=1e-6, maxiter=1000, maxfun=100000)
    # r_el = fmin(minfunc, r, args=[x, offset_el,order], xtol=1e-6, ftol=1e-6, maxiter=1000, maxfun=100000)
    # print r
    # print r_sp
    # print r_el
    # foffset = func(r, x, order)
    # foffset_sp = func(r_sp, x, order)
    # foffset_el = func(r_el, x, order)
    # print foffset
    # print foffset_sp
    # print foffset_el
    # print max(foffset - offset)
    # print max(foffset_sp - offset_sp)
    # print max(foffset_el - offset_el)
    # pyplot.plot(w, foffset, 'ok')
    # pyplot.plot(w, foffset_sp, 'ob')
    # pyplot.plot(w, foffset_el, 'or')
    fig = pyplot.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(w, offset, 'x-k')
    ax1.plot(w, offset_sp, '*-b')
    ax1.plot(w, offset_el, '.-r')
    ax1.plot(w, offset_combined, '+-g')
    ax2 = ax1.twiny()
    ax1.set_xlim(2000, 23000)
    ax1.set_ylim(3.5, -2.5)
    ax1.set_ylabel('mag offset $$')
    ax2.set_xlabel('wavelength, $\AA$')
    ax2.set_xlim(2000, 23000)
    ax1.set_xticks(w)
    ax1.set_xticklabels(['$'+i+'$' for i in bands])
    pyplot.setp(ax1.get_xticklabels(), va='baseline')
    pyplot.setp(ax1.get_xaxis().get_major_ticks(), pad=20.)
    pyplot.setp(ax1.get_yaxis().get_major_ticks(), pad=8.)
    ax2.xaxis.labelpad = 12
    fig.subplots_adjust(top=0.85)
    if extremes:
        fig.savefig('extreme_seds.pdf')
    else:
        fig.savefig('seds.pdf')
    print
    print '### Figures for models'
    print 'Overall SED:',
    print galfit_format(offset + 15)
    print 'Overall SED (half flux):',
    print galfit_format(offset + 15 - 2.5*numpy.log10(fcomb))
    print 'Disk SED (split flux, fraction=%f):'%fcomb,
    print galfit_format(offset_sp + 15 - 2.5*numpy.log10(fcomb))
    print 'Spheroid SED (split flux, fraction=%f):'%(1-fcomb),
    print galfit_format(offset_el + 15 - 2.5*numpy.log10(fcomb))
    print 'Disk+Spheroid combined SED:',
    print galfit_format(offset_combined + 15)
    print 'Disk+Spheroid combined SED (half flux):',
    print galfit_format(offset_combined + 15 - 2.5*numpy.log10(fcomb))
    print
    print '### Figures for model C'
    for fcomb in numpy.arange(0.1, 0.95, 0.1):
        print 'Disk SED (split flux, fraction = %f):'%fcomb,
        print galfit_format(offset_sp + 15 - 2.5*numpy.log10(fcomb))
        print 'Spheroid SED (split flux, fraction = %f):'%(1-fcomb),
        print galfit_format(offset_el + 15 - 2.5*numpy.log10(1-fcomb))
    

def galfit_format(x):
    return repr([round(i, 3) for i in x]).replace(' ', '')

def mag_size():
    m = d.field('GAL_MAG_R')
    re = d.field('GAL_RE_R')
    n = d.field('GAL_INDEX_R')
    ok = (0 < m) & (m < 30) & (re > 0) & (re < 10) & (n > 0.5) & (n < 6)
    fig = pyplot.figure()
    pyplot.scatter(m[ok], re[ok], s=1, lw=0)
    fig.savefig('mag_size.pdf')

if __name__ =='__main__':
    calc_seds()
