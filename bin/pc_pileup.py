#!/usr/bin/env python

import pyfits
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import optimize, integrate

def king_prof(r,A):
  r_c = 5.8
  beta = 1.55
  return A * ( 1 + (r/r_c)**2 )**(-1*beta)


fits = pyfits.open(sys.argv[1])
image = fits[0].data
exposure = fits[0].header['EXPOSURE']
fits.close()

binsize = 1

ceny, cenx = np.unravel_index(np.argmax(image),np.shape(image))
y, x = np.indices(image.shape) 

r = np.hypot(x - cenx, y - ceny)

nbins = (np.round(r.max() / binsize)+1)
maxbin = nbins * binsize
bins = np.linspace(0,maxbin,nbins+1)

bin_centers = (bins[1:]+bins[:-1])/2.0
whichbin = np.digitize(r.flat,bins)

radial_prof = np.array([image.flat[whichbin==b].sum() for b in xrange(1,nbins+1)])
radial_prof_err = np.sqrt(radial_prof)

areas = ( ( bins[1:]*binsize )**2 - ( bins[:-1]*binsize )**2 ) * (2.36/60.0)**2 # in sq arcmin

prof_r = bin_centers[:15]*2.36 # in arcsec
prof_cnts = radial_prof[:15]/areas[:15]/exposure # in counts/sq arcmin/sec
prof_err = radial_prof_err[:15]/areas[:15]/exposure 

A_init = 1.0
errfunc = lambda p, x, y, err: (y - king_prof(x, p)) / err
output = optimize.leastsq(errfunc, A_init, args=(prof_r[7:], prof_cnts[7:], prof_err[7:]), full_output=1)

p1 = output[0]
cov = output[1]


plt.errorbar(prof_r,prof_cnts,prof_err,binsize*2.36/2.0,fmt="o")
rplot = np.arange(prof_r[0], prof_r[-1], 0.1)
plt.plot(rplot,king_prof(rplot,p1) )

king_area = integrate.quad(king_prof,binsize*2.36, prof_r[-1],args=p1)[0]
psf_area = np.sum(prof_cnts[1:])*binsize*2.36

print 'Pile-up Fraction:', (king_area - psf_area) / king_area

plt.xscale('log')
plt.yscale('log')
plt.xlabel('R (arcsec)')
plt.ylabel('counts / sq arcmin / sec')

plt.show()
