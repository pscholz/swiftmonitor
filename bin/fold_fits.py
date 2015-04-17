#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import swiftmonitor.utils as smu
from swiftmonitor import model_profile
from optparse import OptionParser

parser = OptionParser("Usage: %prog [options] fitsfile",version="%prog 1.0")
parser.add_option("-p", "--par",
		  dest="parfile", type='string',
		  help="Name of par file with pulsar ephemeris.",
		  default=None)
parser.add_option("-n", "--nbins",
		  dest="nbins", type='int',
		  help="Number of bins to fold. Default is 16.",
		  default=16)
parser.add_option("-c", "--cycles",
		  dest="ncycles", type='int',
		  help="Number of cycles of pulsar to plot. Default is 2.",
		  default=2)
parser.add_option("-H", "--H-test",
      action="store_true", dest="H_test",
      help="If True, compute H-test statistic, plot best profile.")
 
(options,args) = parser.parse_args()


 
bins, folded = smu.fold_fits(args[0],options.parfile,nbins=options.nbins)
if options.H_test:
    H, Nharm, Hfpp=smu.h_test_obs(args[0],options.parfile)
    prof_mod = model_profile.makeProfileModel("fourier", np.array([np.arange(0, len(folded)),folded]).T, nharm=Nharm)

ax=plt.subplot(111)
plot_bins = np.array([])
for i in range(options.ncycles):
    plot_bins = np.append(plot_bins,i+bins)

plot_fold = np.tile(folded,options.ncycles) 

plt.step(plot_bins,plot_fold,where='mid',c='k')
plt.errorbar(plot_bins,plot_fold,np.sqrt(plot_fold),fmt='ko')
if options.H_test:
    x = np.linspace(0, options.ncycles, 5*options.ncycles*options.nbins)
    plt.plot(x,prof_mod.prof_mod(x)*np.mean(folded), lw=2, color='r') 
    plt.text(0.7, 0.9, '$P_{fa}$='+str(Hfpp), transform = ax.transAxes )
plt.xlabel('Phase')
plt.ylabel('Counts')

plt.show()
