#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import swiftmonitor.utils as smu
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
parser.add_option("--scope",
		  dest="scope", type='string',
		  help="Event files are from this telescope, default swift.",
		  default='swift')
parser.add_option("--Emin",
		  dest="emin", type='float',
		  help="Minimum energy of events to use (works only for Swift, Nustar).",
		  default=None)
parser.add_option("--Emax",
		  dest="emax", type='float',
		  help="Maximum energy of events to use (works only for Swift, Nustar).",
		  default=None)
 
(options,args) = parser.parse_args()

bins, folded = smu.fold_fits(args[0],options.parfile,nbins=options.nbins,scope=options.scope,Emin=options.emin,Emax=options.emax)

plot_bins = np.array([])
for i in range(options.ncycles):
    plot_bins = np.append(plot_bins,i+bins)

plot_fold = np.tile(folded,options.ncycles) 

plt.step(plot_bins,plot_fold,where='mid',c='k')
plt.errorbar(plot_bins,plot_fold,np.sqrt(plot_fold),fmt='ko')

plt.xlabel('Phase')
plt.ylabel('Counts')

plt.show()
