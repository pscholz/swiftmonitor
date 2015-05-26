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
parser.add_option("-H", "--H-test",
    action="store_true", dest="H_test",
    help="If True, compute H-test statistic, plot best profile.")
parser.add_option("-S", "--save-plot",
    dest="plot_fn", type='string',
    help="Save plot to filename given.")      
parser.add_option("-o", "--output",
    dest="prof_fn", type='string',
    help="Save profile in ascii format to filename given.")      
parser.add_option("-l", "--list",
    dest="list", type='string',
    help="File with list of event files to get events from.",
    default=None)      
 
(options,args) = parser.parse_args()

if options.list:
    EVTs = np.loadtxt(options.list, dtype='S')#.T[0]
    phases = np.zeros(0)
    for evt in EVTs:                     
        phases = np.append(phases,smu.fits2phase(evt,options.parfile, 
                           scope=options.scope, Emin=options.emin, 
                           Emax=options.emax))

else:                      
    phases = smu.fits2phase(args[0],options.parfile, scope=options.scope,
                            Emin=options.emin, Emax=options.emax)
bins, folded = smu.fold_phases(phases,nbins=options.nbins)

if options.H_test:
    H, Nharm, Hfpp=smu.h_test(phases)
    prof_mod = model_profile.makeProfileModel("fourier",
               np.array([np.arange(0, len(folded)),folded]).T, nharm=Nharm)

ax=plt.subplot(111)
     
plot_bins = np.array([])
for i in range(options.ncycles):
    plot_bins = np.append(plot_bins,i+bins)

plot_fold = np.tile(folded,options.ncycles) 

plt.step(plot_bins,plot_fold,where='post',c='k')
plt.errorbar(plot_bins+0.5/options.nbins,plot_fold,np.sqrt(plot_fold),fmt='ko')

if options.H_test:
    H, Nharm, Hfpp=smu.h_test(phases)
    prof_mod = model_profile.makeProfileModel("fourier", \
               np.array([np.arange(0, len(folded)),folded]).T, nharm=Nharm)
    x = np.linspace(0, options.ncycles, 10*options.ncycles*options.nbins)
    plt.plot(x,prof_mod.prof_mod(np.roll(x,5))*np.mean(folded), lw=2, color='r') 
    plt.text(0.99, 0.99, "$N_{harm}$=%d\n$P_{fa}$=%.2E" % (Nharm, Hfpp), 
             horizontalalignment="right", verticalalignment="top",
             transform = plt.gca().transAxes)

plt.xlabel('Phase')
plt.ylabel('Counts')

if options.prof_fn:
    out_array = np.array([np.arange(len(folded)),folded,np.sqrt(folded)]).T
    np.savetxt(options.prof_fn,out_array,fmt="%-4d%6d%8.2f")

if options.plot_fn:
    plt.savefig(options.plot_fn)
else:  
    plt.show()
