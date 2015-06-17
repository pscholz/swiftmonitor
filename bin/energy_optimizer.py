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
    default=0.0)
parser.add_option("--Emax",
    dest="emax", type='float',
    help="Maximum energy of events to use (works only for Swift, Nustar).",
     default=None)		  	
parser.add_option("-s", "--steps",
    dest="steps", type='int',
    help="Number of steps in Enengy optimization.",
    default=16)     	
parser.add_option( "--2-D",
    dest="twoD", action='store_true',
    help="Optimize in 2-D, default, only cut off low. ")      
      
#parser.add_option("-H", "--H-test",
#    action="store_true", dest="H_test",
#    help="If True, compute H-test statistic, plot best profile.")
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
if options.emin ==None or options.emax == None:
    print('Must give Enegy bounds')
    raise SystemExit
if options.list:
    EVTs = np.loadtxt(options.list, dtype='S')#.T[0]
    phases = np.zeros(0)
    Es = np.zeros(0)
    for evt in EVTs:                     
        phase, E = smu.fits2phase(evt,options.parfile, 
                           scope=options.scope, Emin=options.emin, 
                           Emax=options.emax, give_t_E=True)
        phases=np.append(phases, phase)
        Es=np.append(Es, E)                   

else:                      
    phases, Es = smu.fits2phase(args[0],options.parfile, 
                           scope=options.scope, Emin=options.emin, 
                           Emax=options.emax, give_t_E=True) 

ECuts=np.linspace(options.emin, options.emax, options.steps)
Cuts = smu.energy2chan(ECuts, scope=options.scope)
if options.twoD:
    hs=np.zeros((options.steps, options.steps))
    for i in range(0, options.steps):
        for j in range(i, options.steps):
           filt=np.where((Es>Cuts[i])*(Es<Cuts[j])) 
           hs[i][j]=smu.h_test(phases[filt])[0]
    hsm=np.ma.masked_where(np.isnan(hs),hs)  
    plt.pcolor(ECuts, ECuts, hsm, cmap='Greens')
    maxhs=np.nanmax(hs)
    bestEmin=ECuts[np.nanargmax(hs)/options.steps]
    bestEmax=ECuts[np.nanargmax(hs)%options.steps]
    plt.text(0.01, 0.99, "Best $E_{min}$ = %.2f keV\nBest $E_{max}$ = %.2f keV\nBest h-score = %.2E" % (bestEmin,bestEmax, maxhs), 
             horizontalalignment="left", verticalalignment="top",
             transform = plt.gca().transAxes)
    plt.colorbar()
    plt.xlabel('$E_{max}$ (keV)')
    plt.ylabel('$E_{min}$ (keV)')
         
else:
    hs = np.zeros((0))
    for cut in Cuts:
        filt=np.where(Es>cut)
        hs=np.append(hs, smu.h_test(phases[filt])[0])    
    plt.plot(ECuts, hs, 'ko')
    plt.xlabel('E (keV)')
    plt.ylabel('h-score')
    maxhs=np.nanmax(hs)
    bestEmin=ECuts[np.nanargmax(hs)]
    plt.text(0.99, 0.99, "Best Emin = %.2f keV\nBest h-score = %.2E" % (bestEmin, maxhs), 
             horizontalalignment="right", verticalalignment="top",
             transform = plt.gca().transAxes)
   

if options.plot_fn:
    plt.savefig(options.plot_fn)
else:  
    plt.show()    

