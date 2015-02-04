#!/usr/bin/env python

import numpy as np
from swiftmonitor import ml_toa, model_profile
from optparse import OptionParser

models = model_profile.list_models()

models_str = "Available Models:  "
for model in models:
  models_str += model + '  '

parser = OptionParser("Usage: %prog [options] fitsfiles\n OR -l file_list", epilog=models_str ,version="%prog 1.0")
parser.add_option("-m", "--model",
		  dest="model", type='string',
		  help="Name of model to fit to profile.",
		  default=None)
parser.add_option("-f", "--profile",
		  dest="profile", type='string',
		  help="Name of profile to fit model to.",
		  default=None)
parser.add_option("-p", "--par",
		  dest="parfile", type='string',
		  help="Name of parfile with timing solution.",
		  default=None)
parser.add_option("-l", "--list",
		  dest="list", type='string',
		  help="File with list of event files to get TOAs from.",
		  default=None)
parser.add_option("--chandra",
		  dest="chandra", action='store_true',
		  help="Event files are from chandra, not swift.",
		  default=False)
parser.add_option("--xmm",
		  dest="xmm", action='store_true',
		  help="Event files are from xmm, not swift.",
		  default=False)
parser.add_option("--xte",
		  dest="xte", action='store_true',
		  help="Event files are from xte, not swift.",
		  default=False)
parser.add_option("--offsets",
		  dest="offsets", type='string',
		  help="File to write phase offset used for TOA.",
		  default=None)
parser.add_option("--periodogram",
		  dest="periodogram", type='string',
		  help="Name of periodogram.py output file to use as input.",
		  default=None)
parser.add_option("--sim",
		  dest="sim", action='store_true',
		  help="Determine the error in the TOA using simulations.",
		  default=False)
parser.add_option("--gauss-err",
		  dest="gauss_err", action='store_true',
		  help="Determine the error by fitting a gaussian.",
		  default=False)
parser.add_option("--plot-dist",
		  dest="plot_dist", action='store_true',
		  help="Plot the offset probability distribution for debugging.",
		  default=False)
parser.add_option("--bg-counts",
		  dest="bg_counts", type='float',
		  help="Number of counts in background.",
		  default=0)
parser.add_option("--correct-pf",
		  dest="correct_pf", action='store_true',
		  help="Rescale the template profile to have the same pulsed " \
                       + "fraction as the observation for which the TOA is being calculated.",
		  default=0)
parser.add_option("--Emin",
		  dest="emin", type='float',
		  help="Minimum energy of events to use (works only for Swift).",
		  default=None)
parser.add_option("--Emax",
		  dest="emax", type='float',
		  help="Maximum energy of events to use (works only for Swift).",
		  default=None)
parser.add_option("--tempo2",
		  dest="tempo2", action='store_true',
		  help="Print TOA in tempo2 format.",
		  default=False)
 
(options,args) = parser.parse_args()

profile = np.loadtxt(options.profile)

prof_mod = model_profile.makeProfileModel(options.model, profile)

if options.periodogram:
  flist = np.loadtxt(options.periodogram,usecols=[0],dtype='S')
  epoch, frequency = np.loadtxt(options.periodogram,usecols=[1,2],dtype='float',unpack=True)
  for i,fitsfile in enumerate(flist):
    ml_toa.get_ml_toa(fitsfile, prof_mod, None, chandra=options.chandra, xmm=options.xmm, xte=options.xte, bg_counts=options.bg_counts, \
                      print_offs=options.offsets, frequency=frequency[i], epoch=epoch[i], sim=options.sim, \
                      Emin=options.emin, Emax=options.emax, gauss_err=options.gauss_err, tempo2=options.tempo2, \
                      debug=options.plot_dist, correct_pf=options.correct_pf)

else:
  if options.list:
    flist = np.loadtxt(options.list,dtype='S') 
  else:
    flist = args[0:]

  for fitsfile in flist:
    ml_toa.get_ml_toa(fitsfile, prof_mod, options.parfile, chandra=options.chandra, xmm=options.xmm, xte=options.xte, \
                      print_offs=options.offsets, sim=options.sim, bg_counts=options.bg_counts, \
                      Emin=options.emin, Emax=options.emax, gauss_err=options.gauss_err, tempo2=options.tempo2, \
                      debug=options.plot_dist,correct_pf=options.correct_pf)

