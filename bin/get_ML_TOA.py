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
		  help="Name of model to fit to profile. Default is fourier.",
		  default='fourier')
parser.add_option("-n", "--num-harmonics",
		  dest="num_harmonics", type='int',
		  help="Number of fourier components to use in fourier model. Default is 5.",
		  default=5)
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
parser.add_option("--scope",
		  dest="scope", type='string',
		  help="Event files are from this telescope, default swift.",
		  default='swift')
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
		  default=False)
parser.add_option("--Emin",
		  dest="emin", type='float',
		  help="Minimum energy of events to use (works only for Swift, Nustar).",
		  default=None)
parser.add_option("--Emax",
		  dest="emax", type='float',
		  help="Maximum energy of events to use (works only for Swift, Nustar).",
		  default=None)
parser.add_option("--Ntoas",
		  dest="ntoas", type='int',
		  help="Number of TOAs to extract per observation.",
		  default=None)
parser.add_option("--split_days",
		  dest="split_n_days", type='float',
		  help="If given a float, will read in multiple fits files and extract one TOA per split_n_days days. Note: must provide list of files with -l option",
		  default=None)
parser.add_option("--orbits",
		  dest="orbits", action='store_true',
		  help="Extract one TOA per orbit.",
		  default=False)
parser.add_option("--tempo2",
		  dest="tempo2", action='store_true',
		  help="Print TOA in tempo2 format.",
		  default=False)
parser.add_option("--writefile",
		  dest="writefile", type='string',
		  help="Print TOA to file (ONLY FOR tempo2 FORMAT)",
		  default=False)		  
 
(options,args) = parser.parse_args()

profile = np.loadtxt(options.profile)

prof_mod = model_profile.makeProfileModel(options.model, profile, nharm=options.num_harmonics)
if options.tempo2:
    print("FORMAT 1")
if options.periodogram:
  flist = np.loadtxt(options.periodogram,usecols=[0],dtype='S')
  epoch, frequency = np.loadtxt(options.periodogram,usecols=[1,2],dtype='float',unpack=True)
  for i,fitsfile in enumerate(flist):
    ml_toa.get_ml_toa(fitsfile, prof_mod, None, scope=options.scope, bg_counts=options.bg_counts, \
                      print_offs=options.offsets, frequency=frequency[i], epoch=epoch[i], sim=options.sim, \
                      Emin=options.emin, Emax=options.emax, gauss_err=options.gauss_err, tempo2=options.tempo2, \
                      debug=options.plot_dist, correct_pf=options.correct_pf, split_orbits=options.orbits, split_num=options.ntoas)

elif options.list and options.split_n_days:
    ml_toa.get_ml_toa(options.list, prof_mod, options.parfile, scope=options.scope, \
                      print_offs=options.offsets, sim=options.sim, bg_counts=options.bg_counts, \
                      Emin=options.emin, Emax=options.emax, gauss_err=options.gauss_err,
                      tempo2=options.tempo2, debug=options.plot_dist,
                      correct_pf=options.correct_pf,  split_n_days=options.split_n_days,
                      split_orbits=options.orbits, split_num=options.ntoas,
                      writefile=options.writefile)

else:
  if options.list:
    flist = np.loadtxt(options.list,dtype='S') 
  else:
    flist = args[0:]

  for fitsfile in flist:
    ml_toa.get_ml_toa(fitsfile, prof_mod, options.parfile, scope=options.scope, \
                      print_offs=options.offsets, sim=options.sim, bg_counts=options.bg_counts, \
                      Emin=options.emin, Emax=options.emax, gauss_err=options.gauss_err,
                      tempo2=options.tempo2, debug=options.plot_dist,
                      correct_pf=options.correct_pf,  split_n_days=options.split_n_days,
                      split_orbits=options.orbits, split_num=options.ntoas,
                      writefile=options.writefile)

