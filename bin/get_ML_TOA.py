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
parser.add_option("--offsets",
		  dest="offsets", action='store_true',
		  help="Print phase offset used for TOA.",
		  default=False)

 
(options,args) = parser.parse_args()

profile = np.loadtxt(options.profile)

prof_mod = model_profile.makeProfileModel(options.model, profile)

if options.list:
  flist = np.loadtxt(options.list,dtype='S') 
else:
  flist = args[0:]

for fitsfile in flist:
  ml_toa.get_ml_toa(fitsfile, prof_mod, options.parfile, chandra=options.chandra, print_offs=options.offsets)
