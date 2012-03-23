#!/usr/bin/env python
import swiftmonitor.observation as swift_obs
import numpy as np
import os
import sys

from optparse import OptionParser

parser = OptionParser("Usage: %prog pickle_files\n\t OR %prog -l pickle_list", version="%prog 1.0")
parser.add_option("-l", "--list",
                  dest="list", type='string',
                  help="List of pickle files to re-extract spectra.",
                  default=None)
parser.add_option("-m", "--mode",
                  dest="mode", type='string',
                  help="Mode of observation: pc or wt.",
                  default=None)
parser.add_option("--grade0",
                  dest="grade0", action='store_true',
                  help="Extract a grade 0 spectrum.",
                  default=False)
parser.add_option("--backscal",
                  dest="backscal", type='float',
                  help="Value for BACKSCAL keyword.",
                  default=1.0)


(options,args) = parser.parse_args()

if options.list:
  pickle_list = np.loadtxt(options.list, dtype='S')
else:
  pickle_list = args

if options.grade0:
  grade='0'
else:
  grade=None

for pickle_file in pickle_list:

  ob = swift_obs.load(pickle_file)

  ob.extract_region()
  ob.extract_spectrum(grade=grade)
  if options.mode == 'wt':
    ob.correct_backscal(value=options.backscal)
  ob.save()
