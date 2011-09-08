#!/usr/bin/env python
#Paul Scholz, June 2010
#Converts an ascii file of swift events into XTE MET times

import scipy.io
import sys
from optparse import OptionParser

infile = sys.argv[1]

XTEREFI = 49353.0 # integer part of XTE reference epoch in mjd
XTEREFF = 0.000696574074 # fractional part of XTE reference epoch in mjd
SWREFI = 51910 # integer part of Swift reference epoch in mjd
SWREFF = 7.4287037E-4 # fractional part of Swift reference epoch in mjd

times = scipy.io.read_array(infile)

met = (SWREFI + SWREFF - XTEREFI - XTEREFF)*86400.0 + times

for t in met:
  print '%10.12f' % t

