#!/usr/bin/env python

from xspec import *
import sys
import numpy as np
import pyfits


obslist = sys.argv[1]
outfile = sys.argv[2]

specs = np.loadtxt(obslist,dtype='S200')


SWREFI = 51910 # integer part of Swift reference epoch in mjd
SWREFF = 7.4287037E-4 # fractional part of Swift reference epoch in mjd


tstarts = np.array([])
tstops = np.array([])

i = 1
loadstr = ''
for spec in specs:
  loadstr += str(i) + ':' + str(i) + " " + spec + " "
  i += 1

  fits = pyfits.open(spec)
  tstarts = np.append(tstarts,float(fits[0].header['TSTART']))
  tstops = np.append(tstops,float(fits[0].header['TSTOP']))

mjd_starts = SWREFI + SWREFF + tstarts/86400.0
mjd_stops = SWREFI + SWREFF + tstops/86400.0

AllData(loadstr)

AllData.ignore("bad")
AllData.ignore("**:**-0.5 10.0-**")

model = Model("ph(po+bb)")

errorstr = '1.0 1 '
# untie non nH pars and set initial kT = 1
for i in range(5*len(specs)):
  imod = i % 5
  if imod:
    AllModels(i/5+1)(i+1).untie() 
  if imod == 3:
    AllModels(i/5+1)(i+1).values = 1
  if not imod == 0 and not imod == 2 and not imod == 4:
    errorstr += '%d ' % (i+1) 

AllModels.show()

Fit.query = 'yes'
Fit.perform()
Fit.perform()

Plot.xAxis = "keV"
Plot.device = "/ps"
Plot("data","resid")
Plot()

Fit.error(errorstr)
AllModels.calcFlux("1. 10. err")

outf = open(outfile,'w')
nH = AllModels(1)(1)
outf.write('# nH: %f (%f - %f)\tChisq: %f\tDOF: %f\n' % ( nH.values[0], nH.error[0], nH.error[1], Fit.statistic, Fit.dof) )
outf.write('# obsid\t\t\tmjdstart\tmjdstop\tflux\tflux_low\tflux_high\tgamma\tgamma_low\tgamma_high\tkT\tkT_low\tkT_high\n')

for i in range(len(specs)):
  flux = AllData(i+1).flux
  output = '%s\t%f\t%f\t%E\t%E\t%E\t%f\t%f\t%f\t%f\t%f\t%f\n' %\
           ( specs[i], mjd_starts[i], mjd_stops[i], flux[0], flux[1], flux[2], AllModels(i+1)(5*i+2).values[0], 
             AllModels(i+1)(5*i+2).error[0], AllModels(i+1)(5*i+2).error[1], AllModels(i+1)(5*i+4).values[0], 
             AllModels(i+1)(5*i+4).error[0], AllModels(i+1)(5*i+4).error[1] )             

  outf.write(output)  
  

outf.close()
