#!/usr/bin/env python

import pyfits
import numpy as np
import sys
import os.path

def writeinf(fits,fnroot,nbins):
  header=fits[1].header
  RAdeg=header['RA_NOM']
  RAh=RAdeg/360.*24
  RAmin=(RAh-int(RAh))*60
  RAsec=(RAmin-int(RAmin))*60
  RA=str(int(RAh))+':'+str(int(RAmin))+':'+str(RAsec)
  
  DECdeg=header['DEC_NOM']
  DECmin=abs((DECdeg-int(DECdeg))*60.)
  DECsec=(DECmin-int(DECmin))*60
  DEC=str(int(DECdeg))+':'  +str(int(DECmin))+':'+str(DECsec)

  stop=header['TSTOP']/86400.0
  start=header['TSTART']/86400.0

  #Write the .inf file# assuming Swift XRT in wt mode
  g=open(fnroot + '.inf','w+')
  g.write(' Data file name without suffix          =  '+ fnroot 
	  + '\n Telescope used                         =  Swift'
	  + '\n Instrument used                        =  XRT'
	  + '\n Object being observed                  =  ' + header['OBJECT'] 
	  + '\n J2000 Right Ascension (hh:mm:ss.ssss)  =  ' + RA
	  + '\n J2000 Declination     (dd:mm:ss.ssss)  =  ' + DEC
	  + '\n Data observed by                       =  swift'
	  + '\n Epoch of observation (MJD)             =  ' + str(header['MJDREFF'] + header['MJDREFI']+ (start))
	  + '\n Barycentered?           (1=yes, 0=no)  =  1'
	  + '\n Number of bins in the time series      =  ' + str(nbins) 
	  + '\n Width of each time series bin (sec)    =  ' + str(binsize)
	  + '\n Any breaks in the data? (1=yes, 0=no)  =  0'
	  + '\n Type of observation (EM band)          =  X-ray'
	  + '\n Field-of-view diameter (arcsec)        =  1440.0'
	  + '\n Central energy (kev)                   =  ---'
	  + '\n Energy bandpass (kev)                  =  ---'
	  + '\n Data analyzed by                       =  me'
          + '\n Any additional notes:'
	  )

def events2dat(fname,binsize):

  fnroot = os.path.splitext(fname)[0]
  datfname = "%s.dat" % fnroot

  fits = pyfits.open(fname)

  times = fits[1].data['Time']

  bins = np.arange(times[0],times[-1] + binsize, binsize)

  binned = np.histogram(times,bins=bins)[0]
  nbins = len(binned)

  # make sure its divisible by 2 for FFT
  if nbins % 2 != 0:
     binned = np.append(binned,0.0)

  np.asarray(binned, dtype='float32').tofile(datfname)

  writeinf(fits,fnroot,nbins)
  return datfname

if __name__=='__main__':
  fname = sys.argv[1]
  binsize = sys.argv[2]
  events2dat(fname,binsize)

