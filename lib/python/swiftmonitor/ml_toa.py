import numpy as np
from scipy import optimize, integrate
import matplotlib.pyplot as plt
import pyfits
import sys
import psr_utils
from psr_constants import SECPERDAY
from swiftmonitor.utils import randomvariate, events_from_binned_profile
import time

sys.setrecursionlimit(100000)

calcprobtime = 0
logsumtime = 0
integratetime = 0


def logsumexp(array):
  if len(array) == 2:
    return np.logaddexp(array[0],array[1])
  return np.logaddexp(array[0],logsumexp(array[1:]))

def sw2mjd(times):
  SWREFI = 51910.0 # integer part of Swift reference epoch in mjd
  SWREFF = 7.4287037E-4 # fractional part of Swift reference epoch in mjd

  mjds = SWREFI + SWREFF + times/86400.0
  return mjds

def chandra2mjd(times):
  CHANDRAREF = 50814.0

  mjds = CHANDRAREF + times/86400.0
  return mjds

def calc_prob(phases, offset, prof_mod):
  """
  Calculates the probability of an offset given a profile model, and a list of phases.
  """
  probs = prof_mod( phases - offset ) 

  loglike = np.sum( np.log(probs) )
  return loglike


def readParfile(parname):
    """ 
    Read a tempo .par file and output a list containing:
       [F0,F1,F2,F3,F4,PEPOCH]
    """
    F = np.loadtxt(parname, dtype='S')
    freqs = np.zeros(6)
    epoch = 0
    for line in F:
      for i in range(0,5):
        if line[0].startswith('F'+str(i)):
          freqs[i] = line[1].replace('D','E')
        if line[0].startswith('PEPOCH'):
          epoch = line[1].replace('D','E')
        freqs[5] = np.double(epoch)     
    return freqs 

class PSRpar:
  '''This class contains all the relevant information for 
     extracting a TOA, which is extracted from a .par file.
     Currently, this supports up to 4 frequncy derivatives
     but NO GLITCH PARAMETERS (yet).'''
  def __init__(self,parfile):
    A = readParfile(parfile)
    self.epoch = A[5]
    self.f0 = A[0]
    self.fdots = A[1:5]
  def __repr__(self):
    return repr((self.f0, self.fdots, self.epoch))    

def getErrMid50(offsets, Prob, del_off):
    '''Centers the Probability distribution at where 50% of the 
       probability is at 0.5. IT then integrates the Probability on either side of
       the max, returns average distance to 34.1% of area covered'''
    center=np.argmax(Prob)
    CProb = np.roll(Prob, len(offsets)/2 - np.argmax(Prob))
    Area=integrate.cumtrapz(CProb,dx=del_off)
    center2=(np.argwhere(Area>0.5)[0])[0]-len(offsets)/2
    newCProb=np.roll(CProb,-center2)
    halfN1=newCProb[:len(offsets)/2]
    halfN2=newCProb[len(offsets)/2:]
    AreaN1=integrate.cumtrapz(halfN1,dx=del_off)
    AreaN2=integrate.cumtrapz(halfN2,dx=del_off)
    i=0
    init=10000
    a=init
    j=0
    b=init       
    while a==init and i<len(AreaN2):
      if AreaN2[i]>0.341:
        a=i
      i=i+1
    while b==init and j<len(AreaN1):
      if AreaN1[j]>0.341:
       b=len(AreaN1)-j
      j=j+1  
    offset=(center2+center)*del_off
    sigma=(a+b)/2.0*del_off
    #plt.plot(offsets, CProb) 
    #plt.vlines(offsets[len(CProb)/2+center2], 0, max(CProb)*1.05,linestyles='dotted' )
    #plt.errorbar(offsets[len(CProb)/2+center2],0.5*max(CProb), xerr=sigma,fmt='o')
    #plt.show()
    return offset, sigma 

def getErr(offsets, Prob, del_off):
    '''Integrates the Probability on either side of
       the max, returns average distance to 34.1% of area covered'''

    maxoff = offsets[np.argmax(Prob)]
    CProb = np.roll(Prob, len(offsets)/2 - np.argmax(Prob))
    half1=CProb[:len(offsets)/2]
    half2=CProb[len(offsets)/2:]
    Area1=integrate.cumtrapz(half1,dx=del_off)
    Area2=integrate.cumtrapz(half2,dx=del_off)
    i=0
    a=0
    j=0
    b=0       
    while a==0:
      if Area2[i]>0.341:
        a=i
      i=i+1
    while b==0:
      if Area1[j]>0.341:
        b=len(Area1)-j
      j=j+1  
    sigma=(a+b)/2.0*del_off
    #plt.errorbar(offsets[np.argmax(CProb)],0.5*max(CProb), xerr=sigma,fmt='o')
    #plt.plot(offsets, CProb)
    #plt.axvline(0.5,ls='dotted')
    #plt.show()
    return maxoff, sigma 

def simErr(prof_mod,N_counts,phases,from_template=True):
  sim_offsets = []
  N_sim = 1000
  N_bins = 32
  run = 0
  if not from_template:
    bins = np.linspace(0,1,N_bins+1)
    folded = np.histogram(phases,bins)[0] 

  while run < N_sim:
    sys.stderr.write("Sim %d %% Complete \r" % (run*100.0/N_sim))
    sys.stderr.flush()
    if from_template:
      simmed_phases = randomvariate(prof_mod,N_counts) 
    else:
      folded = np.random.poisson(folded)
      simmed_phases = events_from_binned_profile(folded)
    sim_offset, sim_error = calc_toa_offset(simmed_phases,prof_mod,no_err=True)
    sim_offsets.append((sim_offset+0.5) % 1.0)
    run += 1

  #median = np.median(sim_offsets)
  #mad = np.median(np.abs(sim_offsets-median))
  #plt.hist(sim_offsets)
  #print "Std dev offsets",np.std(sim_offsets)
  #print "MAD offsets",mad
  #plt.show()
  sys.stderr.write("Sim 100%% Complete \n")
  return np.std(sim_offsets)


def calc_toa_offset(phases,prof_mod,sim_err=False,no_err=False, bg_counts=0):
  global calcprobtime
  global logsumtime 
  global integratetime

  probs = []
  del_off = 0.0001
  offsets = np.arange(0,1,del_off)
  offsets = np.append(offsets,1.0)

  starttime = time.time()
  for offset in offsets:
    prob =  calc_prob(phases, offset, prof_mod)  
    probs.append(prob)
  calcprobtime += time.time() - starttime

  starttime = time.time()
  probs = probs - logsumexp(probs)
  logsumtime += time.time() - starttime

  starttime = time.time()
  probs = np.exp(np.array(probs))
  area = integrate.trapz(probs,dx=del_off)
  probs_norm = probs/area
  integratetime += time.time() - starttime

  if sim_err:
    maxoff = offsets[np.argmax(probs)]
    error = simErr(prof_mod,len(phases)-bg_counts,phases,from_template=True) 
  elif no_err:
    maxoff = offsets[np.argmax(probs)]
    error = None
  else:
    maxoff, error = getErrMid50(offsets, probs_norm, del_off)
  
  return maxoff, error


def get_ml_toa(fits_fn, prof_mod, parfile, chandra=False, xmm=False, print_offs=False, frequency=None, epoch=None, sim=False, bg_counts=0):

  fits = pyfits.open(fits_fn)
  swift_t = fits[1].data['Time']
  if not chandra:
    exposure = fits[0].header['EXPOSURE']
  obsid = fits[0].header['OBS_ID']

  if chandra and xmm:
    raise ValueError('Data can only be from one of Chandra and XMM!')
  elif chandra or xmm:
    t = chandra2mjd(swift_t) # XMM and chandra use same MJDREF
  else:
    t = sw2mjd(swift_t)

  if frequency and epoch:
    par = lambda: None
    par.epoch = epoch
    par.f0 = frequency
    par.fdots = [0.0,0.0,0.0,0.0]
  else:
    par = PSRpar(parfile)

  sys.stderr.write('Measuring TOA for %s\n' % obsid)

  phases = psr_utils.calc_phs(t, par.epoch, par.f0, par.fdots[0], par.fdots[1], 
                                 par.fdots[2], par.fdots[3]) 

  maxoff, error = calc_toa_offset(phases,prof_mod,sim_err=sim,bg_counts=bg_counts)

  if chandra or xmm:
    midtime = ( chandra2mjd(fits[0].header['TSTART']) + chandra2mjd(fits[0].header['TSTOP']) ) / 2.0
  else:
    midtime = ( sw2mjd(fits[0].header['TSTART']) + sw2mjd(fits[0].header['TSTOP']) ) / 2.0
  p_mid = 1.0/psr_utils.calc_freq(midtime, par.epoch, par.f0, par.fdots[0])

  t0 = psr_utils.calc_t0(midtime, par.epoch, par.f0, par.fdots[0])
  t0i = int(t0)
  t0f = t0 - t0i

  toaf = t0f + maxoff*p_mid / SECPERDAY
  newdays = int(np.floor(toaf))
  psr_utils.write_princeton_toa(t0i+newdays, toaf-newdays, error*p_mid*1.0e6, 0000, 0.0, name=obsid) 

  if print_offs:
    print "\t",error*p_mid*1.0e6,"\t",exposure
    print obsid,"\tOffset:",maxoff,"+/-",error 

  fits.close()

  global calcprobtime
  global logsumtime 
  global integratetime

  sys.stderr.write('\tCalc Prob: %f s\n' % calcprobtime)
  sys.stderr.write('\tLog Sum: %f s\n' % logsumtime)
  sys.stderr.write('\tIntegrate Norm: %f s\n' % integratetime)
