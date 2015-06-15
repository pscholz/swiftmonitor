import numpy as np
import astropy.io.fits as pyfits
import sys 


SECPERDAY=86400.0

string_pars = ['PSR','PSRJ','EPHEM','CLK','BINARY','JUMP',\
               'UNITS','TZRSITE','TIMEEPH','T2CMETHOD',\
               'CORRECT_TROPOSPHERE','PLANET_SHAPIRO','DILATEFREQ',\
               'RAJ', 'DECJ']
  
def calc_freq(MJD, refMJD, *args):
    """
    calc_freq(MJD, refMJD, *args):
        Return the instantaneous frequency at an MJD (can be an array)
            given a reference MJD and the rotational freq (f0) and
            optional freq derivs (f1...) as ordered in the *args
            list (e.g. [f0, f1, f2, ...]).
            Note: from psr_utils (PRESTO)
    """
    t = (MJD-refMJD)*SECPERDAY
    n = len(args) # polynomial order
    taylor_coeffs = np.concatenate(([1.0],
                                     np.cumprod(1.0/(np.arange(float(n-1))+1.0))))
    p = np.poly1d((taylor_coeffs * args)[::-1])
    return p(t)

               
def calc_t0(MJD, refMJD, *args):
    """
    calc_t0(MJD, refMJD, *args):
        Return the closest previous MJD corresponding to phase=0 of the pulse.
        *args are the spin freq (f0) and optional freq derivs (f1...)
         Note: from psr_utils (PRESTO)
    """
    phs = calc_phs(MJD, refMJD, *args)
    p = 1.0 / calc_freq(MJD, refMJD, *args)
    return MJD - phs*p/SECPERDAY

def write_tempo2_toa(toa_MJDi, toa_MJDf, toaerr, freq, dm, obs='@', name='unk', flags=""):
    """
    Write Tempo2 format TOAs.
    Note that first line of file should be "FORMAT 1"
    TOA format is "file freq sat satErr siteID <flags>"
    Note: from psr_utils (PRESTO)
    """
    toa = "%5d"%int(toa_MJDi) + ("%.13f"%toa_MJDf)[1:]
    if dm != 0.0:
        flags += "-dm %.4f" % (dm,)
    print "%s %f %s %.2f %s %s" % (name,freq,toa,toaerr,obs,flags)

def write_princeton_toa(toa_MJDi, toa_MJDf, toaerr, freq, dm, obs='@', name=' '*13):
    """
    Princeton Format

    columns     item
    1-1     Observatory (one-character code) '@' is barycenter
    2-2     must be blank
    16-24   Observing frequency (MHz)
    25-44   TOA (decimal point must be in column 30 or column 31)
    45-53   TOA uncertainty (microseconds)
    69-78   DM correction (pc cm^-3)
    Note: from psr_utils (PRESTO)
    """
    # Splice together the fractional and integer MJDs
    toa = "%5d"%int(toa_MJDi) + ("%.13f"%toa_MJDf)[1:]
    if dm!=0.0:
        print obs+" %13s %8.3f %s %8.2f              %9.4f" % \
              (name, freq, toa, toaerr, dm)
    else:
        print obs+" %13s %8.3f %s %8.2f" % \
              (name, freq, toa, toaerr)

def calc_phs(MJD, refMJD, *args):
    """
    calc_phs(MJD, refMJD, *args):
        Return the rotational phase (0-1) at MJD (can be an array)
            given a reference MJD and the rotational freq (f0) and
            optional freq derivs (f1...) as ordered in the *args
            list (e.g. [f0, f1, f2, ...]).
            Note: from psr_utils (PRESTO)
    """
    t = (MJD-refMJD)*SECPERDAY
    n = len(args) # polynomial order
    nargs = np.concatenate(([0.0], args))
    taylor_coeffs = np.concatenate(([0.0],
                                     np.cumprod(1.0/(np.arange(float(n))+1.0))))
    p = np.poly1d((taylor_coeffs * nargs)[::-1])
    return p(t)% 1.0
               

class Par:
    def __init__(self,name,value,error=None,fit=None):
        self.name = name
        self.value = value
        self.error = error
        self.fit = fit
    def __repr__(self):
        return "Name: " + self.name + " Value: " + str(self.value) + \
               " Error: " + str(self.error) + " Fit: " + str(self.fit)

def read_parfile(parfn):
    pars = {}
    parf = open(parfn,'r')
    for line in parf.readlines():
        split = line.split()
        if split[0] == 'C': #its a comment
            continue
        if len(split) == 1: #has no value
            continue
            
        value = split[1] if split[0] in string_pars else float(split[1])
        
        if len(split) == 2:
            pars[split[0]] = Par(split[0],value)
        elif len(split) == 3:
            if (split[2] == '0') or (split[2] == '1'):
                pars[split[0]] = Par(split[0],value,fit=int(split[2]))
            else:
                pars[split[0]] = Par(split[0],value,error=float(split[2]))
        elif len(split) == 4:
            pars[split[0]] = Par(split[0],value,error=float(split[3]), \
                                 fit=split[2])
    return pars

def energy2chan(E, scope='swift'):
    """Takes a given Energy or array of Energies in keV, and converts them
       into the channel of either the 'PI' or 'PHA' fits column.
       INPUTS: 
              E - energy or energies to convert
              scope - which telescope to use (default 'swift')
       OUTPUTS:
              chans - 'PI' or 'PHA' fits column values
    """
    if scope == 'swift':
        chans = E * 100.0
    elif scope == 'nustar':
        chans = (E - 1.6) / 0.04
    else:
        sys.stderr.write('Warning: scope not found, assuming channels!\n')
        chans = E
    return chans
  

def fits2times(evtname,scope='swift',Emin=None, Emax=None, give_t_E=False):
    """Given a FITS file, this will read the reference epochs,
       and convert MET into MJD
       INPUTS:
              evtname - name of FITS file to read
       OUTPUTS:
             t - Event arrival times in Modified Julian Dates    
             if  give_t_E: - E - Energy of phases
       
    """
    fits = pyfits.open(evtname)
    t = fits[1].data['time']
    t = t 
    t = t / 86400.0

    try:
        t = t + fits[1].header['MJDREFI'] + fits[1].header['MJDREFF']
 
    except (KeyError):
        t = t + fits[1].header['MJDREF'] 

    if "PI" in fits[1].columns.names:
      Echans = fits[1].data['PI']
    elif "PHA" in fits[1].columns.names:
      Echans = fits[1].data['PHA']
    else:
        sys.stderr.write('No Energy Column\n')
        Emin, Emax = None, None

    if (Emin and Emax):
        PI_min = energy2chan(Emin, scope)
        PI_max = energy2chan(Emax, scope)
        t = t[(Echans < PI_max) & (Echans > PI_min)]
    elif Emin:
        PI_min = energy2chan(Emin, scope)
        t = t[(Echans > PI_min)]
    elif Emax:
        PI_max = energy2chan(Emax, scope)
        t = t[(Echans < PI_max)]
    else:
        sys.stderr.write('No Energy Filter\n')

    fits.close()
    if give_t_E:
        return t, Echans
    else:    
        return t  

def fits2phase(fits_fn, par_fn, scope='swift',Emin=None, Emax=None, give_t_E=False):
    """Given a FITS file and a parfile, this will read the reference epochs,
       and convert into phases
       INPUTS:
           fits_fn - name of FITS file to read
       OUTPUTS:
           phase - Pulsar phase, from 0-1.     
       
    """
    if  give_t_E:
        t, E = fits2times(fits_fn,Emin=Emin,Emax=Emax,scope=scope, give_t_E=True)
    else:    
        t = fits2times(fits_fn,Emin=Emin,Emax=Emax,scope=scope, give_t_E=False)    
    phases = times2phases(t, par_fn)
    if  give_t_E:
        return phases, E
    else:   
        return phases


def fold_fits(fits_fn, par_fn, nbins=32, scope='swift', Emin=None, Emax=None):
    """Given a FITS file and a parfile, this will convert the events to phases
       and fold them, returning a histogram of folded phases.
       INPUTS:
           fits_fn - name of FITS file to read
       OUTPUTS:
           bins - the left bin edges for each bin
           folded - the number of events in each bin  
       
    """
    phases = fits2phase(fits_fn, par_fn, scope=scope, Emin=Emin, Emax=Emax)

    return fold_phases(phases, nbins=nbins)

def fold_phases(phases, nbins=32):
    """Given list of phase (e.g. from fits2phases), this will bin them 
       into a histogram with nbins between 0 and 1.
       INPUTS:
           phases - a list or array of phases
       OUTPUTS:
           bins - the left bin edges for each bin
           folded - the number of events in each bin  
    """

    bins = np.linspace(0,1,nbins+1) # add an extra bin for np.histogram's 
                                    # rightmost limit
    folded = np.histogram(phases,bins)[0]

    return bins[:-1],folded

def times2phases(t, par_fn):
    """Given an array of times and a parfile, this will read the reference epoch
       and frequency parameters, and convert into phases
       INPUTS:
           t -an array of photon arrival times in MJD
           par_fn - 
       OUTPUTS:
           phase - Pulsar phase, from 0-1.     
       
    """
    par = read_parfile(par_fn)

    phs_args = [ t, par['PEPOCH'].value, par['F0'].value ]

    for i in range(12):
        fdot_name = 'F' + str(i+1)
        if fdot_name in par.keys():  
            phs_args.append(par[fdot_name].value)
        else:
            phs_args.append(0.0)

    phases = calc_phs(*phs_args) 
    return phases
    
def events_from_binned_profile(profile): 
    binsize = 1.0 / len(profile)
    phases = np.array([])
    for i,counts in enumerate(profile):
      phases = np.append(phases, np.random.rand(counts)*binsize + i*binsize) 
    return phases

def randomvariate(pdf,n=1000,xmin=0,xmax=1,zero_min=True):
  """ Generate random numbers from an arbitrary distribution using the rejection
  method. See Bevington pg. 83.
    Inputs: pdf - probability distribution function from which you want to generate random numbers
            n - number of random values to output
            xmin, xmax  - range of random numbers
            zero_min - (hard to explain)
    Output: array of random values drawn from input PDF
  """

  # Calculate the minimal and maximum values of the PDF in the desired interval. 
  x = np.linspace(xmin,xmax,1000)  
  y = pdf(x)  
  pmin = 0 if zero_min else y.min()
  pmax = y.max()  

  x = np.random.uniform(xmin,xmax,n)
  y = np.random.uniform(pmin,pmax,n)

  reject = True 
  while np.any(reject):
    reject = y>pdf(x)
    x[reject] = np.random.uniform(xmin,xmax,len(reject))  
    y[reject] = np.random.uniform(pmin,pmax,len(reject)) 

  return x

def randomvariate_old(pdf,n=1000,xmin=0,xmax=1):  
  """  
  Rejection method for random number generation  
  ===============================================  
  Uses the rejection method for generating random numbers derived from an arbitrary   
  probability distribution. For reference, see Bevington's book, page 84. Based on  
  rejection*.py.  
    
  Usage:  
  >>> randomvariate(P,N,xmin,xmax)  
   where  
   P : probability distribution function from which you want to generate random numbers  
   N : desired number of random values  
   xmin,xmax : range of random numbers desired  
     
  Returns:   
   the sequence (ran,ntrials) where  
    ran : array of shape N with the random variates that follow the input P  
    ntrials : number of trials the code needed to achieve N  
    
  Here is the algorithm:  
  - generate x' in the desired range  
  - generate y' between Pmin and Pmax (Pmax is the maximal value of your pdf)  
  - if y'<P(x') accept x', otherwise reject  
  - repeat until desired number is achieved  
    
  Rodrigo Nemmen  
  Nov. 2011  
  """  
  # Calculates the minimal and maximum values of the PDF in the desired  
  # interval. The rejection method needs these values in order to work  
  # properly.  
  x=np.linspace(xmin,xmax,1000)  
  y=pdf(x)  
  #pmin=y.min()  
  pmin=0
  pmax=y.max()  
   
  # Counters  
  naccept=0  
  ntrial=0  
   
  # Keeps generating numbers until we achieve the desired n  
  ran=[] # output list of random numbers  
  while naccept<n:  
    x=np.random.uniform(xmin,xmax) # x'  
    y=np.random.uniform(pmin,pmax) # y'  
   
    if y<pdf(x):  
      ran.append(x)  
      naccept=naccept+1  
    ntrial=ntrial+1  
    
  ran=np.asarray(ran)  
    
  return ran,ntrial  

def h_test(phases, max_harmonic=20):
    """Apply the H test for uniformity on [0,1).
    The H test is an extension of the Z_m^2 or Rayleigh tests for
    uniformity on the circle. These tests estimate the Fourier coefficients
    of the distribution and compare them with the values predicted for
    a uniform distribution, but they require the user to specify the number
    of harmonics to use. The H test automatically selects the number of
    harmonics to use based on the data. The returned statistic, H, has mean
    and standard deviation approximately 2.51, but its significance should
    be evaluated with the routine h_fpp. This is done automatically in this
    routine.
    Arguments
    ---------
    events : array-like
        events should consist of an array of values to be interpreted as
        values modulo 1. These events will be tested for statistically
        significant deviations from uniformity.
    Returns
    -------
    H : float
        The raw score. Larger numbers indicate more non-uniformity.
    M : int
        The number of harmonics that give the most significant deviation
        from uniformity.
    fpp : float
        The probability of an H score this large arising from sampling a
        uniform distribution.
    Reference
    ---------
    de Jager, O. C., Swanepoel, J. W. H, and Raubenheimer, B. C., "A
    powerful test for weak periodic signals of unknown light curve shape
    in sparse data", Astron. Astrophys. 221, 180-190, 1989.
    
    Updated false alarm rate  to match Jager, Busching 2010
    """
    ev = np.reshape(phases, (-1,))
    cs = np.sum(np.exp(2.j*np.pi*np.arange(1,max_harmonic+1)*ev[:,None]),axis=0)/len(ev)
    Zm2 = 2*len(ev)*np.cumsum(np.abs(cs)**2)
    Hcand = (Zm2 - 4*np.arange(1,max_harmonic+1) + 4)
    M = np.argmax(Hcand)+1
    H = Hcand[M-1]
    fpp =np.exp(-0.4*H) 
    return (H, M, fpp)

def h_test_obs(fits_fn, par_fn):
    '''Given a fits file name, and a par filename, will return the H-score 
       and false alarm probability. 
    '''
    par = read_parfile(par_fn)
    times = fits2times(fits_fn)
    fits=pyfits.open(fits_fn)
    phs_args = [ times, par['PEPOCH'].value, par['F0'].value ]
    
    for i in range(12):
        fdot_name = 'F' + str(i+1)
        if fdot_name in par.keys():  
            phs_args.append(par[fdot_name].value)
        else:
            phs_args.append(0.0)

    phases = calc_phs(*phs_args) 
    H, M, fpp=h_test(phases)
    
    return (H, M, fpp)
    
    
class SwiftMonError(Exception):
    """
    A generic exception to be thrown by the swiftmonitor software.
    """
    pass
