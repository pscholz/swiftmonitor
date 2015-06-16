import numpy as np
import astropy.io.fits as pyfits
import sys 
import scipy.optimize

import swiftmonitor.utils as smu
import swiftmonitor.config

try:
    import xspec
    XSPEC=True 
    xspec.Xset.abund='wilm'
    xspec.Xset.xsect='vern'
    xspec.Fit.statMethod = "cstat"
except(Exception):
    print('Could not import xspec')
    XSPEC=False
        
def Hardness_ratio(fits_fn, Emin=0.70, Ecut=2.0, scope='swift'):
    '''As a quick spectral, this Function reads in a FITS file, and computes
       the X-ray hardness ratio for that observation. Defaults to Swift
       Inputs:
              fits_fn - name of FITs file
              Emin  - minimum Energy [keV]
              Ecut  - Enery cut fot Hard/soft split
              scope - which telescope's Energy conversion to use
       Returns:
              HR - Hardness ratio = (HardCounts)/(SoftCounts)
              
                            
    '''
    Chan_cut=smu.energy2chan(Ecut, scope=scope) 
    t, E = smu.fits2times(fits_fn,scope=scope,Emin=Emin, Emax=None, give_t_E=True)
    SC = sum(E<Chan_cut)
    HC = sum(E>Chan_cut)
    HR = float(HC)/float(SC)
    HRerr = np.sqrt((HC)/(SC)**2+((HC*np.sqrt(SC))/(SC)**2)**2)
    return HR, HRerr
    
def prof_chisq(scale,prof, template):
    '''Compares a profile to a template using CHI2 minimization
      Called by prof_compare
      Inputs:
            scale - scale value
            prof  - profile from one observation
            template - template with which to compare  
    '''
    meanprofile = sum(prof)/float(len(prof))
    meantemplate = sum(template)/float(len(template))
    noisy = scipy.array(prof)
    errorbars = np.sqrt(meanprofile*scale+meantemplate)
    noisy = noisy - meanprofile
    template = template - meantemplate
    chisq = sum((  (template*abs(scale)-noisy)/errorbars  )**2)
    redchisq = chisq/float(len(template)-2.)
    return redchisq
    
def prof_compare(fits_fn,par_fn, templatenm,  Emin=0.7, Emax=10.23, scope='swift'):
    '''Compares a profile to a template using CHI2 minimization
       Inputs:
            fits_fn - FITS file to fold
            par_nm  - par file with which to fold
            templatenm - template to compare
       Output : 
               best - Minimal reduced Chi2
    '''
    template = np.loadtxt(templatenm).T[1]
    Nbins = len(template)
    profile = smu.fold_fits(fits_fn, par_fn, nbins=Nbins, scope='swift', Emin=Emin, Emax=Emax)[1]
    chi = []
    guessscale=sum(profile)/float(sum(template))
    best=100
    for j in range(0, Nbins):
        scalevalue = scipy.optimize.optimize.fmin(prof_chisq,guessscale, 
                           args=(template,np.roll(profile,j)), disp=0)
        value = prof_chisq(scalevalue, template,np.roll(profile,j))

        if best>value:
            best = value
    return best
    
class TOA:  
    '''Class to store TOAs
    '''                                                                   
    def __init__(self,toa_string):                                             
        self.toa = toa_string                                                  
        self.parse_toa()                                                       
                                                                               
    def parse_toa(self):                                                       
        split = self.toa.split()                                               
        self.fname = split[0]                                                  
        self.freq = float(split[1])                                            
        self.mjd_toa = float(split[2])                                         
        self.error = float(split[3])                                           
        self.scope = split[4]
        if len(split) > 5 and split[5] == "-pn":
            self.pn = int(split[6])
        elif len(split) > 7 and split[7] == "-pn":
            self.pn = int(split[8])
           
                                                                               
    def __str__(self):                                                         
       return self.toa                                                         
                                                                               

def read_tim_file(fname):
    '''Given a .tim file, will return an array of TOAs.
    '''
    toas = []
    toa_file = open(fname,"r")
   
    skipping = False
   
    for toa_str in toa_file.readlines():                                           
        if toa_str.startswith('END'):
            break
        if toa_str.startswith("SKIP"):
            skipping = True
        if toa_str.startswith("NOSKIP"):
            skipping = False
            continue
        if not skipping and not (toa_str.startswith('C') or toa_str.startswith('FORMAT')):
            toa = TOA(toa_str)
            toas.append(toa)
                                                                               
    toa_file.close()
    return np.array(toas)    

def Roughfit(obsname,Sourcename, obs_suffix='xwtw2po_cl_seporb.pha',NH=1.0,
             PL=2, kT=0.5, outfile='Roughfit.spec', PLOT=True):
    '''Given the observation number file of an
    observation, the Sourcename, and the name
    of an output file, this will preform a
    rough spectral fit to get a handle on the
    sources flux.'''
    spec = obsname+'/sw'+ obsname + obs_suffix
    fits = pyfits.open(spec)
    xspec.AllData(spec)
    xspec.AllData.ignore("bad")
    xspec.AllData.ignore("**-0.5 10.0-**")
    model = xspec.Model("tbabs(po+bb)")
    errorstr = '1.6 2 3 4 5'
    model.TBabs.nH = NH
    model.TBabs.nH.frozen = True
    model.powerlaw.PhoIndex.values=[PL, 0.01, -3.0, -2.0, 9.0, 10.0]
    model.bbody.kT.values=[kT, 0.01, 0.0001, 0.01, 100.0, 200.0]
    xspec.Fit.query = 'yes'
    xspec.Fit.perform()
    xspec.Fit.perform()
    xspec.Plot.xAxis = "keV"
    xspec.Plot.device = "/xw"
    xspec.Plot("data","resid")
    xspec.Plot()
    xspec.Fit.error(errorstr)
    try:
        xspec.AllModels.calcFlux("0.5 10. err 1000 90")
    except (RuntimeError, TypeError, NameError, 'Flux Command Error'):
        pass
    outf = open(outfile+'.spec','w')
    nH = xspec.AllModels(1)(1)
    outf.write('# nH: %f (%f - %f)\tChisq: %f\tDOF: %f\n' % ( nH.values[0], nH.error[0], nH.error[1],   xspec.Fit.statistic, xspec.Fit.dof) )
    outf.write('# obsid\tflux\tflux_low\tflux_high\tgamma\tgamma_low\tgamma_high\tkT\tkT_low\tkT_high\tPLnorm/kTnorm\n')

    flux = xspec.AllData(1).flux
    output = '%s\t%E\t%E\t%E\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' %\
             ( spec, flux[0], flux[1], flux[2],xspec.AllModels(1)(2).values[0], 
             xspec.AllModels(1)(2).error[0], xspec.AllModels(1)(2).error[1], xspec.AllModels(1)(4).values[0], 
             xspec.AllModels(1)(4).error[0], xspec.AllModels(1)(4).error[1], 
             xspec.AllModels(1)(3).values[0]/xspec.AllModels(1)(5).values[0] )             
  
    outf.write(output)   
    outf.close()
    xspec.AllData.clear()
    xspec.AllModels.clear() 


def Roughfit_PL(obsname,Sourcename, obs_suffix='xwtw2po_cl_seporb.pha', NH=1.0,
                PL=2.,  outfile='Roughfit'):
    '''Given the observation number file of an
     observation, the Sourcename, and the name
     of an output file, this will preform a
     rough spectral fit to get a handle on the
     sources flux.'''
    spec = obsname+'/sw'+obsname+obs_suffix

    fits = pyfits.open(spec)

    xspec.AllData(spec)
    xspec.AllData.ignore("bad")
    xspec.AllData.ignore("**-0.5 10.0-**")
    model = xspec.Model("tbabs(po)")
    errorstr = '1.6 2 3'
    model.TBabs.nH = NH
    model.TBabs.nH.frozen = True
    model.powerlaw.PhoIndex.values=[PL, 0.01, -3.0, -2.0, 9.0, 10.0]
    
    xspec.Fit.query = 'yes'
    xspec.Fit.perform()
    xspec.Fit.perform()
    xspec.Plot.xAxis = "keV"
    xspec.Plot.device = "/xw"
    xspec.Plot("data","resid")
    xspec.Plot()
    xspec.Fit.error(errorstr)
    try:
        xspec.AllModels.calcFlux("0.5 10. err 1000 90")
    except (RuntimeError, TypeError, NameError, 'Flux Command Error'):
        pass
    outf = open(outfile+'.spec','w')
    nH = xspec.AllModels(1)(1)
    outf.write('# nH: %f (%f - %f)\tChisq: %f\tDOF: %f\n' % ( nH.values[0], nH.error[0], nH.error[1],   xspec.Fit.statistic, xspec.Fit.dof) )
    outf.write('# obsid\tflux\tflux_low\tflux_high\tgamma\tgamma_low\tgamma_high\tkT\tkT_low\tkT_high\tPLnorm/kTnorm\n')

    flux = xspec.AllData(1).flux
    output = '%s\t%E\t%E\t%E\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' %\
             ( spec, flux[0], flux[1], flux[2],xspec.AllModels(1)(2).values[0], 
             xspec.AllModels(1)(2).error[0], xspec.AllModels(1)(2).error[1], 0, 
             0, 0, xspec.AllModels(1)(3).values[0] )             
  
    outf.write(output)   
    outf.close()
    xspec.AllData.clear()
    xspec.AllModels.clear() 

def read_pulsar_spectrum_file(name, specfile=swiftmonitor.config.monitor_code_path+'/pulsar_spectral_params.txt'):       
    ''' Input:
                name - Pulsar Name
        Output:
                PulsarName, NH, PLindex, kT
    '''
    psfile=np.loadtxt(specfile, dtype='S').T
    pulsarloc=np.where(psfile[0]==name)
    if len(pulsarloc)==0:
        print("pulsar not found, using dummy params")
        return ['dummy', 1, 1, 0]
    else:    
        return psfile.T[pulsarloc]
         
