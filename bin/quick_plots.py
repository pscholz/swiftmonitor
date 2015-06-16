import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import scipy.stats 
import swiftmonitor.utils as smu
import swiftmonitor.quicklook_tools as qlt
from swiftmonitor import ml_toa, model_profile

from optparse import OptionParser

parser = OptionParser("Usage: %prog [options] fitsfile",version="%prog 1.0")
parser.add_option("-p", "--par",
    dest="parfile", type='string',
    help="Name of par file with pulsar ephemeris.",
    default=None)
parser.add_option("-n", "--name",
    dest="name", type='string',
    help="Name of  pulsar.",
    default=None)    
parser.add_option("--Emin",
    dest="emin", type='float',
    help="Minimum energy of events to use (works only for Swift, Nustar).",
    default=0.7)
parser.add_option("--Emax",
    dest="emax", type='float',
    help="Minimum energy of events to use (works only for Swift, Nustar).",
    default=10.0)
parser.add_option("--Ecut",
    dest="ecut", type='float',
    help="Maximum energy of events to use (works only for Swift, Nustar).",
     default=2.0)		  		    
parser.add_option("-l", "--list",
    dest="list", type='string',
    help="File with list of event files to get events from.",
    default=None) 
parser.add_option("-f", "--profile",
    dest="profile", type='string',
    help="Name of profile to fit model to.",
    default=None)
parser.add_option("--scope",
    dest="scope", type='string',
    help="Event files are from this telescope, default swift.",
    default='swift')
parser.add_option("--num-harmonics",
    dest="num_harmonics", type='int',
    help="Number of fourier components to use in fourier model. Default is 2.",
		  default=2)
parser.add_option("--writefile",
    dest="writefile", type='string',
    help="Print TOA to file (ONLY FOR tempo2 FORMAT)",
    default='toa.tim')
parser.add_option("--correct-pf",
    dest="correct_pf", action='store_true',
    help="Rescale the template profile to have the same pulsed " \
               + "fraction as the observation for which the TOA is being calculated.",
    default=False)		
parser.add_option( "--quicklookdataname",
    dest="quicklookdataname", type='string',
    help="Name of file to store quicklook results.",
    default='quicklookdata.qld')     
parser.add_option( "--fluxstorage",
    dest="fluxstorage", type='string',
    help="Name of file to store Flux results.",
    default='Roughfit')  		  
(options,args) = parser.parse_args()

profile = np.loadtxt(options.profile)
prof_mod = model_profile.makeProfileModel('fourier', profile, nharm=options.num_harmonics)
specparams = qlt.read_pulsar_spectrum_file(options.name)[0]#[name, NH, Gamma, kT]    
freq = smu.read_parfile(options.parfile)['F0'].value

if options.list:
    flist = np.loadtxt(options.list,dtype='S')
    plt.figure(figsize=[8.5, 11])
    plt.subplots_adjust(wspace=0.05, hspace = 0.02, left=0.15, bottom=0.1, right=0.95, top=0.9)
    for obs in flist:
 
        fits_fn = obs+'/sw'+obs+'xwtw2po_cl_bary_reg.evt'
        QLDs=glob.glob(obs+'/'+options.quicklookdataname) 
        if len(QLDs)==0:
            try:
                HR, HRerr = qlt.Hardness_ratio(fits_fn, Emin=options.emin, Ecut=options.ecut,
                               scope=options.scope)
                redchi = qlt.prof_compare(fits_fn,options.parfile, options.profile,  Emin=options.emin, Emax=options.emax,
                              scope=options.scope)
                outf = open(obs+'/'+options.quicklookdataname,'w')
                outf.write('#obsid\tHR\tHRerr\tRedChi\n')
                outf.write('%s\t%E\t%E\t%E\n' %(obs, HR, HRerr, redchi))              
            except (Exception):
                print('MANUALLY CHECK HR, Chi ', obs)                   
                           
              
        SPECs=glob.glob(obs+'/'+options.fluxstorage+'.spec') 
        if len(SPECs)==0:    
            try:
                if specparams[3]=='0':
                    qlt.Roughfit_PL(obs,specparams[0], obs_suffix='xwtw2po_cl_seporb.pha',
                            NH=float(specparams[1]), PL=float(specparams[2]),  outfile=obs+'/'+options.fluxstorage)
                else:
                    qlt.Roughfit(obs,specparams[0],  NH=float(specparams[1]), PL=float(specparams[2]), kT=float(specparams[3]), outfile=obs+'/'+options.fluxstorage)        
            except (Exception):
                print('MANUALLY CHECK SPECTRUM ', obs)                 
            
        TIMs=glob.glob(obs+'/toa.tim') 
        if len(TIMs)==0:  
             try:      
                #os.system('python ~/swiftmonitor/bin/get_ML_TOA.py --tempo2 -m fourier -n 5 -f '
                #         +options.profile+' -p '+options.parfile + ' '+ fits_fn + ' > '+ obs+ '/toa.tim')
                ml_toa.get_ml_toa(fits_fn, prof_mod, options.parfile, scope=options.scope, \
                      Emin=options.emin, Emax=options.emax,  tempo2=True, \
                      correct_pf=options.correct_pf, writefile=obs+'/'+options.writefile)
             except (Exception):
                print('MANUALLY CHECK TOA', obs)
        BPs=glob.glob(obs+'/32nd_s.burst') 
        if len(BPs)==0:  
             try:      
                os.system('python ~/Scripts/swburst_find.py -b 0.03125 -o '+obs+'/32nd_s.burst ' +fits_fn )
             except (Exception):
                print('MANUALLY CHECK BURST ', obs)           
        else:      
            try:              
                HR, HRerr, redchi = np.loadtxt(obs+'/'+options.quicklookdataname, usecols=(1,2,3))
                rawtoa = qlt.read_tim_file(obs+'/toa.tim')
                toa=rawtoa[0].mjd_toa
                toa_err=rawtoa[0].error/1E6
                
                F, Flo, Fhi= np.loadtxt(obs+'/'+options.fluxstorage+'.spec', usecols=(1,2,3))*1E11
                BM=np.loadtxt(obs+'/32nd_s.burst').T
                ax1=plt.subplot(511)
                plt.errorbar(toa, F, yerr=[(F-Flo, Fhi-F)],fmt='ko' )
              
                plt.ylabel('0.5-10 keV Flux\n [$10^{-11}$erg s$^{-1}$cm$^{-2}$]')
                ax2=plt.subplot(512,sharex=ax1)
              
                plt.errorbar(toa, redchi, fmt='ko', label=None )
                plt.ylabel('Profile $\chi^2_\\nu$')
                
                ax3=plt.subplot(513,sharex=ax1)
            
                plt.errorbar(toa, HR, yerr=HRerr,fmt='ko' )
                plt.ylabel('Hardness Ratio')
                
                ax4=plt.subplot(514,sharex=ax1)
                plt.errorbar(toa, (smu.times2phases(toa, options.parfile)-0.5)%1, yerr=toa_err*freq,fmt='ko')
                plt.ylabel('Timing Residuals')
            
                ax5=plt.subplot(515,sharex=ax1)
                if len(BM)!=0 :
                    if len(np.shape(BM))==1:
                        plt.errorbar(toa, -np.log10(BM[2]), fmt='ko' )
                    else:
                        plt.errorbar(np.tile(toa, len(BM[2])), -np.log10(BM[2]), fmt='ko' )
                    plt.ylabel('Burst Occurance\n [1 in 10$^{X}$]')
                
            except (Exception):
                print('MANUALLY CHECK ', obs)  
    for ax in (ax1, ax2, ax3, ax4):
        plt.setp(ax.get_xticklabels(), visible=False)     
    Nbins=len(profile.T[0])
    ax2.axhline(scipy.stats.chi2.interval(0.997, Nbins)[1]/Nbins, ls='dashed', color='r' , label='$3-\\sigma$')
    ax2.axhline(scipy.stats.chi2.interval(0.95, Nbins)[1]/Nbins, ls='dotted', color='r', label='$2-\\sigma$')
    ax2.legend(loc=2)
    ax5.axhline(2.5228, ls='dashed', color='r' , label='$3-\\sigma$') 
    ax5.legend(loc=2)
    plt.xlabel('Modified Julian Date')
    plt.savefig('test.png')        

        
