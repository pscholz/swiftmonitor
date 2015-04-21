import numpy as np
from scipy import optimize, integrate, stats
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import sys
import os.path
import psr_utils
from psr_constants import SECPERDAY
import swiftmonitor.utils as smu
from swiftmonitor import model_profile
import time

sys.setrecursionlimit(100000)

calcprobtime = 0
logsumtime = 0
integratetime = 0


def logsumexp(array):
    """
    Recursive algorithm that sums numbers using the log of the exponentials of the input array.
        This is used because the probabilities can be very small in the likelihoods.
        Unfortunately its pretty slow.
    """
    if len(array) == 2:
        return np.logaddexp(array[0],array[1])
    else:
        return np.logaddexp(array[0],logsumexp(array[1:]))
        
def calc_prob(phases, offset, prof_mod):
    """
    Calculates the probability of an offset given a profile model and a list of phases.
    """
    probs = prof_mod( phases - offset ) 

    loglike = np.sum( np.log(probs) )
    return loglike

class PSRpar:
    """
    This class contains all the relevant information for 
        extracting a TOA, which is extracted from a .par file.
        Currently, this supports up to 12 frequency derivatives
        but NO GLITCH PARAMETERS (yet).
    """
    def __init__(self,parfile):
        pars = smu.read_parfile(parfile)
        self.epoch = pars['PEPOCH'].value
        self.f0 = pars['F0'].value

        self.fdots = np.zeros(12)
        for i in range(12):
            fdot_name = 'F' + str(i+1)
            if fdot_name in pars.keys():  
                self.fdots[i] = pars[fdot_name].value

    def __repr__(self):
        return repr((self.f0, self.fdots, self.epoch))    

def get_error(offsets, prob, del_off, debug=False):
    """
    Integrates from either side of the most probable offset to 34.1% 
        of area under distribution. Takes half of the distance between
        those two points as the sigma of the distribution.
    """
    maxoff = offsets[np.argmax(prob)]
    prob_centered = np.roll(prob, len(offsets)/2 - np.argmax(prob))

    cum_prob = np.cumsum(prob_centered)*del_off
    area_to_center = cum_prob[len(offsets)/2]

    lower_area = area_to_center - 0.341
    upper_area = area_to_center + 0.341

    a = np.argmin(np.abs(cum_prob - lower_area))
    b = np.argmin(np.abs(cum_prob - upper_area))

    sigma = (b - a) / 2.0 * del_off

    if debug:
        plt.errorbar(offsets[np.argmax(prob_centered)],0.5*max(prob_centered), xerr=sigma,fmt='o')
        plt.plot(offsets, prob_centered)
        plt.axvline(0.5,ls='dotted')
        plt.show()

    return maxoff,sigma

def get_error_gaussfit(offsets, prob, del_off, debug=False):
    """
    Fits a gaussian to the probability distribution 
        to get the phase-offset error (the width of the gaussian)
    """

    maxoff = offsets[np.argmax(prob)]
    prob_centered = np.roll(prob, len(offsets)/2 - np.argmax(prob))
    
    p0 = 0.05
    errfunc = lambda p, x, y: (y - stats.norm.pdf(x, loc=0.5, scale=p))
    output = optimize.leastsq(errfunc, p0, args=(offsets,prob_centered))
    sigma = output[0]
    
    if debug:
        plt.errorbar(offsets[np.argmax(prob_centered)],0.5*max(prob_centered), xerr=sigma,fmt='o')
        plt.plot(offsets, prob_centered)
        plt.plot(offsets, stats.norm.pdf(offsets, loc=0.5, scale=output[0]), "r--")
        plt.axvline(0.5,ls='dotted')
        plt.show()

    return maxoff, sigma

def sim_error(prof_mod,N_counts,phases,from_template=True, debug=False):
    """
    Estimate the error using simulations.
        Pulls N_counts random phases from the template where N_counts
        is the number of source counts from the observation.

        This is done N_sim number of times, and an offset between
        the template and each simulated profile is measured.

        The standard deviation of the offsets is returned.
    """
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
            simmed_phases = smu.randomvariate(prof_mod,N_counts) 
        else:
            folded = np.random.poisson(folded)
            simmed_phases = smu.events_from_binned_profile(folded)

        sim_offset, sim_error = calc_toa_offset(simmed_phases,prof_mod,no_err=True)
        sim_offsets.append((sim_offset+0.5) % 1.0)
        run += 1

    if debug:
        median = np.median(sim_offsets)
        mad = np.median(np.abs(sim_offsets-median))
        plt.hist(sim_offsets)
        print "Std dev offsets",np.std(sim_offsets)
        print "MAD offsets",mad
        plt.show()

    sys.stderr.write("Sim 100%% Complete \n")
    return np.std(sim_offsets)

def correct_model(phases,prof_mod):
    """
    Correct the model profile to match the pulsed fraction of the data.
        The should be used when the pulsed fraction of the source is varying
        substantially (e.g. during a magnetar outburst).
    """
    import fluxtool

    nbins = 32
    harmonics = 5
    folded = np.histogram(phases,nbins)
    uncertainties = np.sqrt(folded[0])

    # model pulsed flux from fourier components
    model_pulsed_flux = np.sqrt(np.sum(np.abs(prof_mod.comp[1:harmonics+1])**2)*2)/prof_mod.nbins
    model_total_flux = np.abs(prof_mod.comp[0]/prof_mod.nbins)
    model_pulsed_fraction = model_pulsed_flux / model_total_flux
    old_prof = prof_mod.prof_mod

    total_flux = np.mean(folded[0])
    rms_value, rms_uncertainty = fluxtool.rms_estimator(harmonics)(folded[0], uncertainties)
    pulsed_fraction = rms_value/total_flux

    if rms_value == 0:
        sys.stderr.write('PF Correction: Measured zero pulsed flux.\n\t Low S/N observation? \n\t Not correcting for PF.\n')
        return prof_mod.prof_mod, prof_mod.prof_mod, folded
        

    # correct zeroth fourier component to match pulsed fraction
    prof_mod.comp[0] = prof_mod.comp[0] * model_pulsed_fraction / pulsed_fraction
    model2_total_flux = np.abs(prof_mod.comp[0]/prof_mod.nbins)
    model2_pulsed_fraction = model_pulsed_flux / model2_total_flux

    prof_mod.recalc_profile()

    #print pulsed_fraction,model_pulsed_fraction,model2_pulsed_fraction
    return old_prof, prof_mod.prof_mod, folded


def calc_toa_offset(phases, prof_mod, sim_err=False, no_err=False, gauss_err=False, \
                    bg_counts=0, debug=False):
    """
    Calculate an offset between the observation pulse profile and the template pulse profile.
       This is done using the raw events (as phases) and a continuous model of the template
       (prof mod) which is used as a probability distribution.

       The error in the offset can be determined by integrating the resulting likelihood
       distribution or by using simulations (by setting sim_err=True). 

       The simulations use the total number of source counts, which for low S/N the contribution
       from the background can be large, so bg_counts (set to number of background counts
       expected in the source extraction region) can be used to correct for that.
    """
    global calcprobtime
    global logsumtime 
    global integratetime

    probs = []
    del_off = 0.001
    offsets = np.arange(0,1,del_off)
    offsets = np.append(offsets,1.0)

    starttime = time.time()
    for offset in offsets:
        prob =  calc_prob(phases, offset, prof_mod)  
        probs.append(prob)
    calcprobtime += time.time() - starttime

    # normalise as likelihood with logsumexp 
    starttime = time.time()
    logsum = logsumexp(probs)
    probs = probs - logsum - np.log(del_off)
    probs_norm = np.exp(np.array(probs))
    logsumtime += time.time() - starttime

    ## OR normalise as prob using numpy integration (doesn't do well with lots of events, like chandra)
    #starttime = time.time()
    #probs = np.exp(np.array(probs))
    #area = integrate.trapz(probs,dx=del_off)
    #probs_norm = probs/area
    #integratetime += time.time() - starttime

    if sim_err:
        maxoff = offsets[np.argmax(probs)]
        error = sim_error(prof_mod,len(phases)-bg_counts,phases,from_template=True, debug=debug) 
    elif no_err:
        maxoff = offsets[np.argmax(probs)]
        error = None
    elif gauss_err:
        maxoff, error = get_error_gaussfit(offsets, probs_norm, del_off, debug=debug)
    else:
        maxoff, error = get_error(offsets, probs_norm, del_off, debug=debug)
    
    return maxoff, error

def get_ml_toa(fits_fn, prof_mod, parfile, scope='swift', print_offs=None, frequency=None, epoch=None, \
               sim=False, bg_counts=0, Emin=None, Emax=None, gauss_err=False, tempo2=False, debug=False, \
               correct_pf=False, split_num=None, split_orbits=False):

    print_timings = False # if want to print summary of runtime

    fits = pyfits.open(fits_fn)
    if scope!='xte':
      Echans = fits[1].data['PI']
    else:
      Echans = fits[1].data['PHA']
    t = smu.fits2times(fits_fn)
    if (Emin and Emax):
        PI_min = smu.energy2chan(Emin, scope)
        PI_max = smu.energy2chan(Emax, scope)
        t = t[(Echans < PI_max) & (Echans > PI_min)]
    elif Emin:
        PI_min = smu.energy2chan(Emin, scope)
        t = t[(Echans > PI_min)]
    elif Emax:
        PI_max = smu.energy2chan(Emax, scope)
        t = t[(Echans < PI_max)]
    else:
        sys.stderr.write('No Energy Filter\n')

    if scope != 'chandra':
        exposure = fits[0].header['EXPOSURE']

    try:
        obsid = fits[0].header['OBS_ID']
    except KeyError:
        obsid = os.path.basename(fits_fn)

    if bg_counts < 0:
        bg_scale = -1.0*bg_counts
        bg_fits_fn = fits_fn.replace('reg','bgreg')
        bg_fits = pyfits.open(bg_fits_fn)
        bg_counts = int(bg_fits[1].header['NAXIS2'] * bg_scale)
        print 'BG Counts:',bg_counts
        bg_fits.close()
    if frequency and epoch:
        par = lambda: None
        par.epoch = epoch
        par.f0 = frequency
        par.fdots = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    else:
        par = PSRpar(parfile)

    # split times into multiple arrays if needed
    if split_orbits:
        dt = t[1:] - t[:-1]
        splits = np.where(dt > 0.0116)[0] # 1 ks in days
        if len(splits):
            ts = [ t[:splits[0]] ]
            for i in range(len(splits)-1):
                ts.append(t[splits[i]+1:splits[i+1]])
            ts.append(t[splits[-1]+1:])
        else:
            ts = np.atleast_2d(t)

    elif split_num:
        remainder = len(t) % split_num
        if remainder:
            sys.stderr.write("Warning: Number of events in %s not divisable by %d. " \
                             "Dropping last %d events.\n" % (obsid, split_num, remainder))
            ts = np.split(t[:-remainder],split_num)
        else:
            ts = np.split(t,split_num)

    else:
        ts = np.atleast_2d(t)

    if len(ts) > 1 and debug:
        plt.figure()
        for t in ts:
            nbins = int((t[-1] - t[0]) * 8640.0)
            hist = np.histogram(t,bins=nbins)
            plt.plot(hist[1][:-1],hist[0],c='b')
            plt.axvline(t[0],ls='--',c='k',lw=2)
            plt.axvline(t[-1],ls='-',c='k',lw=2)
        plt.show()
           

    for i,t in enumerate(ts):
        sys.stderr.write('Measuring TOA #%d for %s\n' % (i+1,obsid))

        phases = psr_utils.calc_phs(t, par.epoch, par.f0, par.fdots[0], par.fdots[1], par.fdots[2], par.fdots[3],
                                       par.fdots[4], par.fdots[5], par.fdots[6], par.fdots[7], par.fdots[8]) % 1.0 

        if correct_pf:
            old_model, new_model, corr_folded = correct_model(phases,prof_mod)
        maxoff, error = calc_toa_offset(phases,prof_mod.prof_mod,sim_err=sim,bg_counts=bg_counts, gauss_err=gauss_err, debug=debug)
        midtime = (t[-1]+t[0])/2.0
        p_mid = 1.0/psr_utils.calc_freq(midtime, par.epoch, par.f0, par.fdots[0], par.fdots[1], par.fdots[2], par.fdots[3],
                                        par.fdots[4], par.fdots[5], par.fdots[6], par.fdots[7], par.fdots[8]) 

        t0 = psr_utils.calc_t0(midtime, par.epoch, par.f0, par.fdots[0], par.fdots[1], par.fdots[2], par.fdots[3],
                               par.fdots[4], par.fdots[5], par.fdots[6], par.fdots[7], par.fdots[8]) 
        t0i = int(t0)
        t0f = t0 - t0i

        toaf = t0f + maxoff*p_mid / SECPERDAY
        newdays = int(np.floor(toaf))
        
 
        if tempo2:
            psr_utils.write_tempo2_toa(t0i+newdays, toaf-newdays, error*p_mid*1.0e6, 0000, 0.0, name=obsid) 
        else:
            psr_utils.write_princeton_toa(t0i+newdays, toaf-newdays, error*p_mid*1.0e6, 0000, 0.0, name=obsid) 

        if print_offs:
            offs_file = open(print_offs,'a')
            #print "\t",error*p_mid*1.0e6,"\t",exposure # this was for checking uncertainties scaling with exposure time
            offs_file.write(fits_fn + "\t" + str(maxoff) + "\t" + str(error) + "\n")
            #print obsid,"\tOffset:",maxoff,"+/-",error 
            offs_file.close()

        fits.close()

        
        #double check PF correction with measuring binned model pulsed fraction
        if correct_pf and debug:
            plt.figure()
            nbins = len(corr_folded[0])
            uncertainties = np.sqrt(corr_folded[0])
            area = np.sum(corr_folded[0],dtype='float')/nbins
            plt.step(corr_folded[1][:-1],np.roll(corr_folded[0]/area,int(1.0-maxoff*nbins)),where='mid')
            plt.errorbar(corr_folded[1][:-1],np.roll(corr_folded[0]/area,int(1.0-maxoff*nbins)),uncertainties/area,fmt='ko')
            model_x = np.linspace(0,1,100)
            plt.plot(model_x,old_model(model_x),label='uncorrected')
            plt.plot(model_x,new_model(model_x),label='corrected')
            plt.legend()
            plt.show()

    if print_timings:
        global calcprobtime
        global logsumtime 
        global integratetime

        sys.stderr.write('\tCalc Prob: %f s\n' % calcprobtime)
        sys.stderr.write('\tLog Sum: %f s\n' % logsumtime)
        sys.stderr.write('\tIntegrate Norm: %f s\n' % integratetime)
