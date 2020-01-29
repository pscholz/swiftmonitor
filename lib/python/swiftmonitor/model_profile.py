import numpy as np
from scipy import optimize, integrate
from scipy.optimize import curve_fit
import types
import pickle
import matplotlib.pyplot as plt

verbose=False

def find_model(model_name):
  """
  Returns a model class when given a string 'name' of a model
  """
  for objname in globals():
    obj = eval(objname)
    if isinstance(obj, type) and issubclass(obj,Profile_Model):
      if obj.get_name() == model_name:
        return obj

  raise ProfileModelError("Model %s does not exist" % model_name) 

def list_models():
  """
  Lists all available models in a python list.
  """
  models = []
  for objname in globals():
    obj = eval(objname)
    if isinstance(obj, type) and issubclass(obj,Profile_Model):
      model_name = obj.get_name()
      if model_name:
        models.append(model_name)

  return models

class Profile_Model():
  @classmethod
  def get_name(cls):
    """
    Class method that returns the string name of the model.
    """
    return None

  def __init__(self, profile):
    self.profile = profile
    self.prof_mod = None

  def fit_profile(self, profile):
    """
    Fits a profile with a model.
      - takes a profile (array of counts per bin) as an argument.
      - sets self.prof_mod to fitted profile model
      - returns profile model
    """
    raise NotImplementedError("The fit_profile() function must be defined"\
                              " by a subclass of Profile_Model")

  def save(fname):
    """
    Saves the profile as a pickled object.
    """
    if self.prof_mod:
      f = open(fname, 'w')
      pickle.dump(self.prof_mod, f)
      f.close()
    else:
      raise ProfileModelError("Profile has not yet been fitted."\
                              " Cannot save.")
       
class Fourier_Model(Profile_Model):
    """
    Given a grand Profile, this returns a 
      Probability distribution calculated from the
      the first N fourier coeffs
    """
    @classmethod
    def get_name(cls):
        return 'fourier'

    def fit_profile(self, profile=None, N=5):

        if not profile:
            profile = self.profile

        self.nharm = N

        phase = profile.T[0] / len(profile)
        grand = profile.T[1]
        err = np.sqrt(grand)
        self.nbins = len(profile)

        comp = np.fft.rfft(grand)
        self.comp = comp

        icomps = np.zeros(1000,dtype='c16')
        phi = np.linspace(0,1,1000)

        for i in range (1,self.nharm+1):
            icomps+=2*comp[i]*np.exp(2j * np.pi * phi * i)
        icomps += comp[0]
        icomps /= len(grand)
        icomps = np.abs(icomps)
        icomp = lambda x: np.interp(x%1, phi, icomps)

        norm = integrate.quad(icomp, 0, 1)
        prof_mod = lambda x: np.interp(x%1, phi, icomps) / norm[0]
        
        #plt.step(phase,grand,where='mid')
        #plt.errorbar(phase,grand,err,fmt='bo')
        #plt.plot(phi,icomp(phi))
        #plt.show()

        self.prof_mod = prof_mod

        return prof_mod

    def recalc_profile(self):
        """
        Recalc the model after changing the fourier components (self.comp)
        """

        comp = self.comp
        icomps = np.zeros(1000,dtype='c16')
        phi = np.linspace(0,1,1000)

        for i in range (1,self.nharm+1):
            icomps+=2*comp[i]*np.exp(2j * np.pi * phi * i)
        icomps += comp[0]
        icomps /= self.nbins
        icomps = np.abs(icomps)
        icomp = lambda x: np.interp(x%1, phi, icomps)

        norm = integrate.quad(icomp, 0, 1)
        prof_mod = lambda x: np.interp(x%1, phi, icomps) / norm[0]

        self.prof_mod = prof_mod


class Gaussian_Model(Profile_Model):
    """
    Given a grand Profile, this returns a 
      Probability distribution calculated 
      from the sum of N Gaussians.
    """

    @classmethod
    def get_name(cls):
        return 'gaussian'

    def fit_profile(self, profile=None, N=1):

        if not profile:
            profile = self.profile

        self.ngauss = N


        phase = profile.T[0] 
        self.nbins = len(phase)
        phase /= float(self.nbins)
        print('phase = ', phase)
        print('nbins = ', self.nbins)
        grand = profile.T[1]
        err = np.sqrt(grand)

        # Store Fourier components for later possible use
        #comp = np.fft.rfft(grand)
        #self.comp = comp


        # Roll the profile so the peak in centred
        n_roll_bins = int(self.nbins/2) - np.argmax(grand)
        profile_centre = np.roll(grand, n_roll_bins)


        # Make successive Gaussian parameter guesses
        best_guess = []
        profile_resid = profile_centre.copy()
        for i_gauss in np.arange(self.ngauss):
            prof_mean = np.mean(profile_resid)
            # Guess A from height of peak
            A_guess = np.max(profile_resid) - prof_mean
            # Guess mean from index of peak
            mu_guess = phase[np.argmax(profile_resid)]
            # Guess width from approx. width based on nearest point to sigma 
            # width, based on FWHM
            sig_guess = np.abs((phase[np.argmax(profile_resid)]- \
                phase[np.argmin(np.abs(0.5*A_guess + prof_mean - \
                profile_resid))]))/np.sqrt(2.0*np.log(2))
            base_guess = prof_mean
            #guess = [A_guess, mu_guess, sig_guess, base_guess]
            #A, mu, sig, base = fitgauss(bins, profile_resid, yerr=np.sqrt(profile_centre), p0=guess)
            if (sig_guess > 0.): 
                best_guess.append(A_guess)
                best_guess.append(mu_guess)
                best_guess.append(sig_guess)
                if(verbose):
                    print('A, mu, sig, base (guesses): ', A_guess, mu_guess, \
                        sig_guess, base_guess) 
            else:
                self.ngauss -= 1
                    #print 'A, mu, sig, base (fit):     ', A, mu, sig, base
            print('base_guess = ', base_guess)
            print('base type  = ', np.shape(base_guess))
            profile_resid -= self.gaussian(phase, A_guess, mu_guess, 
                sig_guess, prof_mean)

        # Guess base value from mean of profile
        base_guess = np.mean(profile_centre)
        best_guess.append(base_guess)
        if(verbose):
            print('best_guess = ', best_guess)

        # Now do proper fit with ngauss Gaussian components, now that I have 
        # the best automated guess vector
        A_fit, mu_fit, sig_fit, base_fit = self.fitgauss_multi(phase, 
            profile_centre, yerr=np.sqrt(profile_centre), p0=best_guess)

        if(verbose):
            print('Fit parameters:')
            print('')
            print('    A: ', A_fit)
            print('   mu: ', mu_fit)
            print('  sig: ', sig_fit)
            print(' base: ', base_fit)

        # Roll back the fit profile model vi fit mean, using phases
        ### mu_fit -= float(n_roll_bins)/float(self.nbins)
        print('rolled_bins = ', mu_fit*self.nbins)
        print('   mu: ', mu_fit)


        #n_fit_bins=1000
        #x_fit = np.linspace(0, 1, n_fit_bins)#, endpoint=False)
        
#        prof_fit = gaussian_multi(x_fit, A_fit, mu_fit, sig_fit, base_fit)
        prof_norm = lambda x: np.sum(self.gaussian_multi(x%1, A_fit, mu_fit, 
            sig_fit, base_fit)) / float(self.nbins)

        # Make a function and normalized to have unit area under profile.
        # Assume equal bin size from 0 to 1
        ### prof_mod = lambda x: self.gaussian_multi(x%1, A_fit, mu_fit, sig_fit, \
        ###     base_fit) / (np.sum(self.gaussian_multi(x%1, A_fit, mu_fit, \
        ###     sig_fit, base_fit)) / float(len(x)))
        

        # Make a function and normalized to have unit area under profile.
        # Assume equal bin size from 0 to 1.  Roll back to match orginal profile
        prof_mod = lambda x: np.roll(self.gaussian_multi(x%1, A_fit, mu_fit, \
            sig_fit, base_fit) / (np.sum(self.gaussian_multi(x%1, A_fit, \
            mu_fit, sig_fit, base_fit)) / float(len(x))), \
            int(float(-n_roll_bins)*(float(len(x))/float(self.nbins))*np.max(x)))
        #  *int(float(len(x))/float(self.nbins)))
       
        x_prof = np.linspace(0,1,1024)
        comp = np.fft.rfft(prof_mod(x_prof))
        self.comp = comp

        self.prof_mod = prof_mod

        return prof_mod


    def recalc_profile(self):
        """
        Recalc the model after changing the fourier components (self.comp).
        In the case of a Gaussian model, just inverse FFT back to the model 
        profile to get the original function but with the changed component(s)
        (usually comp[0], the DC component)
        """

        comp = self.comp
        icomps = np.zeros(1024,dtype='c16')
        phi = np.linspace(0,1,1024)
        self.nharm=len(comp)-1 # hardcoded for now -- use max number of components

        for i in range (1,self.nharm+1):
            icomps+=2*comp[i]*np.exp(2j * np.pi * phi * i)
        icomps += comp[0]
        icomps /= self.nbins
        icomps = np.abs(icomps)
        icomp = lambda x: np.interp(x%1, phi, icomps)

        norm = integrate.quad(icomp, 0, 1)
        prof_mod = lambda x: np.interp(x%1, phi, icomps) / norm[0]

        self.prof_mod = prof_mod


    # Perform a fit of an x-y data set to a MULTIPLE gaussian profile.
    # Add in a chi-sq calculation as well to return back to the parent 
    # routine.
    def fitgauss_multi(self, xdata, ydata, yerr=None, p0=None):
        sig = None
        if(yerr!=None):
            sig = 1.0/yerr**2.0
        if (verbose):
            print('p0 IN = ', p0)
            print('\nxdata = ', xdata)
            print('ydata = ', ydata)
            print('yerr = ', sig)
            print('')
        popt, pcov = curve_fit(self.gaussian_multi_func, xdata, ydata, 
            sigma=sig, p0=p0)
        if (verbose):
            print('popt = ', popt)
            print('pcov = ', pcov)
        ##print 'n_popt = ', len(popt)
        ngauss = (len(popt) - 1)/3
        print('ngauss = ', ngauss)
#        A = [np.zeros(ngauss, dtype=float)]
#        mean = np.zeros_like(A)
#        width = np.zeros_like(A)
        A = []
        mean = []
        width = []
##        for i_gauss in np.arange(ngauss):
            #print i_gauss*(ngauss+1)
##            A[i_gauss]     = popt[0 + (i_gauss)*(ngauss+1)]
##            mean[i_gauss]  = popt[1 + (i_gauss)*(ngauss+1)]
##            width[i_gauss] = popt[2 + (i_gauss)*(ngauss+1)]
        for i_gauss in np.arange(0, len(popt)-1, 3):
            A.append(popt[i_gauss])
            mean.append(popt[i_gauss + 1])
            width.append(popt[i_gauss + 2])
        base = popt[-1]

        A = np.array(A)
        mean = np.array(mean)
        width = np.array(width)
        return A, mean, width, base

    # Multiple gaussians.  NEED to pass A, mu, and sigma as arrays, 
    # even if length 1.  Base is a scalar float
    def gaussian_multi_func(self, x, *params):
        # Check that A, mu, and sigma have the same array length
        gauss_out = np.zeros_like(x)
        ngauss = (len(params) - 1)/3
        #print 'ngauss = ', ngauss
        ######if(len(A) == len(mu) == len(sigma) == ngauss):  
        #     A, mu, sigma = p
###        for i_gauss in np.arange(ngauss):
###            A   = params[0 + i_gauss*(ngauss)]
###            mu  = params[1 + i_gauss*(ngauss)]
###            sig = params[2 + i_gauss*(ngauss)]
        for i_gauss in np.arange(0, len(params)-1, 3):
            #print 'i_gauss =', i_gauss, '   ngauss = ', ngauss
            A   = params[i_gauss]
            mu  = params[i_gauss + 1]
            sig = params[i_gauss + 2]
            gauss_out += A*np.exp(-(x-mu)**2/(2.0*sig**2))

        # last number is params sequence is the baseline value
        gauss_out += params[-1] 
        return gauss_out

         
    # Multiple gaussians.  NEED to pass A, mu, and sigma as arrays, 
    # even if length 1.  Base is a scalar float
    def gaussian_multi(self, x, A, mu, sigma, base):
        # Check that A, mu, and sigma have the same array length
        gauss_out = np.zeros_like(x)
        ngauss = len(A)

        if(len(A) == len(mu) == len(sigma) == ngauss):  
        #     A, mu, sigma = p
            for i_gauss in np.arange(ngauss):
              gauss_out += A[i_gauss]*np.exp(-(x-mu[i_gauss])**2/ \
                (2*sigma[i_gauss]**2))

            gauss_out += base
        else:
            print('Error: A, mu, and sigma, need to be of same array ')
            print('length (= ngauss)')

        return gauss_out

    # Calculate a Gaussian given parameters
    def gaussian(self, x, A, mu, sigma, base):
    #     A, mu, sigma = p
        return base + A*np.exp(-(x-mu)**2/(2*sigma**2))





class Sine_Model(Fourier_Model):
    @classmethod
    def get_name(cls):
        return 'sine'

    def fit_profile(self, profile=None):
        return Fourier_Model.fit_profile(self,profile=profile,N=1) 

class Lin_Interp(Profile_Model):
  @classmethod
  def get_name(cls):
    return 'lin_interp'

  def fit_profile(self, profile=None):
    '''Given a grand Profile, this returns a 
    Probability distribution by linear interpolation'''

    if not profile:
      profile = self.profile

    phase = profile.T[0] / len(profile)
    grand = profile.T[1]
    phase = np.append(phase,1)  
    grand = np.append(grand, grand[0])
    area = integrate.trapz(grand, x=phase)
    grand = grand/area
    prof_mod = lambda x: np.interp(x%1, phase, grand)
    self.prof_mod = prof_mod

    return prof_mod

def makeProfileModel(model, profile, n=None):
    mod_type = find_model(model)
    mod = mod_type(profile)  
    if ((model == 'fourier') | (model == 'gaussian')) and n:
        mod.fit_profile(N=n)
    else:
        mod.fit_profile()
    
    return mod


def ProfileModelError(Exception):
    pass

if __name__ == '__main__':

    from optparse import OptionParser

    models = list_models()

    models_str = "Available Models:  "
    for model in models:
      models_str += model + '  '

    parser = OptionParser("Usage: %prog [options] profile", epilog=models_str ,version="%prog 1.0")
    parser.add_option("-m", "--model",
          	    dest="model", type='string',
          	    help="Name of model to fit to profile.",
          	    default=None)
    parser.add_option("-l", "--linear_interp",
          	    dest="lin_interp", action='store_true',
          	    help="Use linear interpolation.",
          	    default=False)
     
    (options,args) = parser.parse_args()

    profile = np.loadtxt(args[0])
    mod = makeProfileModel(options.model,profile)
    x = np.arange(0,1,0.001)
    plt.plot(x,mod.prof_mod(x))
    plt.show()
