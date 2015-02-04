import numpy as np
from scipy import optimize, integrate
import types
import pickle
import matplotlib.pyplot as plt

def find_model(model_name):
  """
  Returns a model class when given a string 'name' of a model
  """
  for objname in globals():
    obj = eval(objname)
    if type(obj) == types.ClassType and issubclass(obj,Profile_Model):
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
    if type(obj) == types.ClassType and issubclass(obj,Profile_Model):
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

def makeProfileModel(model, profile, nharm=None):
    mod_type = find_model(model)
    mod = mod_type(profile)  
    if model == 'fourier' and nharm:
        mod.fit_profile(N=nharm)
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
