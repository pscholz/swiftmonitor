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
       

class Sine_Model(Profile_Model):
  @classmethod
  def get_name(cls):
    return 'sine'

  def sine_fit(self, p, x):
    x = x % 1.0
    y = p[0]*np.sin(2*np.pi*x +p[1]) + p[2]
    return y

  def fit_profile(self, profile=None):

    if not profile:
      profile = self.profile

    x = profile[:,0]  / len(profile)
    y = profile[:,1]
    err = np.sqrt(y)

    N = np.sum(y)
    pfrac = ( np.max(y) - np.min(y) ) / N

    p0 = [ N * pfrac, 0.5, N * (1 - pfrac) ]

    errfunc = lambda p, x, y, err: (y - self.sine_fit(p, x)) / err
    p1, success = optimize.leastsq(errfunc, p0, args=(x, y, err))

    norm = integrate.quad(lambda x: self.sine_fit(p1,x) , 0, 1)
    prof_mod = lambda x: self.sine_fit(p1,x) / norm[0]

    self.prof_mod = prof_mod 

    return prof_mod

class Fourier_Model(Profile_Model):
  '''Given a grand Profile, this returns a 
     Probability distribution by the first N fourier coeffs'''
  @classmethod
  def get_name(cls):
    return 'fourier'


  def fit_profile(self, profile=None, N=5):

    if not profile:
      profile = self.profile

    phase = profile.T[0] / len(profile)
    grand = profile.T[1]
    phase = np.append(phase,1)
    grand = np.append(grand, grand[0])
    err = np.sqrt(grand)
    comp=np.fft.fft(grand)
    icomps=np.zeros(1000)
    p=np.linspace(0,1,1000)
    for i in range (1,N+1):
      icomps+=comp[i]*np.exp(2j * np.pi * p*i)
    icomp = lambda x: np.interp(x%1, p, icomps)

    def fourier_fit(p,x):
      x = x % 1.0
      y = p[0]*icomp(x) + p[1]
      return y

    p0 = [1,2]
    errfunc = lambda p, x, y, err: (y - fourier_fit(p, x)) / err
    p1, success = optimize.leastsq(errfunc, p0, args=(phase, grand, err))
    norm = integrate.quad(lambda x: fourier_fit(p1,x) , 0, 1)
    prof_mod = lambda x: fourier_fit(p1,x) / norm[0]
    return prof_mod

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
    return prof_mod

def makeProfileModel(model, profile):
  mod_type = find_model(model)
  mod = mod_type(profile)  
  
  return mod.fit_profile()


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
  prof_mod = makeProfileModel(options.model,profile)
  x = np.arange(0,1,0.001)
  plt.plot(x,prof_mod(x))
  plt.show()
