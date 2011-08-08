import numpy as np
import swiftmonitor.observation
from scipy import optimize, integrate
import matplotlib.pyplot as plt
import pyfits
import sys

def sw2mjd(times):
  SWREFI = 51910 # integer part of Swift reference epoch in mjd
  SWREFF = 7.4287037E-4 # fractional part of Swift reference epoch in mjd

  mjds = SWREFI + SWREFF + times/86400.0
  return mjds

def get_phase(f, fdots, t, epoch):
  phase = f * (t - epoch) 
  n = 2
  for fdot in fdots:
    phase += fdot * (t - epoch)**n
    n += 1

  return phase 

def sine_fit(p, x):
  y = p[0]*np.sin(p[1]*x +p[2]) + p[3]
  return y

def calc_prob(f0, fdots, t, epoch, offset):
  phases = get_phase(f0, fdots, t, epoch)

  probs = prof_mod( phases - offset ) 
  print len(probs)
  print min(probs)
  print max(probs)
  prob_product = np.prod(probs)
  return "%.12E" % prob_product

def fourier_coeffs(phases, n):
  alphas = []
  betas = []
  for i in range(n):
    k = i+1
    alphas.append( np.sum( np.cos(k * phases) ) / len(phases) )
    betas.append( np.sum( np.sin(k * phases) ) / len(phases) )
  return alphas, betas 

def profile_from_coeffs(phases, alphas, betas, n):
  f = 1
  for i in range(n):
    k = i+1
    f += 2 * ( alphas[i] * np.cos(k*phases) + betas[i] * np.sin(k*phases) )
  return f

def events_from_binned_profile(profile):
  profile = profile[:len(profile)//2]    
  binsize = 1.0 / len(profile[:,1])
  phases = np.array([])
  i = 0
  for counts in profile[:,1]:
    phases = np.append(phases, np.random.rand(counts)*binsize + i*binsize) 
    i += 1
  return phases
  


grand_profile = np.loadtxt('grand_prof.txt')
grand_phases = events_from_binned_profile(grand_profile)

x = grand_profile[:,0]* 2 / len(grand_profile)
y = grand_profile[:,1]
err = np.sqrt(y)

plt.errorbar(x - (1.0/len(grand_profile)),y,err,fmt='k.')
plt.step(x,y,'k',where='pre')

hist = np.histogram(grand_phases,bins=len(grand_profile[:,1]) / 2)
print hist

plt.step(hist[1][:-1],hist[0],'r',where='post')

n=32
alphas, betas = fourier_coeffs(grand_phases, n)

xplot = np.arange(0,1,0.01)
yplot = profile_from_coeffs(xplot,alphas,betas,n)

plt.figure()
plt.plot(xplot,yplot)

plt.show()

"""
x = profile[:,0]* 2 / len(profile)
y = profile[:,1]
err = np.sqrt(y)
p = [5000,3.14,0.5, 12000]

errfunc = lambda p, x, y, err: (y - sine_fit(p, x)) / err
p1, success = optimize.leastsq(errfunc, p, args=(x, y, err))

print p1

xplot = np.arange(x[0],x[-1],0.01)

#plt.errorbar(x,y,err,fmt='k.')
#plt.step(x,y,'k',where='mid')

prof_mod = lambda a: sine_fit(p1,a)
ything = prof_mod(xplot)
#norm = integrate.quad(lambda x: sine_fit(p1,x) , 0, 1)
#prof_mod = lambda x: sine_fit(p1,x) / *norm[0]
#norm2 = integrate.quad(lambda x: sine_fit(p1,x) / norm[0], 0, 1)
prof_mod = lambda a: sine_fit(p1,a) / np.max(ything)

yplot = prof_mod(xplot)

plt.plot(xplot,yplot)
plt.show()
"""

