import numpy as np

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
