import matplotlib.pyplot as plt
import matplotlib.image as mpim
from matplotlib import cm
import swiftmonitor.observation as swift_obs
import pyfits
import numpy as np
import sys


def plot_region(region_file,ax,color):

  region_file.readline()
  region_file.readline()
  region_file.readline()
  line = region_file.readline()

  line = line.split('(')
  reg_type = line[0]
  coords = line[1][:-1].split(',')
  
  if reg_type == 'circle':
    cir = plt.Circle((float(coords[0]),float(coords[1])),radius = float(coords[2]),fill=False,ec=color)
    ax.add_patch(cir)
  elif reg_type == 'annulus':
    inner_cir = plt.Circle((float(coords[0]),float(coords[1])),radius = float(coords[2]),fill=False,ec=color)
    outer_cir = plt.Circle((float(coords[0]),float(coords[1])),radius = float(coords[3]),fill=False,ec=color)
    ax.add_patch(inner_cir)
    ax.add_patch(outer_cir)
  else:
    print "You shouldn't be here!"
  
def get_countrate(ob):
  fits = pyfits.open(ob.path + ob.reg_obsfile)

  countrate = fits[1].header['NAXIS2'] / fits[1].header['EXPOSURE']
  fits.close()
  return countrate
  


def plot(ob):

  if not ob.imagefile:
    ob.extract(ob.obsroot, infile= ob.path + ob.obsfile, image=True,events=False)
    ob.imagefile = ob.obsroot + '.img'

  cenx = float(ob.centroidx)
  ceny = float(ob.centroidy)

  fits = pyfits.open(ob.path + ob.imagefile)
  image = fits[0].data
  fits.close()
  
  fig = plt.figure(figsize=(11,8.5))
  im_ax = plt.axes([0.05,0.65,0.4,0.4])
  im_map = im_ax.imshow(image,interpolation='nearest',cmap=cm.gist_yarg,origin='lower',extent=(1,1000,1,1000)) 
  im_ax.plot(cenx,ceny,"rx")

  plot_region(file(ob.path + ob.back_region),im_ax,'b')
  plot_region(file(ob.path + ob.src_region),im_ax,'g')
  
  im_ax.set_ylim((ceny - 20 , ceny + 20 ))
  im_ax.set_xlim((cenx - 50 , cenx + 50 ))

  countrate = get_countrate(ob) 
  #im_ax.set_position([0.1,0.70,0.4,0.4])
  fig.text(0.55,0.80,"Countrate: %.2f counts/s" % countrate)
  fig.text(0.55,0.85,"File: %s" % ob.reg_obsfile)
 
  if ob.spec_fit:
    spec_ax = plt.axes([0.05,0.10,0.5,0.5])
    spec_img = mpim.imread(ob.path + ob.obsroot + '_xspecfit.png')
    spec_img = np.rot90(spec_img,3)
    spec_ax.imshow(spec_img,aspect='equal',interpolation='nearest')
    
    spec_ax.yaxis.set_visible(False)
    spec_ax.xaxis.set_visible(False)
   

  plt.show()




if __name__ == '__main__':
  ob = swift_obs.load(sys.argv[1])
  plot(ob)
