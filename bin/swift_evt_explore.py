#!/usr/bin/env python
# swift_evt_explore.py Written by Paul Scholz Oct 2012

import os
import tempfile
import matplotlib.pyplot as plt
import numpy as np
import pyfits
from matplotlib import cm
import matplotlib.ticker
from optparse import OptionParser
from swiftmonitor import ftools

class ImageExplorer:
  def __init__(self,ax1,ax2,evt_fn,out_fn=None,src_reg=None,back_reg=None):

    self.evt_fn = evt_fn
    fits = pyfits.open(evt_fn)
    self.data = fits['EVENTS'].data
    self.orig_gtis = list(fits['GTI'].data)
    self.gtis = self.orig_gtis
    fits.close()

    self.lcax = ax2.axes # light curve     
    self.imax = ax1.axes # image
    self.cbax = None # colorbar
    self.press = None
    self.fig = lcax.figure
    self.t1, self.t2 = None, None
    self.excludes = []

    hist, binedges = np.histogram(self.data['TIME'],bins=np.arange(np.min(self.data['TIME']),\
                                  np.max(self.data['TIME']),1))
    self.lc_x, self.lc_y = binedges[:-1], hist
    self.gti_draw = True
    self.src_reg = src_reg
    self.back_reg = back_reg
    self.out_fn = out_fn

    self.plot_lc()
    self.plot_img()
    self.print_usage()

    #connect to all the events we need
    self.cidpress = self.fig.canvas.mpl_connect(
        'button_press_event', self.on_press)
    self.cidrelease = self.fig.canvas.mpl_connect(
        'button_release_event', self.on_release)
    self.cidmotion = self.fig.canvas.mpl_connect(
        'motion_notify_event', self.on_motion)
    self.cidkeypress = self.fig.canvas.mpl_connect(
        'key_press_event', self.on_keypress)

  def print_usage(self):
    usage = "\nSwift Event File Explorer\n\n"\
            "\tSelect time intervals for which to plot image by right clicking\n"\
            "\ton start of desired time interval and dragging to end of desired\n"\
            "\ttime interval.\n\n"\
            "Hotkeys:\n"\
            "\th - print this help message.\n"\
            "\ti - toggle plotting of GTIs and excluded intervals.\n"\
            "\te - exclude selected (plotted) interval from GTIs.\n"\
            "\tc - clear exlcuded intervals and plot entire observation.\n"\
            "\tw - write out new event file with events from new GTIs.\n"\
            "\tq - quit.\n"

    print usage  

  def on_press(self, event):
    if (event.inaxes != self.lcax) or (event.button != 3): return
    if self.press:
      self.plot_lc()
    self.lcax.axvline(event.xdata,c='r')
    self.press=event.xdata
    self.fig.canvas.draw()

  def on_motion(self,event):
    return

  def on_release(self,event):
    if (event.inaxes != self.lcax) or (event.button != 3): return
    self.lcax.axvline(event.xdata,c='r')
    self.t1, self.t2 = self.press, event.xdata
    self.plot_img()
    self.fig.canvas.draw()

  def on_keypress(self,event):
    if (event.key == 'q'):
      exit(1)

    if (event.key == 'h'):
      self.print_usage()

    if (event.key == 'e'):
      exclude = ( self.t1, self.t2 )
      if exclude not in self.excludes and np.all(exclude):
        self.excludes.append(exclude)
      self.exclude_intervals()
      self.plot_lc()
      
    if (event.key == 'i'):
      self.gti_draw = not self.gti_draw
      self.plot_lc()
 
    if (event.key == 'w'):
      self.write_gtis()

    if (event.key == 'c'):
      self.t1, self.t2 = None, None 
      self.excludes = []
      self.gtis = self.orig_gtis
      self.gti_draw = True

      self.plot_lc()
      self.plot_img()
      self.fig.canvas.draw()

  def plot_lc(self):
    self.lcax.cla()
    self.lcax.plot(self.lc_x,self.lc_y,'b')
    if self.gti_draw:
      self.draw_excludes()
      self.draw_gtis()
    self.fig.canvas.draw()

  def plot_img(self):

    if self.t1 and self.t2:
      if self.t1 > self.t2:
        t1, t2 = self.t1,self.t2
        self.t1, self.t2 = t2, t1
      data_subset = self.data[(self.data['TIME'] < self.t2) &\
                              (self.data['TIME'] > self.t1)]
    else:
      data_subset = self.data

    self.imax.cla()
    if len(data_subset):
      bins = np.arange(1024)
      image,xedges,yedges = np.histogram2d(data_subset['Y'],data_subset['X'],bins)
      img = self.imax.imshow(np.log10(image),interpolation='nearest',cmap=cm.gist_yarg,\
                             origin='lower',extent=(1,1024,1,1024),aspect=1)
      src_xlim, src_ylim = plot_region(file(self.src_reg),imax,'b')
      back_xlim, back_ylim = plot_region(file(self.back_reg),imax,'g')
      xlim = [min(src_xlim[0],back_xlim[0]),max(src_xlim[1],back_xlim[1])]
      ylim = [min(src_ylim[0],back_ylim[0]),max(src_ylim[1],back_ylim[1])]
      imax.set_ylim(ylim)
      imax.set_xlim(xlim)

      log_format = matplotlib.ticker.FuncFormatter(lambda x,pos: ('%.1f' % 10**x))
      if self.cbax: 
        self.cbax.cla()
        self.fig.colorbar(img,cax=self.cbax,format=log_format)
      else:
        cb = self.fig.colorbar(img,ax=self.imax,format=log_format)
        self.cbax = cb.ax 
    else:
      print "NO EVENTS IN SELECTION"

  def draw_excludes(self,excludes=None):
    if not excludes:
      excludes = self.excludes
    for exclude in excludes:
      self.lcax.axvspan(exclude[0],exclude[1],alpha=0.2,color='r')

  def draw_gtis(self):
    for gti in self.gtis:
      self.lcax.axvspan(gti[0],gti[1],alpha=0.2,color='g')

  def exclude_intervals(self):
    merged = merge_intervals(self.excludes)
    self.excludes = merged
    new_gtis = remove_intervals(self.gtis,merged)
    self.gtis = new_gtis
    self.plot_lc()

  def write_gtis(self):
    fdesc, gti_fn = tempfile.mkstemp(suffix='.txt',prefix='gti')
    f = os.fdopen(fdesc,'w')
    for gti in self.gtis:
      f.write('%12.5f\t%12.5f\n' % gti) 
    f.close()
    
    out_fn = self.out_fn if self.out_fn else os.path.splitext(self.evt_fn)[0] + '_clean'
    ftools.extract(out_fn,self.evt_fn,gtifile=gti_fn)
    os.remove(gti_fn)
    print "\nExtraction done. Output event file:\n\t%s.evt" % out_fn 
    
def merge_intervals(intervals):
  if not intervals:
    return []
  intervals.sort()
  result = []
  (a, b) = intervals[0]
  for (x, y) in intervals[1:]:
    if x <= b:
      b = max(b, y)
    else:
      result.append((a, b))
      (a, b) = (x, y)
  result.append((a, b))
  return result

def remove_intervals(intervals,to_remove):
  if not to_remove:
    return intervals
  intervals.sort()
  result = []
  for (a,b) in intervals:
    for (x,y) in to_remove:
      if (y < b) and (y > a):
        if (x > a):
          result.append((a,x))
        (a,b) = (y,b)
      elif (x < b) and (x > a):
        (a,b) = (a,x)
      elif (x < a) and (y > b):
        break
    else:
      result.append((a,b))

  return result

def plot_region(region_file,ax,color):

  region_file.readline()
  region_file.readline()
  region_file.readline()
  line = region_file.readline()

  line = line.split('(')
  reg_type = line[0]
  coords = line[1].rstrip().rstrip(')').split(',')
  
  if reg_type == 'circle':
    x, y, r = float(coords[0]),float(coords[1]),float(coords[2])
    cir = plt.Circle((x,y),radius = r,fill=False,ec=color)
    ax.add_patch(cir)
    return [x-r, x+r],[y-r, y+r]
  elif reg_type == 'annulus':
    x, y, r1, r2 = float(coords[0]),float(coords[1]),float(coords[2]),float(coords[3])
    inner_cir = plt.Circle((x,y),radius = r1,fill=False,ec=color)
    outer_cir = plt.Circle((x,y),radius = r2,fill=False,ec=color)
    ax.add_patch(inner_cir)
    ax.add_patch(outer_cir)
    return [x-r2, x+r2],[y-r2,y+r2]
  else:
    print "You shouldn't be here!"

if __name__=='__main__':

  parser = OptionParser("Usage: %prog [options] input_file", version="%prog 1.0")
  parser.add_option("-s",
                  dest="src_reg", type='string',
                  help="Name of source region file. Default='source.reg'.",
                  default='source.reg')
  parser.add_option("-b",
                  dest="back_reg", type='string',
                  help="Name of background region file. Default='back.reg'.",
                  default='back.reg')
  parser.add_option("-o",
                  dest="out_root", type='string',
                  help="Root of output event file. '.evt' will be appended."\
                       "Default=input_root + '_clean'",
                  default=None)

  (options, args) = parser.parse_args()


  plt.figure(figsize=(8,16))
  
  imax = plt.subplot(211)
  lcax = plt.subplot(212)
  
  ImageExplorer(imax,lcax,args[0],out_fn=options.out_root,src_reg=options.src_reg,\
                back_reg=options.back_reg)
  
  plt.show()

