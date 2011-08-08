import os, sys, time
import numpy as np
import subprocess, re, pickle
import pyfits
import swiftmonitor.config

class TableDownError(Exception):
  def __str__(self):
    return 'swiftxrlog table is down.'


def timed_execute(cmd): 
  """
  Execute the command 'cmd' after logging the command
    to STDOUT.  Return the wall-clock amount of time
    the command took to execute.
  """
  sys.stdout.write("\n'"+cmd+"'\n")
  sys.stdout.flush()
  start = time.time()
  os.system(cmd)
  end = time.time()
  return end - start

def load(pickle_file):
  return pickle.load(file(pickle_file))

class Observation:
  
  def __init__(self,obsid,path,mode='wt',pulsar=None):
    self.obsid = obsid
    self.path = path
    self.mode = mode
    self.pulsar = pulsar
    self.baryfile = None
    self.imagefile = None
    self.centroidx = None
    self.centroidy = None
    self.spec_fit = None
  
    if self.pulsar:
      self.ra, self.dec = self._get_pulsar_position()
 
    if path[-1] != '/':
      self.path = path + '/'

  def _get_pulsar_position(self):
    """
    Private method that reads the coordinates of a pulsar from 
      the pulsar_coords.txt file. This position is used in barycentering
      and when extracting the default region. Returns the RA and Dec in 
      decimal degrees format.
    """
    f = open(swiftmonitor.config.monitor_code_path + "pulsar_coords.txt",'r')

    ra, dec = None, None
    for line in f:
      line_split = line.split() 
      if line_split[0] == self.pulsar:
        ra, dec = line_split[1], line_split[2]
        break

    if ra and dec:
      split_ra = ra.split(':')
      split_dec = dec.split(':')
      dec_sign = float(split_dec[0]) / abs(float(split_dec[0]))
      ra_deg = ( float(split_ra[0]) + float(split_ra[1])/60.0 + float(split_ra[2])/3600.0 ) * 15.0  
      dec_deg = abs(float(split_dec[0])) + float(split_dec[1])/60.0 + float(split_dec[2])/3600.0   
      dec_deg *= dec_sign
    else:
      print "Pulsar %s not in coords file." % self.pulsar
      ra_deg, dec_deg = ra, dec

    f.close()
    return ra_deg, dec_deg


  def _query_heasarc(self):
    """
    Query heasarc archive for the location and names of the cleaned 
      observation file and the orbit file. Returns the filename of 
      the observation, the filename of the orbit file, and the 
      path to the obsid directory.
    """

    obsid = self.obsid

    if self.mode == 'wt':
      cmd = subprocess.Popen(['browse_extract_wget.pl', 'table=swiftxrlog', 'position=none', 'operation_mode=WINDOWED', 
                              'pointing_mode=pointing', 'param=obsid,' + obsid ],stdout=subprocess.PIPE)
      point_re = re.compile("^\|.*WINDOWED\|pointing(?!_mode)")
    elif self.mode == 'pc':
      cmd = subprocess.Popen(['browse_extract_wget.pl', 'table=swiftxrlog', 'position=none', 'operation_mode=PHOTON', 
                              'pointing_mode=pointing', 'param=obsid,' + obsid ],stdout=subprocess.PIPE)
      point_re = re.compile("^\|.*PHOTON\|pointing(?!_mode)")
  
    date = None
    window = None

    tabledown_re = re.compile('Tableswiftxrlogdoesnotseemtoexist!')

    for line in cmd.stdout:
      line = line.replace(' ','')
      if tabledown_re.search(line):
        raise TableDownError
      if point_re.match(line):
	split_line = line.split('|')
	new_date = split_line[5][:10]
	new_window = split_line[7]
	if new_date != date and date != None:
	  print "Not all on same date."
	if new_window != window and window != None:
	  print "More than one windowsize!"
	date = new_date
	window = new_window

    cmd.stdout.close()

    if self.mode == 'wt':
      if window == '100': 
	size='1' 
      elif window == '200':
	size='2' 
      elif window == '300':
	size='3' 
      elif window == '400': 
	size='4' 
      elif window == '500': 
	size='5' 
      else: 
	size='0' 

      obsfile = "sw" + obsid + "xwtw" + size + "po_cl.evt.gz"

    elif self.mode == 'pc':
      if window == '490x490':
	size='1' 
      elif window == '500x500':
	size='2' 
      elif window == '600x600':
	size='3' 
      elif window == '480x480':
	size='4' 
      else: 
	size='0'

      obsfile = "sw" + obsid + "xpcw" + size + "po_cl.evt.gz"

    date = date.split('-')
    date = date[0] + "_" + date[1]

    orbitfile = "sw" + obsid + "sao.fits.gz"

    obsdir_url = "http://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/" + date + "//" + obsid + "/"

    return obsfile, orbitfile, obsdir_url

  def download(self):
    """
    Download the observation. 
    """
    print "Querying HEASARC...\n"
  
    obsfile, orbitfile, obsurl = self._query_heasarc()

    print "Downloading observation...\n"

    cmd = 'wget -q -N -O ' + self.path + '/' + obsfile + " " + obsurl + '/xrt/event/' + obsfile
    wget_time = timed_execute(cmd)
    cmd = 'wget -q -N -O ' + self.path + '/' + orbitfile + " " + obsurl + '/auxil/' + orbitfile
    wget_time += timed_execute(cmd)

    self.orbitfile = orbitfile
    self.obsfile = obsfile
    self.obsroot = obsfile.split('.')[0]

  def manual_add(self, obsfile, orbitfile):
    """
    Set up for an observation for which the orbitfile and obsfile were downloaded manually. 
      Result is the same as download() without downloading the files or querying HEASARC.
    """
    self.orbitfile = orbitfile
    self.obsfile = obsfile
    self.obsroot = obsfile.split('.')[0]

  def barycentre(self, RA=None, Dec=None):
    """
    Barycentre the observation. If not provided a RA and Dec, will use pulsar RA and Dec
      if a pulsar is defined for the Observation object, otherwise the RA and Dec from
      the fits header will be used.
    """
    print "Barycentreing observation...\n"

    self.baryfile = self.obsroot + '_bary.evt'
    if RA and Dec:
      cmd = 'barycorr infile=%s/%s outfile=%s/%s orbitfiles=%s/%s ra=%s dec=%s clobber=yes' %\
            (self.path, self.obsfile, self.path, self.baryfile, self.path, self.orbitfile, RA, Dec)
    elif self.pulsar:
      cmd = 'barycorr infile=%s/%s outfile=%s/%s orbitfiles=%s/%s ra=%s dec=%s clobber=yes' %\
            (self.path, self.obsfile, self.path, self.baryfile, self.path, self.orbitfile, self.ra, self.dec)
    else:
      print "No RA and Dec given. Using RA and Dec of target in fits header..."
      cmd = 'barycorr infile=%s/%s outfile=%s/%s orbitfiles=%s/%s clobber=yes' %\
            (self.path, self.obsfile, self.path, self.baryfile, self.path, self.orbitfile)
      
    bary_time = timed_execute(cmd)

  def extract(self,outroot,infile=None,events=True,image=False,pha=False,lc=False,region=None,chanlow=0,chanhigh=1023):
    """
    Wrapper for extractor ftool. If infile is None will use baryfile or obsfile as input.
   
      Arguments:
        - outroot: root of the output files. Will attach extension depending on
                   type of output (.img, .evt, etc.)

      Optional Arguments:
        - infile: Input file to extract from.
                  Default=None.
        - events: Boolean of whether or not to extract events file. 
                  Default=True.
        - image: Boolean of whether or not to extract image file. 
                 Default=False.
        - pha: Boolean of whether or not to extract spectrum. 
               Default=False.
        - lc: Boolean of whether or not to extract binned light curve. 
              Default=False. 
        - region: Filename (with path) of region file to extract from region. 
                  Default=False.
       
    """

    if infile == None:
      if self.baryfile == None:
        print "Using obsfile as input."
        infile = self.path + self.obsfile
      else:
        print "Using baryfile as input."
        infile = self.path + self.baryfile
   
    args = '[%d:%d] xcolf=X ycolf=Y tcol=TIME ecol=PI timefile=NONE xcolh=X ycolh=Y ' % (chanlow, chanhigh)

    if image == True:
      args += 'imgfile=%s%s.img ' % (self.path, outroot)
    else:
      args += 'imgfile=NONE '
    if pha == True:
      args += 'phafile=%s%s.pha ' % (self.path, outroot)
    else:
      args += 'phafile=NONE '
    if lc == True:
      args += 'fitsbinlc=%s%s.lc ' % (self.path, outroot)
    else:
      args += 'fitsbinlc=NONE '
    if events == True:
      args += 'eventsout=%s%s.evt ' % (self.path, outroot)
    else:
      args += 'eventsout=NONE '
    if region != None:
      args += 'regionfile=%s ' % region
    else:
      args += 'regionfile=NONE '

    args += 'clobber=yes'
    cmd = 'extractor ' + infile + args 
    extract_time = timed_execute(cmd)

  def find_centroid(self,force_redo=False,use_max=True):
    """
    Finds centroid of source (hopefully)
      Returns x, y coordinates of centroid in pixels.
    """

    if not self.imagefile:
      self.extract(self.obsroot, infile= self.path + self.obsfile, image=True,events=False)
      self.imagefile = self.obsroot + '.img'

    if use_max:
      fits = pyfits.open(self.path + self.imagefile)
      image = fits[0].data
      fits.close()
    
      y,x = np.unravel_index(np.argmax(image),np.shape(image))
      x += 1
      y += 1

    else:

      if self.centroidx != None and self.centroidy != None and force_redo == False:
	return self.centroidx, self.centroidy


      if self.pulsar:
	cmd = subprocess.Popen(['ximage', '@/homes/borgii/pscholz/bin/swiftmonitor/wt_centroid_radec.xco',\
				  self.path + self.imagefile, str(self.ra), str(self.dec) ], stdout=subprocess.PIPE) 
      else:  
	cmd = subprocess.Popen(['ximage', '@/homes/borgii/pscholz/bin/swiftmonitor/wt_centroid.xco',\
				  self.path + self.imagefile], stdout=subprocess.PIPE) 

      region_re = re.compile('^[ ]*X/Ypix')
      for line in cmd.stdout:
	if region_re.match(line):
	 split_line = line.split()
	 x,y = split_line[2], split_line[3] 
   
    self.centroidx = x
    self.centroidy = y
    
    return x,y
    
    
  def define_regions(self):
    """
    Create default source and background .reg files. Will extract an image file if one has not
      already been extracted (Not barycentred).
      Default Regions:
	Source: Circle of radius 20 pixels
	Backgound: Annulus with inner radius 40 pixels and outer radius 60 pixels
    """

    x, y = self.find_centroid()

    f = open(self.path + "/source.reg","w")
    region = '# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal"'\
             + ' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nphysical\ncircle(%s,%s,20)' % (x,y)
    f.write(region) 
    f.close()
    
    f = open(self.path + "/back.reg","w")
    region = '# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal"'\
             + ' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nphysical\nannulus(%s,%s,40,60)' % (x,y)
    f.write(region) 
    f.close()
 
    self.src_region = "source.reg"
    self.back_region = "back.reg"

  def correct_backscal(self):
    """
    Corrects the BACKSCAL keyword in the .pha.grp spectrum file and the background spectrum file.
      Does NOT change the BACKSCAL keyword in the ungrouped source spectrum.
      Assumes the default regions which are:
	Source: Circle of radius 20 pixels
	Background: Annulus with inner radius 40 pixels and outer radius 60 pixels
    """
 
    # Set both BACKSCAL keywords to 1.0 as the 'area' along the strip is the same for
    # a circle of radius 20 pixels and an annulus of 'width' 20 pixels

    hdus = pyfits.open(self.path + self.spectrum, mode='update')
    hdus[1].header['BACKSCAL'] = 1.0
    hdus.close()
    
    hdus = pyfits.open(self.path + self.bg_spectrum, mode='update')
    hdus[1].header['BACKSCAL'] = 1.0
    hdus.close()


  def extract_spectrum(self,infile=None,chan_low=None,chan_high=None,energy_low=None,energy_high=None,grouping=20):
    """
    Extract a spectrum.
      If both the PHA channel limits and energy limits are None will extract entire band.

      Optional Arguments:
        - infile: Name of input file. If None will use baryfile or obsfile.
                  Default = None.
        - chan_low, chan_high: limits of PHA channels to extract.
                               Default = None. 
        - energy_low, energy_high: energy limits in keV to extract. 
                                   Will be superceded by chan_low and chan_high.
                                   Default=None.

    """
    print "Extracting spectrum...\n"

    if chan_high == None or chan_low == None:
      if energy_low == None or energy_high == None:
        chan_low = 0
        chan_high = 1023
      else:
        chan_low = energy_low * 10
        chan_high = energy_high * 10

    x, y = self.find_centroid()

    self.extract(self.obsroot + "_source",infile=infile, events=False, pha=True,\
                   region=self.path + self.src_region, chanlow=chan_low, chanhigh=chan_high)  
    self.extract(self.obsroot + "_back",infile=infile, events=False, pha=True,\
                   region=self.path + self.back_region, chanlow=chan_low, chanhigh=chan_high)  

    cmd = "xrtmkarf outfile=%s%s_source.arf phafile=%s%s_source.pha psfflag=yes srcx=%s srcy=%s clobber=yes"%\
          (self.path, self.obsroot, self.path, self.obsroot, x, y) 
    timed_execute(cmd)  

    grppha_comm = "chkey backfile %s%s_back.pha & chkey ancrfile %s%s_source.arf & chkey respfile"%\
                  (self.path, self.obsroot, self.path, self.obsroot)\
                  + " /exports/scratch/software/CALDB/data/swift/xrt/cpf/rmf/swxpc0to12s6_20070901v011.rmf"\
                  + " & group min %d & exit" % grouping

    cmd = "grppha infile=%s%s_source.pha outfile=%s%s_source.pha.grp clobber=yes comm=\"%s\""%\
          (self.path, self.obsroot, self.path, self.obsroot, grppha_comm)
    timed_execute(cmd)

    self.spectrum = "%s_source.pha.grp" % (self.obsroot)
    self.bg_spectrum = "%s_back.pha" % (self.obsroot)

  def fit_spectrum(self, spectrum=None):
    """
    Fit a blackbody plus power-law model to a spectrum using XSPEC.
      Outputs a .txt file with the spectral parameters as well as .ps and .png plots of the spectrum.

      Optional Arguments:
	- spectrum: filename with path of the grouped pha spectrum to extract. 
		    If None will extract the spectrum defined by self.spectrum. 
		    Default=None.
    """

    if not spectrum:
      spectrum = self.path + self.spectrum
 
    f = open('xspec_script.xcm','w')

    xspec_cmds = "source /homes/borgii/pscholz/.xspec/write_out.tcl\ndata " + spectrum +\
		 "\n@/homes/borgii/pscholz/bin/swiftmonitor/default_fit.xcm\nwrite_out " + self.path +\
                 self.obsroot + "_xspecfit.txt\nplot ldata delchi\nexit"

    f.write(xspec_cmds)
    f.close()

    timed_execute('xspec - xspec_script.xcm')
    cmd = 'mv pgplot.ps ' + self.path + self.obsroot + '_xspecfit.ps'
    timed_execute(cmd)
    timed_execute('rm xspec_script.xcm')

    cmd = 'gs -q -sDEVICE=png16m -r288 -dBATCH -dNOPAUSE -dFirstPage=1 -dLastPage=1 -sOutputFile=' +\
	   self.path + self.obsroot + '_xspecfit.png ' + self.path + self.obsroot + '_xspecfit.ps'
    timed_execute(cmd)
    cmd = 'convert %s -trim %s' % ( self.path + self.obsroot + '_xspecfit.png',self.path + self.obsroot + '_xspecfit.png' )
    timed_execute(cmd)

    self.spec_fit = self.obsroot + "_xspecfit.txt"

  def extract_region(self):
    """
    Extract events from the source and background region files.
    """

    print "Extracting events for default region...\n"
    self.extract(self.obsroot + "_bary_reg",infile=self.path + self.baryfile, events=True, region=self.path + self.src_region)  
    self.extract(self.obsroot + "_bary_bgreg",infile=self.path + self.baryfile, events=True, region=self.path + self.back_region)  
    self.reg_obsfile = self.obsroot + "_bary_reg.evt"
    self.bgreg_obsfile = self.obsroot + "_bary_bgreg.evt"

  def get_countrates(self):
    fits = pyfits.open(self.path + self.reg_obsfile)
    bg_fits = pyfits.open(self.path + self.bgreg_obsfile)

    exposure = float(fits[1].header['EXPOSURE'])
    countrate = float(fits[1].header['NAXIS2']) / exposure
    bg_countrate = float(bg_fits[1].header['NAXIS2']) / exposure
    
    fits.close()
    bg_fits.close()
 
    self.exposure = exposure 
    self.countrate = countrate
    self.bg_countrate = bg_countrate

    return countrate, bg_countrate

  def convert_ds(self):
    cmd = 'swevt2ds.csh ' + self.path + self.reg_obsfile + ' ' + self.path + self.obsroot + '_bary_reg.ds -5' 
    timed_execute(cmd)
    self.dsfile = self.obsroot + '_bary_reg.ds'

  def save(self):
    pickle.dump(self,file(self.path + self.obsroot+'.pkl','w'))
    

# for testing
if __name__ == '__main__':

  ob = Observation('00340573000','/exports/data/pscholz', pulsar='1E1547-5408')
  print ob.ra, ob.dec


