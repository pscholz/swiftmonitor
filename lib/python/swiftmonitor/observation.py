import os, sys, time, shutil, glob
import numpy as np
import subprocess, re, pickle
import astropy.io.fits as pyfits
import swiftmonitor.config
from swiftmonitor import ftools, utils

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
    self.ra, self.dec = None, None
    self.expmap=None

    #auxillary files set by injest_auxil
    self.attfile = None
    self.hdfile = None
    self.teldeffile = None
    self.alignfile = None
  
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
      cmd = subprocess.Popen(['browse_extract_wget.pl', 'table=swiftxrlog', 'position=none', 'param=obsid,' + obsid ],
                               stdout=subprocess.PIPE)
      point_re = re.compile("^\|.*WINDOWED\|pointing(?!_mode)")
    elif self.mode == 'pc':
      cmd = subprocess.Popen(['browse_extract_wget.pl', 'table=swiftxrlog', 'position=none', 'param=obsid,' + obsid ],
                               stdout=subprocess.PIPE)
      point_re = re.compile("^\|.*PHOTON\|pointing(?!_mode)")
  
    date = None
    window = None

    tabledown_re = re.compile('Tableswiftxrlogdoesnotseemtoexist!')
    norows_re = re.compile('returns0rows')

    for line in cmd.stdout:
      line = line.replace(' ','')
      print line
      if tabledown_re.search(line) or norows_re.search(line):
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

    obsdir_url = "http://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/" + date + "/" + obsid + "/"

    return obsfile, orbitfile, obsdir_url

  def download(self):
    """
    Download the observation. 
    """
    print "Querying HEASARC...\n"
  
    # attempt 5 times to retrieve data from heasarc (table often returns empty for no reason) 
    attempts = 0
    while attempts < 5:
        try:
            obsfile, orbitfile, obsurl = self._query_heasarc()
        except TableDownError:
            attempts += 1
            time.sleep(30)
        else:
            break
        
    if attempts >= 5:
        raise TableDownError

    print "Downloading observation...\n"

    cmd = 'wget -q -N -O ' + self.path + '/' + obsfile + " " + obsurl + '/xrt/event/' + obsfile
    wget_time = timed_execute(cmd)
    cmd = 'wget -q -N -O ' + self.path + '/' + orbitfile + " " + obsurl + '/auxil/' + orbitfile
    wget_time += timed_execute(cmd)

    self.orbitfile = orbitfile
    self.obsfile = obsfile
    self.obsroot = obsfile.split('.')[0]

  def download_raw(self):
    print "Querying HEASARC...\n"
  
    # attempt 5 times to retrieve data from heasarc (table often returns empty for no reason) 
    attempts = 0
    while attempts < 5:
        try:
            obsfile, orbitfile, obsurl = self._query_heasarc()
        except TableDownError:
            attempts += 1
            time.sleep(30)
        else:
            break
        
    if attempts >= 5:
        raise TableDownError

    print "Downloading observation (all data)...\n"
    
    os.mkdir(self.path + '/raw/')

    cmd = "wget -q -nH --cut-dirs=6 -r -l0 -R 'index*' -np -erobots=off --retr-symlinks -N -P "\
          + self.path + "/raw/" + " " + obsurl + '/xrt/'
    timed_execute(cmd)
    cmd = "wget -q -nH --cut-dirs=6 -r -l0 -R 'index*' -np -erobots=off --retr-symlinks -N -P "\
          + self.path + "/raw/" + " " + obsurl + '/auxil/'
    timed_execute(cmd)

  def manual_add(self, obsfile, orbitfile):
    """
    Set up for an observation for which the orbitfile and obsfile were downloaded manually. 
      Result is the same as download() without downloading the files or querying HEASARC.
    """
    self.orbitfile = orbitfile
    self.obsfile = obsfile
    self.obsroot = obsfile.split('.')[0]

  def injest_auxil(self):
    """
    Finds auxil files and stores locations as Observation attributes.

      Currently injests: - attitude file (attfile)
                         - hd housekeeping file (hdfile)
                         - teldef file, from CALDB (teldeffile)
                         - alignment file, from CALDB (alignfile)

    Dev: can add auxil files here as needed. Also do some caldb through quzcif?
    """

    # find best attitude file available (uat > pat > sat)
    #attexts = ["uat.fits.gz", "pat.fits.gz", "sat.fits.gz"]
    attexts = ["pat.fits.gz", "sat.fits.gz"]

    for attext in attexts:
        attfile = glob.glob(os.path.join(self.path,'raw/auxil/sw*' + attext))
        if len(attfile):
            self.attfile = attfile[0]
            break
    
    if not self.attfile:
       print "No attitude file not found in auxil files."

    hdfile = glob.glob(os.path.join(self.path,'raw/xrt/hk/sw*hd.hk.gz'))

    if len(hdfile):
        self.hdfile=hdfile[0]
    else:
       print "HD file not found in auxil files."

    event_file = glob.glob(os.path.join(self.path,'raw/xrt/event/sw' + self.obsid + \
                                        'x' + self.mode + '??po_cl.evt.gz'))[0]
    fits = pyfits.open(event_file)
    date_obs = fits[0].header['DATE-OBS']

    date_obs_split = date_obs.strip().strip('\'').split("T")

    self.alignfile = ftools.quzcif('ALIGNMENT','NOW','-', instrument='SC')[0][0]
    self.teldeffile = ftools.quzcif('TELDEF',date_obs_split[0],date_obs_split[1])[0][0]
    

  def preprocess(self, raw_dir, out_dir, xrtpipeline_args=""):
    """
    Preprocess the raw unfiltered data using xrtpipeline.
    """
    self.injest_auxil()
    cmd = 'xrtpipeline indir=%s outdir=%s steminputs=sw%s createexpomap=yes %s' %\
          (raw_dir, out_dir, self.obsid, xrtpipeline_args)
    if self.ra and self.dec:
        cmd += ' srcra=%s srcdec=%s' % (self.ra, self.dec)
    if self.attfile:
        cmd += ' attfile=%s' % self.attfile

    cmd += " %s > %s/xrtpipeline.log" % (xrtpipeline_args, self.path)
    timed_execute(cmd)
 
    event_files = glob.glob(out_dir + "/sw" + self.obsid + "x" + self.mode + "*" + "po_cl.evt")
    orbit_files = glob.glob(raw_dir + "/auxil/sw" + self.obsid + "sao.fits*")
    expmap_files = glob.glob(out_dir + "/sw" + self.obsid + "x" + self.mode + "*" + "po_ex.img")
    
    if not event_files or len(event_files) > 1:
      print "No or more than one cleaned event file output in %s" % out_dir
    if not orbit_files or len(orbit_files) > 1:
      print "No or more than one orbit file exists in %s/auxil/" % raw_dir
    if not expmap_files or len(expmap_files) > 1:
      print "No or more than one exposure map file exists in %s" % out_dir

    shutil.copy(event_files[0], self.path)
    shutil.copy(orbit_files[0], self.path)
    shutil.copy(expmap_files[0], self.path)

    self.obsfile = os.path.basename(event_files[0])
    self.orbitfile = os.path.basename(orbit_files[0])
    self.expmap = os.path.basename(expmap_files[0])
    self.obsroot = self.obsfile.split('.')[0]
    

  def barycentre(self, RA=None, Dec=None):
    """
    Barycentre the observation. If not provided a RA and Dec, will use pulsar RA and Dec
      if a pulsar is defined for the Observation object, otherwise the RA and Dec from
      the fits header will be used.
    """
    self.baryfile = self.obsroot + '_bary.evt'

    outfile = os.path.join(self.path,self.baryfile)
    infile = os.path.join(self.path,self.obsfile)
    orbitfile = os.path.join(self.path,self.orbitfile)

    if RA and Dec:
      ftools.barycentre(infile, outfile, orbitfile, RA=RA, Dec=Dec)
    elif self.pulsar:
      ftools.barycentre(infile, outfile, orbitfile, RA=self.ra, Dec=self.dec)
    else:
      ftools.barycentre(infile, outfile, orbitfile)

  def extract(self,outroot,infile=None,events=True,image=False,pha=False,lc=False,region=None,\
              grade=None,gtifile=None,chanlow=0,chanhigh=1023):
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

    full_outroot = os.path.join(self.path,outroot)

    ftools.extract(full_outroot,infile,events=events,image=image,pha=pha,lc=lc,region=region, \
                   grade=grade,gtifile=gtifile,chanlow=chanlow,chanhigh=chanhigh)
   
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

  def correct_backscal(self, sourcefile=None, bgfile=None, value=1.0):
    """
    Corrects the BACKSCAL keyword in the .pha.grp spectrum file and the background spectrum file.
      Does NOT change the BACKSCAL keyword in the ungrouped source spectrum.
      Assumes the default regions which are:
	Source: Circle of radius 20 pixels
	Background: Annulus with inner radius 40 pixels and outer radius 60 pixels
    """
 
    # Set both BACKSCAL keywords to 1.0 as the 'area' along the strip is the same for
    # a circle of radius 20 pixels and an annulus of 'width' 20 pixels

    if not sourcefile:
      sourcefile = self.path + self.spectrum

    if not bgfile:
      bgfile = self.path + self.bg_spectrum

    hdus = pyfits.open(sourcefile, mode='update')
    hdus[1].header['BACKSCAL'] = value
    hdus.close()
    
    hdus = pyfits.open(bgfile, mode='update')
    hdus[1].header['BACKSCAL'] = value
    hdus.close()


  def extract_spectrum(self,infile=None,chan_low=None,chan_high=None,energy_low=None,energy_high=None,grouping=20,grade=None):
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
        - grouping: the minimum counts per bin for the spectrum.
                    Default=20

    """
    print "Extracting spectrum...\n"

    if chan_high == None or chan_low == None:
      if energy_low == None or energy_high == None:
        chan_low = 0
        chan_high = 1023
      else:
        chan_low = int( energy_low * 100.0 )
        chan_high = int( energy_high * 100.0 )

    x, y = self.find_centroid()

    outroot = self.obsroot 
    if grade:
      outroot += '_g%s' % grade

    self.extract(outroot + "_source",infile=infile, events=False, pha=True,\
                   region=self.path + self.src_region, chanlow=chan_low, chanhigh=chan_high,grade=grade)  
    self.extract(outroot + "_back",infile=infile, events=False, pha=True,\
                   region=self.path + self.back_region, chanlow=chan_low, chanhigh=chan_high,grade=grade)  

    cmd = "xrtmkarf outfile=%s%s_source.arf phafile=%s%s_source.pha psfflag=yes srcx=%s srcy=%s clobber=yes"%\
          (self.path, outroot, self.path, outroot, x, y) 
    if self.expmap:
      cmd += " expofile=%s%s" % (self.path, self.expmap)
    #timed_execute(cmd)  

    pipe = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    xrtmkarf_out = pipe.stdout.read()
    pipe.stdout.close() 
    print xrtmkarf_out

    rmf_re = re.compile("Processing \'(?P<rmf>.*)\.rmf\' CALDB file\.")
    rmf_search = rmf_re.search(xrtmkarf_out)
    if rmf_search:
      rmf = rmf_search.groupdict()['rmf'] + '.rmf'
    else:
      print "ERROR: No rmf filename found from xrtmkarf output."

    if grade and grade != '0':
      print "Grade selection not 0 or default, rmf in 'respfile' keyword may be wrong."

    #if self.mode == 'pc':
    #  rmf = '/exports/scratch/software/CALDB/data/swift/xrt/cpf/rmf/swxpc0to12s6_20010101v013.rmf'
    #elif self.mode == 'wt' and not grade:
    #  rmf = '/exports/scratch/software/CALDB/data/swift/xrt/cpf/rmf/swxwt0to2s6_20010101v014.rmf'
    #elif self.mode == 'wt' and grade == '0':
    #  rmf = '/exports/scratch/software/CALDB/data/swift/xrt/cpf/rmf/swxwt0s6_20010101v014.rmf'


    grppha_comm = "chkey backfile %s%s_back.pha & chkey ancrfile %s%s_source.arf & chkey respfile %s"%\
                  (self.path, outroot, self.path, outroot, rmf)\
                  + " & group min %d & exit" % grouping

    cmd = "grppha infile=%s%s_source.pha outfile=%s%s_source.pha.grp clobber=yes comm=\"%s\""%\
          (self.path, outroot, self.path, outroot, grppha_comm)
    timed_execute(cmd)

    self.spectrum = "%s_source.pha.grp" % (outroot)
    self.bg_spectrum = "%s_back.pha" % (outroot)

  def extract_spectrum_seporb(self,infile=None,chan_low=None,chan_high=None,
                              energy_low=None,energy_high=None,grouping=20,
                              grade=None, badcol_tol=0.0, out_suffix='_seporb'):
    outroot = self.obsroot 
    if grade:
      outroot += '_g%s' % grade

    if infile == None:
        print "Using obsfile as input."
        infile = self.path + self.obsfile

    if not self.attfile or not self.hdfile:
        raise utils.SwiftMonError("No attitude file or no HD file found: attfile=%s, hdfile=%s" % (self.attfile, self.hdfile))
        
    split_files = ftools.split_orbits(infile)

    split_spectra = []
    for split_file in split_files:
        badcol_dist = ftools.locate_bad_columns(self.ra, self.dec, split_file, self.teldeffile, \
                                                self.alignfile, self.attfile)[1][0]

        if badcol_dist < badcol_tol:
            print "\nWARNING: Ignoring orbit %s because %.2f pixels away from bad columns.\n" \
                  "\t Bad column tolerance is set to %.2f pixels.\n" % \
                   (os.path.basename(split_file), badcol_dist, badcol_tol) 

        else:

            ftools.make_expomap(split_file, self.attfile, self.hdfile)
            offaxis_angle = ftools.calc_offaxis_angle(self.ra, self.dec, split_file, self.teldeffile, \
                                                      self.alignfile, self.attfile)

            split_root = os.path.splitext(split_file)[0]
            expomap = split_root + '_ex.img'

            ftools.extract_spectrum(split_root, split_file, expmap=expomap, grade=grade, grouping=None,\
                                    chan_low=chan_low,chan_high=chan_high,energy_low=energy_low,energy_high=energy_high,\
                                    source_region=self.path + self.src_region, back_region=self.path + self.back_region,\
                                    offaxis_angle=offaxis_angle)
                
            split_spectrum = '%s_g%s_source.pha' % (split_root,grade) if grade else '%s_source.pha' % split_root

            split_spectra.append(split_spectrum)
    ftools.add_spectra(split_spectra, self.path + outroot + out_suffix, grouping=grouping)

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
    counts = float(fits[1].header['NAXIS2'])
    bg_counts = float(bg_fits[1].header['NAXIS2'])
    countrate = counts / exposure
    bg_countrate = bg_counts / exposure
 
    fits.close()
    bg_fits.close()
 
    self.counts = counts 
    self.bg_counts = bg_counts 
    self.exposure = exposure 
    self.countrate = countrate
    self.bg_countrate = bg_countrate

    return countrate, bg_countrate

  def convert_ds(self, infile=None, outfile=None, timeres=-5):
    if not infile:
      infile = self.path + self.reg_obsfile

    if not outfile:
      outfile = self.path + self.obsroot +  '_bary_reg.ds'
      self.dsfile = self.obsroot + '_bary_reg.ds'

    cmd = 'swevt2ds.sh ' + infile + ' ' + outfile + ' ' + str(timeres) 
    timed_execute(cmd)

  def save(self):
    pickle.dump(self,file(self.path + self.obsroot+'.pkl','w'))
    

# for testing
if __name__ == '__main__':

  ob = Observation('00340573000','/exports/data/pscholz', pulsar='1E1547-5408')
  print ob.ra, ob.dec


