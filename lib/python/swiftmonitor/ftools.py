import os, sys, time
import pyfits
import numpy as np
import subprocess
import re

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

def extract(outroot,infile,events=True,image=False,pha=False,lc=False,region=None,\
            grade=None,gtifile=None,chanlow=0,chanhigh=1023):
    """
    Wrapper for extractor ftool. If infile is None will use baryfile or obsfile as input.
 
      Arguments:
        - outroot: root of the output files. Will attach extension depending on
                   type of output (.img, .evt, etc.)
        - infile: Input file to extract from.

      Optional Arguments:
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

    args = "'%s[PI = %d : %d]' xcolf=X ycolf=Y tcol=TIME ecol=PI gcol=GRADE xcolh=X ycolh=Y gti='GTI' " %\
            (infile, chanlow, chanhigh)

    if image:
      args += 'imgfile=%s.img ' % outroot
    else:
      args += 'imgfile=NONE '
    if pha:
      args += 'phafile=%s.pha ' % outroot
    else:
      args += 'phafile=NONE '
    if lc:
      args += 'fitsbinlc=%s.lc ' % outroot
    else:
      args += 'fitsbinlc=NONE '
    if events:
      args += 'eventsout=%s.evt ' % outroot
    else:
      args += 'eventsout=NONE '
    if region:
      args += 'regionfile=%s ' % region
    else:
      args += 'regionfile=NONE '
    if gtifile:
      args += 'timefile=%s ' % gtifile
    else:
      args += 'timefile=NONE '
    if grade:
      args += 'gstring=%s ' % grade

    args += 'clobber=yes'
    cmd = 'extractor ' + args 
    extract_time = timed_execute(cmd)

def find_centroid(event_file=None,imagefile=None,force_redo=False,use_max=True):
    """
    Finds centroid of source (hopefully)
      Returns x, y coordinates of centroid in pixels.
    """

    if not imagefile:
      extract("temp", infile=event_file, image=True,events=False)
      imagefile = 'temp.img'

    if use_max:
      fits = pyfits.open(imagefile)
      image = fits[0].data
      fits.close()
    
      y,x = np.unravel_index(np.argmax(image),np.shape(image))
      x += 1
      y += 1

    #else:

    #  if self.centroidx != None and self.centroidy != None and force_redo == False:
    #    return self.centroidx, self.centroidy


    #  if self.pulsar:
    #    cmd = subprocess.Popen(['ximage', '@/homes/borgii/pscholz/bin/swiftmonitor/wt_centroid_radec.xco',\
    #    			  self.path + self.imagefile, str(self.ra), str(self.dec) ], stdout=subprocess.PIPE) 
    #  else:  
    #    cmd = subprocess.Popen(['ximage', '@/homes/borgii/pscholz/bin/swiftmonitor/wt_centroid.xco',\
    #    			  self.path + self.imagefile], stdout=subprocess.PIPE) 

    #  region_re = re.compile('^[ ]*X/Ypix')
    #  for line in cmd.stdout:
    #    if region_re.match(line):
    #     split_line = line.split()
    #     x,y = split_line[2], split_line[3] 
 
    return x,y

def extract_spectrum(outroot,infile,chan_low=None,chan_high=None,energy_low=None,energy_high=None,\
                     grouping=20,grade=None,expmap=None,source_region=None,back_region=None):
    """
    Extract a spectrum.
      If both the PHA channel limits and energy limits are None will extract entire band.

      Arguments:
        - outroot: root of the output files. Will attach extension depending on
                   type of output (.img, .evt, etc.)
        - infile: Input file to extract from.

      Optional Arguments:
        - chan_low, chan_high: limits of PHA channels to extract.
                               Default = None. 
        - energy_low, energy_high: energy limits in keV to extract. 
                                   Will be superceded by chan_low and chan_high.
                                   Default=None.
        - grouping: the minimum counts per bin for the spectrum.
                    Default=20
        - expmap: exposure map file to use in xrtmkarf.
                    Default=None
        - source_region: region file to use to extract source spectrum.
                    Default=None
        - back_region: region file to use to extract background spectrum.
                    Default=None

    """
    print "Extracting spectrum...\n"

    if chan_high == None or chan_low == None:
      if energy_low == None or energy_high == None:
        chan_low = 0
        chan_high = 1023
      else:
        chan_low = int( energy_low * 100.0 )
        chan_high = int( energy_high * 100.0 )

    # FIX THIS 
    x, y = find_centroid(infile)

    if grade:
      outroot += '_g%s' % grade

    extract("temp_source",infile=infile, events=False, pha=True,\
                   region=source_region, chanlow=chan_low, chanhigh=chan_high,grade=grade)  
    extract(outroot + "_back",infile=infile, events=False, pha=True,\
                   region=back_region, chanlow=chan_low, chanhigh=chan_high,grade=grade)  

    cmd = "xrtmkarf outfile=%s_source.arf phafile=temp_source.pha psfflag=yes srcx=%s srcy=%s clobber=yes"%\
          (outroot, x, y) 
    if expmap:
      cmd += " expofile=%s" % (expmap)

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

    grppha_comm = "chkey backfile %s_back.pha & chkey ancrfile %s_source.arf & chkey respfile %s"%\
                  (outroot, outroot, rmf)
    if grouping:
      grppha_comm += " & group min %d" % grouping
    grppha_comm += " & exit"

    cmd = "grppha infile=temp_source.pha outfile=%s_source.pha clobber=yes comm=\"%s\""%\
          (outroot, grppha_comm)
    timed_execute(cmd)

    os.remove('temp_source.pha')

def split_GTI(infile):
    """
    Split an event file into separate event files for each GTI.

      Arguments:
        - infile: Input events file to split.
    """

    fits = pyfits.open(infile)
    rows = float(fits['GTI'].header['NAXIS2'])
    fits.close()

    for i in range(rows):
        
        tempgti_fn = "tempGTI_%d.fits" % (i+1)
        cmd = "fcopy %s[GTI][#row==%d] %s" % (infile, i+1, tempgti_fn)
        timed_execute(cmd)

        outroot = os.path.splitext(infile)[0] + "_s" + str(i+1)  
     
        extract(outroot, infile=infile, events=True, gtifile=tempgti_fn)

        cmd = "fappend %s[BADPIX] %s.evt" % (infile, outroot)
        timed_execute(cmd)

        os.remove(tempgti_fn)
  
def make_expomap(infile, attfile, hdfile, stemout=None, outdir="./"):
    """
    Make an exposure map for an event file. Wraps xrtexpomap.

      Arguments:
        - infile: Name of the input cleaned event FITS file.
        - attfile: Name of the input Attitude FITS file.
        - hdfile: Name of the input Housekeeping Header Packets FITS file.

      Optional Arguments:
        - stemout: Stem for the output files.
                   Default is standard Swift naming convention.
        - outdir: Directory for the output files.
                  Default is the current working dir.
    """

    cmd = "xrtexpomap infile=%s attfile=%s hdfile=%s outdir=%s clobber=yes " % (infile, attfile, hdfile, outdir)
    if stemout:
        cmd += "stemout=%s " % stemout 
    else:
        cmd += "stemout=%s " % os.path.splitext(infile)[0]
    
    timed_execute(cmd)

class region:
    """
    Class containing info from a region file.
    """
    
    def __init__(self, shape, dimensions, location, coords='physical'):
        self.dim = dimensions
        self.shape = shape 
        self.loc = location
        self.coords = coords
        self.region_str = self.get_region_str()

    def get_region_str(self):
        # TODO: add other shapes and assert that shape is available
        if self.shape is 'circle':
            region_str = 'circle(%s,%s,%d)' % (self.loc[0], self.loc[1], self.dim[0])
        if self.shape is 'annulus':
            region_str = 'annulus(%s,%s,%d,%d)' % (self.loc[0], self.loc[1], self.dim[0], self.dim[1])
        return region_str

    def write(self,output_fn):
        f = open(output_fn, 'w')
        region = '# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal"'\
             + ' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
        region += self.coords + "\n"
        region += self.region_str
        f.write(region)
        f.close()
    
    def __str__(self):
        return self.shape + ":\nDimension: " + str(self.dim) + "\nLocation: " + str(self.loc)

def read_region_file(region_fn):
    """
    Reads a region file and returns a 'region' object.
    """
    region_file = open(region_fn,'r')
    lines = region_file.readlines()
    region_file.close()

    for line in lines:
        if line.startswith('#'):
            lines.remove(line)

    coords = lines[1]
    region_str = lines[2]

    reg_split = region_str.split('(')
    reg_type = reg_split[0]
    params = reg_split[1].rstrip().rstrip(')').split(',')
    
    if reg_type == 'circle':
      x, y, r = float(params[0]),float(params[1]),float(params[2])
      reg_obj = region('circle', [r], [x,y], coords)

    elif reg_type == 'annulus':
      x, y, r1, r2 = float(params[0]),float(params[1]),float(params[2]),float(params[3])
      reg_obj = region('annulus', [r1,r2], [x,y], coords)

    return reg_obj

def make_wt_regions(event_file, source_rad, back_rad, source_fn='source.reg', back_fn='back.reg'):
    """
    Create source and background .reg files. Source is a circle of radius=source_rad and ...
    """

    x, y = find_centroid(event_file)

    source_reg = region('circle', [source_rad], [x,y])
    back_reg = region('annulus', [ 100 - back_rad, 100 + back_rad ], [x,y])

    source_reg.write(source_fn)
    back_reg.write(back_fn)

def correct_backscal(source_file, back_file, source_reg_fn, back_reg_fn):
    """
    Corrects the BACKSCAL keyword in the source_file and back_file spectra.
      Assumes the default WT mode regions which are:
	Source: Circle of radius X pixels
	Background: Annulus with inner radius 100-Y pixels and outer radius 100+Y pixels
          where X and Y are determined from the region files.
    """
    source_reg = read_region_file(source_reg_fn)
    back_reg = read_region_file(back_reg_fn)

    if source_reg.shape is 'circle' and back_reg.shape is 'annulus':
        source_backscal = 2*source_reg.dim[0]
        back_backscal = back_reg.dim[1] - back_reg.dim[0] - 1
    else:
        raise ValueError("Regions not 'default' shapes.")

    
    hdus = pyfits.open(source_file, mode='update')
    hdus[1].header['BACKSCAL'] = source_backscal
    hdus.close()
    
    hdus = pyfits.open(back_file, mode='update')
    hdus[1].header['BACKSCAL'] = back_backscal
    hdus.close()
