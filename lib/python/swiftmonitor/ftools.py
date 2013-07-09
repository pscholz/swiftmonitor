import os, sys, time, shutil
import pyfits
import numpy as np
import subprocess
import re
from swiftmonitor import utils

def execute_cmd(cmd, stdout=sys.stdout, stderr=sys.stderr): 
    """
    Execute the command 'cmd' after logging the command
      to STDOUT.  

      stderr and stdout can be sys.stdout/stderr or any file 
      object to log output to file.

      stdout and stderr are returned if subprocess.PIPE is
      provided for their input parameters. Otherwise will 
      return None.
    """
    sys.stdout.write("\n'"+cmd+"'\n")
    sys.stdout.flush()

    pipe = subprocess.Popen(cmd, shell=True, stdout=stdout, stderr=stderr)
    (stdoutdata, stderrdata) = pipe.communicate()

    retcode = pipe.returncode

    if retcode < 0:
        raise utils.SwiftMonError("Execution of command (%s) terminated by signal (%s)!" % \
                                (cmd, -retcode))
    elif retcode > 0:
        raise utils.SwiftMonError("Execution of command (%s) failed with status (%s)!" % \
                                (cmd, retcode))
    else:
        # Exit code is 0, which is "Success". Do nothing.
        pass

    return (stdoutdata, stderrdata)

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
    extract_time = execute_cmd(cmd)

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

def skyradec_to_det(ra, dec, teldeffile, alignfile, attfile, time):
    """
    Runs and captures the output from pointxform to work out the
      detector and sky pixel coord for the given RA, Dec

      Returns a tuple of (detx, dety, skyx, skyy)
    """

    cmd = 'pointxform fromworld=yes from=\'SKY\' to=\'DET\''
    cmd += ' x=' + str(ra)
    cmd += ' y=' + str(dec)
    cmd += ' teldeffile=' + teldeffile
    cmd += ' alignfile=' + alignfile
    cmd += ' attfile= ' + attfile
    cmd += ' time=' + str(time)

    output = execute_cmd(cmd,stdout=subprocess.PIPE)[0]

    det_result = re.compile("DET.*\[", re.M).search(output)
    sky_result = re.compile("SKY.*\[", re.M).search(output)

    if det_result and sky_result:

        det_str = (det_result.group())[3:-1].split(',')
        sky_str = (sky_result.group())[3:-1].split(',')

        detx, dety = float(det_str[0]), float(det_str[1])
        skyx, skyy = float(sky_str[0]), float(sky_str[1])

        return (detx, dety, skyx, skyy)

    else:
        return (None, None, None, None)

def calc_offaxis_angle(ra, dec, evtfile, teldeffile, alignfile, attfile):
    """
    Calculate offaxis angle of source.
    """
    orbits, gtis, orb_inds = define_orbits(evtfile)

    offaxisang_wsum = 0
    ontime_sum = 0
    for i,orbit in enumerate(orbits):
        tmid = (orbit[0] + orbit[1]) / 2.0
        detx, dety, skyx, skyy = skyradec_to_det(ra, dec, teldeffile, alignfile, attfile, tmid)

        # Work out offaxis angle
        # assume optical axis is at (300.5, 300.5)
        dx = detx - 300.5
        dy = dety - 300.5

        dr = (dx * dx + dy * dy)**0.5
        # offaxis angle in arcmin
        offaxisang = dr * 2.3573 / 60.

        ontime = 0
        for gti in gtis[ orb_inds == i ]:
            ontime += gti[1] - gti[0]

        offaxisang_wsum += offaxisang * ontime
        ontime_sum += ontime

    return offaxisang_wsum / ontime_sum
        


def extract_spectrum(outroot,infile,chan_low=None,chan_high=None,energy_low=None,energy_high=None,\
                     grouping=20,grade=None,expmap=None,source_region=None,back_region=None, \
                     offaxis_angle=None):
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
        - offaxis_angle: offaxis angle (in arcmin) of source to feed to xrtmkarf for 
                    the vignetting correction.
                    Default is to let xrtmkarf do its own calculation.

    """
    print "Extracting spectrum...\n"

    if chan_high == None or chan_low == None:
      if energy_low == None or energy_high == None:
        chan_low = 0
        chan_high = 1023
      else:
        chan_low = int( energy_low * 100.0 )
        chan_high = int( energy_high * 100.0 )

    if offaxis_angle:
        x, y = '-1', '-1'
    else:
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
    if offaxis_angle:
      cmd += " offaxis=%f" % (offaxis_angle)

    xrtmkarf_out = execute_cmd(cmd,stdout=subprocess.PIPE)[0]
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
    execute_cmd(cmd)

    os.remove('temp_source.pha')


def add_spectra(spec_list, outroot, grouping=None):
    """
    Add pha files together. Reimplements addspec ftool.
    """

    back_tmp_spec = []
    src_tmp_spec = []
    tmp_arfs = []
    weights = []

    tstarts = []
    tstops = []
    i = 0

    for spec in spec_list:
        fits = pyfits.open(spec)
        back_fn = fits[1].header['BACKFILE']
        ancr_fn = fits[1].header['ANCRFILE']
        resp_fn = fits[1].header['RESPFILE']
        exposure = fits[1].header['EXPOSURE']
        tstarts.append(fits[1].header['TSTART'])
        tstops.append(fits[1].header['TSTOP'])
        fits.close()

        i += 1
        temp_root = 'temp_spec' + str(i)

        # copy source spectra to cwd
        tmp_spec = temp_root + '.pha' 
        shutil.copy(spec,tmp_spec) 
        src_tmp_spec.append(tmp_spec)
        
        # copy back spectra to cwd
        tmp_spec = temp_root+ '.bak' 
        shutil.copy(back_fn,tmp_spec) 
        back_tmp_spec.append(tmp_spec)

        # copy arf file to cwd
        tmp_arf = temp_root + '.arf' 
        shutil.copy(ancr_fn,tmp_arf) 
        tmp_arfs.append(tmp_arf)

        weights.append(exposure)


    src_math_expr = '+'.join(src_tmp_spec)
    back_math_expr = '+'.join(back_tmp_spec)

    weights = weights / np.sum(weights) # response files weighted by exposure time

    f = open('tmp_arfs.list','w')
    for tmp_arf, weight in zip(tmp_arfs,weights):
        print tmp_arf, weight
        f.write(tmp_arf + ' ' + str(weight) + '\n')
    f.close()

    cmd = "addarf @tmp_arfs.list out_ARF=%s" % outroot + '.arf'
    execute_cmd(cmd)

    cmd = "mathpha expr=%s units=C outfil=temp_final_spec.bak exposure=CALC areascal='%%' backscal='%%' ncomment=0" % (back_math_expr)
    execute_cmd(cmd)

    cmd = "mathpha expr=%s units=C outfil=temp_final_spec.pha exposure=CALC areascal='%%' backscal='%%' ncomment=0" % (src_math_expr)
    execute_cmd(cmd)

    #Run grppha to change the auxfile keys and to do grouping if needed
    grppha_comm = "chkey backfile %s.bak & chkey ancrfile %s.arf & chkey respfile %s"%\
                  (outroot, outroot, resp_fn)
    if grouping:
      grppha_comm += " & group min %d" % grouping
    grppha_comm += " & exit"

    cmd = "grppha infile=temp_final_spec.pha outfile=%s.pha clobber=yes comm=\"%s\""%\
          (outroot, grppha_comm)
    execute_cmd(cmd)

    shutil.copy('temp_final_spec.bak', outroot + '.bak')

    # put TSTART and TSTOP keywords into final spectrum
    out_spec_fits = pyfits.open(outroot + '.pha', mode='update')
    out_spec_fits[0].header['TSTART'] = (np.min(tstarts), 'time start')
    out_spec_fits[0].header['TSTOP'] = (np.max(tstops), 'time stop')
    out_spec_fits.close()

    for temp_fn in src_tmp_spec + back_tmp_spec + tmp_arfs:
        os.remove(temp_fn)

    os.remove('temp_final_spec.bak')
    os.remove('temp_final_spec.pha')
    os.remove('tmp_arfs.list')
    

def define_orbits(event_file, sep_time=1000.0):
    """
    Define the start and end times of orbits for a given event file.
      An orbit is defined a set of GTIs separated by more than sep_time seconds. 

      Returns a list of orbits [(orb_start1,orb_stop1),(orb_start2,orb_stop2),...],
      a list of gtis, and a list of index of orbit in orbit list for each gti.
    """
    fits = pyfits.open(event_file)
    gtis = fits['GTI'].data
    fits.close()

    gtistart = gtis.field('START')
    gtistop  = gtis.field('STOP')

    orb_end_inds = np.where(gtistart[1:] - gtistop[:-1] > sep_time)
    orb_ends = gtistop[orb_end_inds]
    orb_starts = gtistart[1:][orb_end_inds]

    # add start and end time to orbit starts and ends
    orb_ends = np.append(orb_ends,[gtistop[-1]])
    orb_starts = np.append([gtistart[0]],orb_starts)

    orbits = zip(orb_starts, orb_ends)

    orb_inds = np.empty(len(gtis)) # which orbit the GTI belongs to
    orb_inds.fill(len(orbits) + 1) # fill with index that will break things if not replaced
    for i,gti in enumerate(gtis):
        for j,orbit in enumerate(orbits):
           if gti[0] >= orbit[0] and gti[1] <= orbit[1]:
               orb_inds[i] = j

    #import matplotlib.pyplot as plt    
    #times = fits['EVENTS'].data['TIME']
    #plt.plot(times, np.ones(len(times)),'k.')
    #plt.plot(orb_ends, np.ones(len(orb_ends)),'ro')
    #plt.plot(orb_starts, np.ones(len(orb_starts)),'go')
    #plt.plot(gtis['START'],np.ones(len(gtis)),'kx')
    #plt.plot(gtis['STOP'],np.ones(len(gtis)),'kx')
    #plt.show()
    #exit()

    fits.close()
    return (orbits, gtis, orb_inds)

def split_orbits(infile):
    """
    Split an event file into separate event files for each orbit.

      Arguments:
        - infile: Input events file to split.
        - orbitfile: Orbit (.sao.fits) file for observation.
    """

    orbits, gtis, orb_inds = define_orbits(infile)
               
    outfiles = []
    snapshot_num = 0
    
    if np.any(orb_inds == len(orbits) + 1):
        print '\t WARNING: SOME GTIS NOT IN ORBITS!!!'

    for i_orb in range(len(orbits)):
        # make list of gtis in the orbit
        gtifile_str = ''
        for gti in gtis[orb_inds == i_orb]:
            gtifile_str += str(gti[0]) + ' ' + str(gti[1]) + '\n'

        if len(gtifile_str):

            snapshot_num += 1
            tempgti_fn = 'tempGTI_%d.txt' % snapshot_num
            f = open(tempgti_fn, 'w')
            f.write(gtifile_str)
            f.close()

            # extract each orbit and append badpix file   
            outroot = os.path.splitext(infile)[0] + "_s" + str(snapshot_num)  
            extract(outroot, infile=infile, events=True, gtifile=tempgti_fn)

            cmd = "fappend %s[BADPIX] %s.evt" % (infile, outroot)
            execute_cmd(cmd)
            outfiles.append(outroot + '.evt')

            os.remove(tempgti_fn)

    return outfiles
  
def make_expomap(infile, attfile, hdfile, stemout=None, outdir=None):
    """
    Make an exposure map for an event file. Wraps xrtexpomap.

      Arguments:
        - infile: Name of the input cleaned event FITS file.
        - attfile: Name of the input Attitude FITS file.
        - hdfile: Name of the input Housekeeping Header Packets FITS file.

      Optional Arguments:
        - stemout: Stem for the output files.
                   Default is same as input file.
        - outdir: Directory for the output files.
                  Default is the current working dir.
        If both stemout and outdir are None will use same dir and stem 
        as input file. 
    """

    cmd = "xrtexpomap infile=%s attfile=%s hdfile=%s clobber=yes " % (infile, attfile, hdfile)

    inf_path, inf_base = os.path.split(infile)
    if stemout and outdir:
        cmd += "stemout=%s outdir=%s" % (stemout, outdir)
    elif stemout:
        cmd += "stemout=%s outdir=./" % stemout
    elif outdir:
        cmd += "stemout=%s outdir=%s" % (os.path.splitext(inf_base)[0], outdir)
    else:
        cmd += "stemout=%s outdir=%s" % (os.path.splitext(inf_base)[0], inf_path)
    
    execute_cmd(cmd)

class region:
    """
    Class containing info from a region file.
      shape(loc[0],loc[1],dim[0],dim[1],...)
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

def quzcif(codename, date, time, mission='SWIFT', instrument='XRT', detector='-', filter='-', expr='-'):
    """
    Wrapper for quzcif ftool. Queries the CALDB for calibration datasets that match a criteria defined
      by the input parameters. Returns a list of files (outlist) and a list of extentsion numbers
      in those files (extlist).
    """
    cmd = "quzcif "
    cmd += " mission=" + mission
    cmd += " instrument=" + instrument
    cmd += " detector=" + detector
    cmd += " filter=" + filter
    cmd += " expr=" + expr
    cmd += " codename=" + codename
    cmd += " date=" + date
    cmd += " time=" + time

    quzcif_out = execute_cmd(cmd, stdout=subprocess.PIPE)[0]

    outlist = []
    extlist = []
    for line in quzcif_out.split("\n"):
        if len(line):
            sline = line.split()
            outlist.append(sline[0])
            extlist.append(sline[1])
    
    return (outlist, extlist)


