import os, sys, time, shutil
import astropy.io.fits as pyfits
import numpy as np
import subprocess
import re
from swiftmonitor.utils import region
from ciao_contrib.runtool import *

def preprocess(indir,outdir='',root=None,badpixel=True,process_events=True,
               destreak=True,set_ardlib=True,check_vf_pha=False,pix_adj='default',
               recreate_tg_mask=False,cleanup=True,clobber=True,verbose=1):
    """
    Preprocess data using ciao's 'chandra_repro' (see: 
    http://cxc.harvard.edu/ciao/ahelp/chandra_repro.html)
    
    Arguments:
        - indir: directory which contains all of the downloaded data files
        - outdir: output directory for reprocessed data files; by default creates
                  a "repro" subdirectory beneath the indir directory
                  
    Optional Arguments:
        - root: root for output file names. Proper extensions will be attached
        - badpixel: boolean; whether or not to create a new bad pixel file
        - process_events: boolean; whether or not to create a new level=1 .evt file
        - destreak: boolean; whether or not to destreak the ACIS-8 chip
        - set_ardlib: boolean; whether or not to set ardlib.par with the bad pixel file
        - check_vf_pha: boolean; whether or not to clean ACIS background in very faint data
        - pix_adj: pixel randomization parameter; default|edser|none|randomize
        - recreate_tg_mask: boolean; whether or not to re-run tgdetect2 and tg_create_mask
                            rather than use level 2 region extension. Depends on pile-up
        - cleanup: bolean; whether or not to cleanup intermediate files upon exit
        - clobber: boolean; whether or not to overwrite existing files with same filename
        - verbose: verbose=0 turns off output, =1 prints status messages, =2-5 print 
                   more information about commands being run
    """
    
    print('Preprocessing Chandra observation using chandra_repro...\n')
    
    #args = locals()
    
    chandra_repro.punlearn()
    
    chandra_repro.indir = indir
    chandra_repro.outdir = outdir
    chandra_repro.root = root
    chandra_repro.badpixel = badpixel
    chandra_repro.process_events = process_events
    chandra_repro.destreak = destreak
    chandra_repro.set_ardlib = set_ardlib
    chandra_repro.check_vf_pha = check_vf_pha
    chandra_repro.pix_adj = pix_adj
    chandra_repro.recreate_tg_mask = recreate_tg_mask
    chandra_repro.cleanup = cleanup
    chandra_repro.clobber = clobber
    chandra_repro.verbose = verbose
    
    #print(chandra_repro)
            
    chandra_repro()
    
def find_centroid(event_file):
    """
    Finds centroid of the source by using the target's RA and DEC from the level 2 evt file.
    
    Arguments:
        - event_file: the event file of which to find the centroid
    """
    
    print('Finding the centroid of the event file...\n')
    
    make_img(event_file,clobber=True)
    
    fits = pyfits.open('temp.fits')
    
    #Previously used the RA and DEC headers to find the centre, now trying a more nuanced
    #max pixel value method
    
    #source_ra = fits[1].header['RA_TARG']
    #source_dec = fits[1].header['DEC_TARG']
    
    #return source_ra,source_dec
    
    data = fits[0].data
    
    #As the data from make_img is 1024x1024 based on the centre of the image, use modulo
    #arithmetic to find the physical x and y coordinates
    
    argmax = np.argmax(data)
    
    x = argmax%1024 + 3584
    y = int(argmax/1024) + 3584
    
    return x,y
    
def make_cc_regions(event_file,source_radius_asec=30,back_rad_1_asec=60,back_rad_2_asec=90,
                    source_fn='source.reg',back_fn='back.reg'):
    """
    Creates the source and background region files for a given event file
    
    Arguments:
        - event_file: the event file to create regions for
    
    Optional Arguments:
        - source_radius_asec: the radius of the source region circle, in pixels
        - back_radius_1_asec: the inner radius of the background region annulus, in pixels
        - back_radius_2_asec: the outer radius of the background region annulus, in pixels
        - source_fn: filename of the source region
        - back_fn: filename of the background region
    """
    
    print('Making source and background region files...\n')
    
    x,y = find_centroid(event_file)
    
    source_reg = region('circle',[source_radius_asec],[x,y],coords='physical')
    back_reg = region('annulus',[back_rad_1_asec,back_rad_2_asec],
    [x,y],coords='physical')
    
    source_reg.write(source_fn)
    back_reg.write(back_fn)

def extract_spectrum(event_file,outroot,source_reg_fn,back_reg_fn,
                     bkgresp=True,weight=False,weight_rmf=False,correctpsf=True,
                     combine=False,grouptype='NUM_CTS',binspec=15,bkg_grouptype='NONE',
                     energy='0.3:11.0:0.01',channel='1:1024:1',energy_wmap='300:2000',
                     binarfcorr='1',binwmap='tdet=8',binarfwmap='1',clobber=False,
                     verbose=0):
    """
    Extracts the spectrum using ciao's 'specextract' (see:
    http://cxc.harvard.edu/ciao/ahelp/specextract.html)
    
    Arguments:
        - event_file: the event file to extract a spectrum from
        - outroot: the prefix to the newly created files
        - source_reg_fn: the filename of the source region
        - back_reg_fn: the filename of the background region
        
    Optional Arguments:
        - bkgresp: boolean; whether or not to write background responses
        - weight: boolean; whether or not to generate weighted ARFs
        - weight_rmf: boolean; whether or not to generate weighted RMFs
        - correctpsf: boolean; whether or not to apply a psf correction to unweighted ARFs
        - combine: boolean; whether or not to combine output spectra and responses
        - grouptype: source spectrum grouping type (see ciao page for detailed info)
        - binspec: source spectrum grouping specification
        - bkg_grouptype: background spectrum grouping type
        - energy: energy grid in keV: (start value):(stop value):(step size)
        - channel: RMF binning specification in pixels: (start):(stop):(step)
        - energy_wmap: energy range for wmaps
        - binarfcorr: detectro pixel bin size for arfcorr psf correction
        - binwmap: binning factor for RMF wmaps
        - binarfwmap: binning factor for ARF wmaps
        - clobber: boolean; whether or not existing output files should be overwritten
        - verbose: integer from 0-5, level of displaying diagnostic messages
        
    """
    
    print('Extracting a spectrum using specextract...\n')
    
    specextract.punlearn()
    
    specextract.infile = event_file + '[sky=region(' + source_reg_fn + ')]'
    specextract.bkgfile = event_file + '[sky=region(' + back_reg_fn + ')]'
    specextract.outroot = outroot
    specextract.bkgresp = bkgresp
    specextract.weight = weight
    specextract.weight_rmf = weight_rmf
    specextract.correctpsf = correctpsf
    specextract.combine = combine
    specextract.grouptype = grouptype
    specextract.binspec = binspec
    specextract.bkg_grouptype = bkg_grouptype
    specextract.energy = energy
    specextract.channel = channel
    specextract.energy_wmap = energy_wmap
    specextract.binarfcorr = binarfcorr
    specextract.binwmap = binwmap
    specextract.binarfwmap = binarfwmap
    specextract.clobber = clobber
    specextract.verbose = verbose
    
    #print(specextract)
    
    specextract()

def correct_cc_backscal(source_file,back_file,source_reg_fn,back_reg_fn):
    """
    Corrects the BACKSCAL keyword in the source_file and back_file spectra.
    Assumes the default CC regions: a circle for the source, and an annulus for the background.
    As the observation is CC mode, sets the background scaling to 1 and the source scaling to
    the ratio of the background's one dimensional extent and the source's one dimensional extent.
    
    Arguments:
        - source_file: the source spectrum file to correct the scaling of
        - back_file: the background spectrum file to correct the scaling of
        - source_reg_fn: the filename of the source .reg file
        - back_reg_fn: the filename of the background .reg file
    """
    
    source_reg = read_region_file(source_reg_fn)
    back_reg = read_region_file(back_reg_fn)
    
    if source_reg.shape is 'circle' and back_reg.shape is 'annulus':
        source_backscal = source_reg.dim[0]/(back_reg.dim[1]-back_reg.dim[0])
    else:
        raise ValueError("Regions are not 'default' shapes.")
        
    fits = pyfits.open(source_file,mode='update')
    fits[1].header['BACKSCAL'] = source_backscal
    fits.close()
    
    fits = pyfits.open(source_file,mode='update')
    fits[1].header['BACKSCAL'] = 1
    fits.close()

def barycentre(infile,orbitfile,outfile,ra=None,dec=None,refframe='INDEF',clobber=False):
    """
    Barycentre the observation using ciao's 'axbary' (see:
    http://cxc.harvard.edu/ciao/ahelp/axbary.html)
    
    Arguments:
        - infile: the name of the input event file to barycentre
        - outfile: the output filename. Set clobber=True to do corrections in situ
        - orbitfile: the input orbit ephemeris file: of the form "orbitf051004864N002_eph1.fits."
    
    Optional Arguments:
        - ra: RA in decimal degrees for use in barycentre corrections
        - dec: DEC in decimal degrees for use in barycentre corrections
        - refframe: reference frame to be used for corrections - 'FK5' and 'ICRS' are allowed values
        - clobber: boolean; whether or not to overwrite the existing file.
    """
    
    print('Barycentre the observation...\n')
    
    #This fixes the issue where axbary would not work through scripts
    
    from ciao_contrib import runtool
    runtool._no_par_file.append("axbary")
    axbary = runtool.make_tool("axbary")
    
    axbary.punlearn()
    
    axbary.infile = infile
    axbary.orbitfile = orbitfile
    axbary.outfile = outfile
    axbary.refframe = refframe
    axbary.clobber = clobber
    
    if ra:
    
        axbary.ra = ra
        
    if dec:
    
        axbary.dec = dec
        
    axbary()
    
def make_img(infile,outfile='temp.fits',option='image',clobber=False):
    """
    Create a temporary image file for use in finding the centroid: makes a 1024x1024 square
    of the centre of the image.
    
    Arguments:
        - infile: the name of the input event file to make a temporary image file of
    Optional Arguments:
        - outfile: the name of the output file to be created
        - option: defaults to image
        - clobber: boolean; whether or not to overwrite a file if outfile exists
    """
    
    print('Making 1024x1024 central image...\n')
    
    dmcopy.punlearn()
    
    dmcopy.infile = infile + '[EVENTS][bin x=3584:4608,y=3584:4608]'
    dmcopy.outfile = outfile
    dmcopy.option = option
    dmcopy.clobber = clobber
    
    dmcopy()
    
def make_expomap(asphistfile,outfile,instmapfile,xygrid,normalize=True,useavgaspect=False,
                 geompar='geom',verbose=0,clobber=False):
    """
    Creates a Chandra imaging exposure map using ciao's mkexpmap. (see:
    http://cxc.harvard.edu/ciao/ahelp/mkexpmap.html)
    
    Arguments:
        - asphistfile: the input aspect histogram file (can be made through asphist)
        - outfile: the name of the output exposure map file
        - instmapfile: the name of the input instrument map file (can be made through mkinstmap)
        - xygrid: the string specifying a region in sky coordinates and resolution of output
        
    Optional Arguments:
        - normalize: selects the physical units of the output exposure map: if true, units of
                     cm**2 counts/photon, if false, units of s * cm**2 counts/photon
        - useavgaspect: if set to true, only the average aspect pointing used to derive exposure map
        - geompar: name of the Pixlib Geometry parameter file
        - verbose: integer from 0-5, level of displaying diagnostic messages
        - clobber: boolean; whether or not existing output files will be overwritten
    """
    
    print('Making an exposure map...\n')
    
    mkexpmap.punlearn()
    
    mkexpmap.asphistfile = asphistfile
    mkexpmap.outfile = outfile
    mkexpmap.instmapfile = instmapfile
    mkexpmap.xygrid = xygrid
    mkexpmap.normalize = normalize
    mkexpmap.useavgaspect = useavgaspect
    mkexpmap.geompar = geompar
    mkexpmap.verbose = verbose
    mkexpmap.clobber = clobber
    
    mkexpmap()
    
def make_3D_hist(infile,outfile,evtfile,dtffile=None,verbose=0,clobber=False):
    """
    Make a 3D histogram of duration vs pointing offset and roll offset using ciao's 'asphist' (see:
    http://cxc.harvard.edu/ciao/ahelp/asphist.html)
    
    Arguments:
        - infile: the name of the input aspect solution file
        - outfile: the name of the file to write
        - evtfile: the input event file name
        
    Optional Arguments:
        - dtffile: the input DTF file name
        - verbose: integer from 0-5, level of displaying diagnostic messages
        - clobber: boolean; whether or not existing output files will be overwritten
    """
    
    print('Making a 3D histogram of duration vs. pointing offset and roll offset...\n')
    
    asphist.punlearn()
    
    asphist.infile = infile
    asphist.outfile = outfile
    asphist.evtfile = evtfile
    asphist.dtffile = dtffile
    asphist.verbose = verbose
    asphist.clobber = clobber
    
    asphist()
    
def instrument_map(outfile,pixelgrid,obsfile,detsubsys='ACIS-S3',spectrumfile=None,
                   monoenergy=1.0,grating=None,maskfile=None,verbose=0,clobber=False):
    """
    Generates a Chandra instrument map using ciao's 'mkinstmap' (see:
    http://cxc.harvard.edu/ciao/ahelp/mkinstmap.html)
    
    Arguments:
        - outfile: name of the output instrument map file
        - pixelgrid: string specifying the detecter region and resolution of instrument map
        - obsfile: name of the event fits file
        - detsubsys: specifies detector subsystem
    
    Optional Arguments:
        - spectrumfile: input spectral weights file
        - monoenergy: the energy, in keV, for which the instrument map will be computed
        - grating: whether to compute a zeroth order instrument map: may be 'HETG' or 'LETG'
        - maskfile: the mask fits file for the observation
        - verbose: integer from 0-5, level of displaying diagnostic messages
        - clobber: boolean; whether or not existing output files will be overwritten
    """
    
    print('Generating an instrument map...\n')
    
    mkinstmap.punlearn()
    
    mkinstmap.outfile = outfile
    #mkinstmap.spectrumfile = spectrumfile
    mkinstmap.monoenergy = monoenergy
    mkinstmap.pixelgrid = pixelgrid
    mkinstmap.obsfile = obsfile
    mkinstmap.detsubsys = detsubsys
    #mkinstmap.grating = grating
    #mkinstmap.maskfile = maskfile
    mkinstmap.verbose = verbose
    mkinstmap.clobber = clobber
    
    mkinstmap()
    
