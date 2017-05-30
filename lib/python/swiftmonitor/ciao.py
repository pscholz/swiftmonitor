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
    
    fits = pyfits.open(event_file)
    
    source_ra = fits[1].header['RA_TARG']
    source_dec = fits[1].header['DEC_TARG']
    
    return source_ra,source_dec
    
def make_cc_regions(event_file,source_radius_asec=10,back_rad_1_asec=20,back_rad_2_asec=40,
                    source_fn='source.reg',back_fn='back.reg'):
    """
    Creates the source and background region files for a given event file
    
    Arguments:
        - event_file: the event file to create regions for
    
    Optional Arguments:
        - source_radius_asec: the radius of the source region circle, in asec
        - back_radius_1_asec: the inner radius of the background region annulus, in asec
        - back_radius_2_asec: the outer radius of the background region annulus, in asec
        - source_fn: filename of the source region
        - back_fn: filename of the background region
    """
    
    print('Making source and background region files...\n')
    
    source_ra,source_dec = find_centroid(event_file)
    
    source_reg = region('circle',[source_radius_asec],[source_ra,source_dec],coords='fk5')
    back_reg = region('annulus',[back_rad_1_asec,back_rad_2_asec],
    [source_ra,source_dec],coords='fk5')
    
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
    
    print(specextract)
    
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
    
