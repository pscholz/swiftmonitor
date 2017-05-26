import os, sys, time, shutil
import numpy as np
import subprocess
import re
#from swiftmonitor import utils
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
    
    print('Preprocessing Chandra observation...\n')
    
    #args = locals()
    
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
            
    chandra_repro()
