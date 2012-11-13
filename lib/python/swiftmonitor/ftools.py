import os, sys, time

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
