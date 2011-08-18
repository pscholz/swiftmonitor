#!/usr/bin/env python

import swiftmonitor.observation as swift_obs
import numpy as np
import os
import sys


ob = swift_obs.load(sys.argv[1])

#os.remove(ob.path + ob.spectrum)
#os.remove(ob.path + ob.reg_obsfile)
#os.remove(ob.path + ob.bg_spectrum)

#spec_split = os.path.splitext(ob.spectrum) # split off the .grp ext
#os.remove(ob.path + spec_split[0])
#spec_split2 = os.path.splitext(spec_split[0]) # split off the .pha ext
#os.remove(ob.path + spec_split2[0] + '.arf')

ob.extract_region()
ob.extract_spectrum()
ob.correct_backscal()
ob.save()
