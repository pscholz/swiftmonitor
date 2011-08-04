import swiftmonitor.observation as swift_obs
import sys
import numpy as np
import matplotlib.pyplot as plt
import xspec

ob = swift_obs.load(sys.argv[1])

f = open('xspec_script.xcm','w')

def fit_spectrum(ob):
  spec = xspec.Spectrum(ob.path + ob.spectrum)
  spec.ignore("**-0.5 10.0-**")
  xspec.AllData.ignore("bad")


  model = xspec.Model("ph(bb+po)")

  xspec.Fit.query = 'yes'
  xspec.Fit.perform()
  
  energy_ranges = np.asarray(spec.energies)
  energy = ( energy_ranges[:,0] + energy_ranges[:,1] ) / 2
  energy_err = energy - energy_ranges[:,0]

  x1 = energy
  ex_1 = energy_err
  y1 = spec.values
  ey1 = np.sqrt(spec.variance)

  xspec.AllData.notice("all")

  energy_ranges = np.asarray(spec.energies)
  energy = energy_ranges[:,0] + energy_ranges[:,1] / 2

  x2 = energy
  y2 = model.folded(1)

  plt.errorbar(x1,y1,ey1,ex1,fmt="o")
  plt.plot(x2,y2)  

  #plt.plot(spec.noticed, spec.values, 'ro', spec.noticed, model.folded(1))


"""
xspec_cmds = "source /homes/borgii/pscholz/.xspec/write_out.tcl\ndata " + ob.path + ob.spectrum +\
             "\n@/homes/borgii/pscholz/bin/swiftmonitor/default_fit.xcm\nwrite_out " + ob.path +\
             ob.obsroot + "_xspecfit.txt\nplot ldata delchi\nexit"

f.write(xspec_cmds)
f.close()

swift_obs.timed_execute('xspec - xspec_script.xcm')
cmd = 'mv pgplot.ps ' + ob.path + ob.obsroot + '_xspecfit.ps'
swift_obs.timed_execute(cmd)
swift_obs.timed_execute('rm xspec_script.xcm')


cmd = 'gs -q -sDEVICE=png16m -r288 -dBATCH -dNOPAUSE -dFirstPage=1 -dLastPage=1 -sOutputFile=' +\
       ob.path + ob.obsroot + '_xspecfit.png ' + ob.path + ob.obsroot + '_xspecfit.ps'
swift_obs.timed_execute(cmd)
"""
