import swiftmonitor.observation as swift_obs
import sys

ob = swift_obs.load(sys.argv[1])

f = open('xspec_script.xcm','w')

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

