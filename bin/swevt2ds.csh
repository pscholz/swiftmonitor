#!/bin/csh

##################
#Converts filtered barycentered swift events lists in fits format to ds format in XTE MET time
# -log2timeres is x from 2^x
#################

if ($#argv != 3) then
   echo "Usage: swevt2ds infile outfile log2timeres"
   exit 1
endif

fdump $1'[events]' prhead=no showcol=no showrow=no showunit=no column='TIME' row=- outfile=temp.dat
sw2xtetimes.py temp.dat > temp2.dat
abin_ds -d -z 0 -t $3 -f temp2.dat > temp3.dat
ds_writekey COL0001 "Placeholder for XTE PCA channels" < temp3.dat > $2

rm temp.dat
rm temp2.dat
rm temp3.dat

exit 0
