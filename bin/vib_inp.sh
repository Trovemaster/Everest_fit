#!/bin/bash
#
# Generation of the input file for trove energies
#
export fname=evib
export kmax=$1

echo $fname


cat<<endb> $fname.inp
 &EVEREST
Memo
  6000
Pot
 3
 2 A' Lambda=0
 User ./xyz_pes_koput_cm_ev.so Init='field1.fit'
 2 A' Lambda=1
 User ./xyz_pes_koput_cm_ev.so Init='field2.fit'
 2 A" Lambda=1
 User ./xyz_pes_koput_cm_ev.so Init='field3.fit'
RTStates
 1
 2 3
Atoms
 H Ca O
Masses
1.007276452321 39.951619264818 15.990525980297
Geome
 BLBA
Kmin
 0
Kmax
* 1
 $kmax
r1grid
 sincdvr 100
 2.6 7
r2grid
 sincdvr 100
 1.1 6
xgrid
 leg 120
 0 180
Dim
 10000
*Full
Dav
Nad
 80
Malg
 0
Acc
 1d-3
JDav
Iter
 100
ERef
 2630.00009537940 18953.96173200754 18953.96173200754
* 2628.28242595333 18030.28537245100 18030.28537245100
* 2628.28182542 18030.2887937 18030.2887937
* 2628.28242595
Emax
* 100 100 100
 10000 10000 10000
Prm
 5
omp
 32
end
endb

