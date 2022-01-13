#!/bin/bash
#
# Generation of the input file for trove energies
#
export fname=erot
export kmax=$1
export jmax=$2

echo $fname


cat<<endb> $fname.inp
 &EVEREST
%Rot
Memo
 5000
Eck
Int
 1
 ls 2 3 msl=1 msr=1 
 User ./xyz_pes_koput_cm_ev.so Init='field4.fit'
Calc
 jtot 1-$jmax
Emb
 0
Jaco
Eref
 2630.00009537940
Dim
 10000
*roots
* 5
EMax
 10000.
kmax
 $kmax
omp
 32
end
endb
