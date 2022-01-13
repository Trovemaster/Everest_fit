#!/bin/csh
#!

setenv wdir $pwd

echo $name
echo $nproc
#!
setenv OMP_NUM_THREADS $nproc
setenv MKL_NUM_THREADS $nproc
#!
setenv KMP_LIBRARY turnaround
setenv KMP_AFFINITY disabled

setenv KMP_STACKSIZE 1gb
setenv OMP_NESTED FALSE
#setenv ppn $(uniq -c "$PBS_NODEFILE" | head --lines=1 | sed -e 's/^ *\([0-9]\+\) .*$/\1/g')
#!
hostname
#!
limit
limit datasize unlimited
limit
#!
cd   $wdir
echo -e "Changed directory to `pwd`.\n"

setenv LAUNCH time  ###"dplace -x2"
setenv TMPDIR $wdir
#!
echo "TMPDIR = " $TMPDIR
echo "USER = " $USER
echo "OMP_NUM_THREADS = " $OMP_NUM_THREADS
echo "wdir" $wdir
#!
##setenv wdir  $PBS_O_WORKDIR
echo "wdir" $wdir
echo "OMP_NUM_THREADS=" $OMP_NUM_THREADS
#!
cd $wdir
#!
echo $wdir
#!
if (-e $name.out) then
    /bin/rm $name.out
endif


setenv JOBID $SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"
#!
$LAUNCH $pwd/$exec < $pwd/$name.inp > $pwd/$name.out
#!
echo "DONE"
