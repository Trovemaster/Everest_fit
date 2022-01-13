#!/bin/bash -l
#
# Check if our working directory is on the central file server
#

export exec=j-ev-fit_0101_ch.x

export pwd=`pwd`
echo $pwd

export name=`echo $1 | sed -e 's/\.inp//'`
echo $name

export JOB=$3
echo $JOB

if [ -e "$name.o" ]; then
   /bin/rm $name.o
fi

if [ -e "$name.e" ]; then
   /bin/rm $name.e
fi

if [ -e "$name.out" ]; then
  if [ -e "$name.tmp" ]; then
    /bin/rm $name.tmp
  fi
  /bin/mv $name.out $name.tmp
fi

export nproc=32


export jobtype="skylake"
export MEM=372gb
export wclim=$2


echo "Nproc=" $nproc



echo "Nnodes=" 1, "Nproc=" $nproc, " Memory = "$MEM, "jobtype = " $jobtype, "wclimit = " $wclim
echo "Working dir is " $pwd

sbatch -A DIRAC-DP060-CPU --nodes=1  --ntasks=$nproc --time=$wclim:00:00  -J $name -o $name.o -e $name.e   \
     --chdir=$pwd --hint=compute_bound --no-requeue  -p skylake-himem --mem=$MEM \
     $pwd/run_everest_fit.csh $nproc $name $exec $pwd
     

#    --workdir=$pwd --hint=compute_bound --no-requeue --mem=$MEM -p skylake \
     
     
     
