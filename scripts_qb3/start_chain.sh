#!/bin/bash

if [ ! -f chains.txt ]; then
    echo "File chains.txt not found."
    exit
fi

chaininfo=`awk '{if (NR=='$SGE_TASK_ID') print $0}' chains.txt`

pdb="$PWD/`basename $1`"
mtz="$PWD/`basename $2`"

chain=`echo $chaininfo | awk '{print $1}'`
first=`echo $chaininfo | awk '{print $2}'`
last=`echo $chaininfo | awk '{print $3}'`
first=$((first+3))
last=$((last-3))

if [ $first -le 0 ]; then
    first=1
fi

if [ ! -d $chain ]; then
    mkdir $chain
fi
cd $chain

tmpstr=`mktemp`
basestr=`basename $tmpstr`
jobbase=`echo $basestr | cut -c 1-10`

basepdb=`basename $pdb .pdb`

# We don't know how long this step will take 
# since we have to launch qfit.sh jobs from within residue.sh,
# which may get stuck in the queue for a long while if a lot of
# other jobs are running,
# which would make residue.sh take a long time while it waits --
# so let's be uber-conservative on runtime here.
qsub \
    -S /bin/bash \
    -l arch=linux-x64,netapp=1G,scratch=1G,mem_free=3G,h_rt=168:0:0 \
    -cwd \
    -o residue_$first-$last.out -e residue_$first-$last.err -j n \
    -N r_${basepdb}${chain}_${first}-${last}_${jobbase} \
    -t $first-$last \
    -V \
    $SCRIPTSDIR/residue.sh $pdb $mtz $chain $first $last

qsub \
    -hold_jid r_${basepdb}${chain}_${first}-${last}_${jobbase} \
    -S /bin/bash \
    -l arch=linux-x64,netapp=1G,scratch=1G,mem_free=3G,h_rt=48:0:0 \
    -cwd -o finish_chain_${chain}.out -e finish_chain_${chain}.err -j n \
    -N fc_${basepdb}${chain}_${jobbase} \
    -V \
    $SCRIPTSDIR/finish_chain.sh $pdb $mtz $chain $first $last
