#!/bin/bash

pdb_scale=$1
mtz=$2
chain=$3
pdb=$4

fragstart=`awk '{if (NR=='$SGE_TASK_ID') print $1}' fragids.txt`
fragend=`awk '{if (NR=='$SGE_TASK_ID') print $NF}' fragids.txt`

jobdir=`awk '{if (NR=='$SGE_TASK_ID') {print $1 "_" $NF}}' fragids.txt`

if [ ! -d $jobdir ]; then
    mkdir $jobdir
fi

cd $jobdir

#$QFITDIR/qFit/bin/reDokkum_frag $pdb_scale $mtz $chain $fragstart $fragend $pdb > reDokkum_frag.log

# Very rough estimates here to get on nodes with <= 2G, <= 4G, or > 4G, respectively
fraglen=$((fragend-fragstart+1))
memfree="1G" # we're fine with nodes with less memory
if [[ $fraglen -ge 7 ]]; then
    memfree="3G" # request a node with > 4 GB
elif [[ $fraglen -ge 10 ]]; then
    memfree="5G" # request a node with > 4 GB
fi

tmpstr=`mktemp`
basestr=`basename $tmpstr`
jobbase=`echo $basestr | cut -c 1-10`

# If we got the memory requirement ~right, 24 hours should be PLENTY of time for this
qsub \
    -S /bin/bash \
    -l arch=linux-x64 \
    -l netapp=1G \
    -l scratch=1G \
    -l mem_free=${memfree} \
    -l h_rt=24:0:0 \
    -cwd \
    -o redokkum_frag_${chain}_${fragstart}_${fragend}.out \
    -e redokkum_frag_${chain}_${fragstart}_${fragend}.err \
    -j n \
    -N rdf_${chain}_${fragstart}_${fragend}_${jobbase} \
    -V \
    $SCRIPTSDIR/redokkum_frag.sh \
        $pdb_scale \
        $mtz \
        $chain \
        $fragstart \
        $fragend \
        $pdb \
        reDokkum_frag.log
