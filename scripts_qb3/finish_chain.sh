#!/bin/bash

pdb=$1
mtz=$2
chain=$3
first=$4
last=$5

rm -f all.pdb
i=$first
while [ $i -le $last ]; do
    # Note: qFit makes Res##.pdb with a capital "R"!
    if [[ ! -f $i/Res$i.pdb && ! -f $i/Res.pdb ]]; then
        echo "qFit output for chain $chain residue $i -- Res$i.pdb -- is missing (!), so taking from input instead for chain $chain final.pdb"
        echo "qFit output for chain $chain residue $i -- Res$i.pdb -- is missing (!), so taking from input instead for chain $chain final.pdb" >> ../qFit.log
        awk '{if( (substr($0,1,4)=="ATOM" || substr($0,1,6)=="HETATM") && substr($0,22,1)=="'$chain'" && substr($0,23,4)+0=='$i' ) print $0}' $pdb >> all.pdb
    else
        echo "qFit output for chain $chain residue $i -- Res$i.pdb -- is present as expected, so using for chain $chain final.pdb"
        awk 'BEGIN {FS=""; OFS=""} {if(substr($0,1,4)=="ATOM" || substr($0,1,6)=="HETATM") {$22="'$chain'"; print $0}}' $i/Res*.pdb >> all.pdb
    fi
    #rm -rf $i # clean up
    i=$((i+1))
done

grep -v END all.pdb > all_.pdb

phenix.reduce -q -trim all_.pdb | grep -v "USER  MOD" > all_trim.pdb

$QFITDIR/TRUNCATE/idmulti 0.1 all_trim.pdb all_idmulti.pdb $chain > fragids.txt

awk '{sub(/ SE /," SD "); print}' all_idmulti.pdb > all_idmulti_sd.pdb

tmpstr=`mktemp`
basestr=`basename $tmpstr`
jobbase=`echo $basestr | cut -c 1-10`

basepdb=`basename $pdb .pdb`
nfrags=`wc fragids.txt | awk '{print $1}'`

# We don't know how long this step will take 
# since we have to launch redokkum_frag.sh jobs from within fragment.sh,
# which may get stuck in the queue for a long while if a lot of
# other jobs are running,
# which would make fragment.sh take a long time while it waits --
# so let's be uber-conservative on runtime here.
qsub \
    -S /bin/bash \
    -l arch=linux-x64,netapp=1G,scratch=1G,mem_free=1G,h_rt=168:0:0 \
    -cwd \
    -o fragment_1-$nfrags.out -e fragment_1-$nfrags.err -j n \
    -N f_${basepdb}${chain}_1-${nfrags}_${jobbase} \
    -t 1-$nfrags \
    -V \
    -sync y \
    $SCRIPTSDIR/fragment.sh \
        ../../${basepdb}_refine_001.pdb \
        ../../${basepdb}_refine_001.mtz \
        $chain \
        ../all_idmulti_sd.pdb

# wait 'til that's done ("-sync y")...

tmpstr=`mktemp`
basestr=`basename $tmpstr`
jobbase=`echo $basestr | cut -c 1-10`

qsub \
    -S /bin/bash \
    -l arch=linux-x64,netapp=1G,scratch=1G,mem_free=1G,h_rt=0:29:0 \
    -cwd \
    -o combine_frags_$chain.out -e combine_frags_$chain.err -j n \
    -N cf_${basepdb}${chain}_${jobbase} \
    -V \
    -sync y \
    $SCRIPTSDIR/combine_frags.sh all_idmulti.pdb $chain

# wait 'til that's done ("-sync y")...

awk '{if(substr($0,22,1)=="'$chain'"&&substr($0,23,4)+0<'$first') print $0}' ../noheader.pdb > final.pdb
grep -v END all_idmulti.pdb >> final.pdb
awk '{if(substr($0,22,1)=="'$chain'"&&substr($0,23,4)+0>'$last') print $0}' ../noheader.pdb >> final.pdb

cd ..  # from chain to structure

#########################  CHECK STATUS OF DAISY CHAIN  #########################

all_chains_done=true
chns=`awk '{printf "%s ", $1 }' chains.txt`
for c in $chns; do
    if [ ! -e $c/final.pdb ]; then
        echo "at least one chain ($c) not done yet"
        all_chains_done=false
        break
    fi
done
if $all_chains_done ; then
    echo "all chains done! ($chain was the last one)"
    
    tmpstr=`mktemp`
    basestr=`basename $tmpstr`
    jobbase=`echo $basestr | cut -c 1-10`
    
    qsub \
        -S /bin/bash \
        -l arch=linux-x64,netapp=1G,scratch=1G,mem_free=3G,h_rt=48:0:0 \
        -cwd -o finish_structure.out -e finish_structure.err -j n \
        -N fs_${basepdb}_${jobbase} \
        -V \
        $SCRIPTSDIR/finish_structure.sh $pdb $mtz
    echo "qsubbed finish_structure.sh at `date`"
fi
