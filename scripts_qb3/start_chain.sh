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

# Decide lists of residues for which a flip WILL be attempted (NOT in secondary structure)
# and will NOT be attempted (IN secondary structure)
# NB: The CalcSecStructure method in MMDB fails if there are any alt confs!
# NB: This QB3 version uses the restoflip executable BEFORE the qFit array
# because we need to decide which nodes (for flip residues) need more memory.
# By contrast, the Stanford version uses the flipornot executable INSIDE the qFit array.
ofn=${pdb%.pdb}_noAlts.pdb
phenix.pdbtools \
    remove_alt_confs=True \
    stop_for_unknowns=False \
    $pdb \
    output.file_name=$ofn
res_numbers_flip=`$QFITDIR/TRUNCATE/restoflip flip $ofn $chain $first $last`
res_numbers_noflip=`$QFITDIR/TRUNCATE/restoflip noflip $ofn $chain $first $last`
num_res_flip=`echo $res_numbers_flip | awk -F"," '{print NF}'`
num_res_noflip=`echo $res_numbers_noflip | awk -F"," '{print NF}'`
echo "Residues for which to try peptide flips: $res_numbers_flip"
echo "Residues for which to NOT try peptide flips: $res_numbers_noflip"

# Launch jobs for residue where we will try a peptide flip, with higher memory requests
# (Requesting 5GB, since that is >4GB; this will give us lots of nodes with 16GB or more)
echo "qsubbing flip residue jobs..."
qsub \
    -S /bin/bash \
    -l arch=linux-x64 \
    -l h_rt=96:0:0 \
    -l netapp=1G \
    -l scratch=1G \
    -l mem_free=5G \
    -cwd \
    -o residue_$first-$last.out -e residue_$first-$last.err -j n \
    -N r_f_${basepdb}${chain}_${first}-${last}_${jobbase} \
    -t 1-$num_res_flip \
    -V \
    $SCRIPTSDIR/residue.sh $pdb $mtz $chain $first $last "flip" $res_numbers_flip
echo "... done qsubbing flip residue jobs"

# Launch jobs for residue where we will NOT try a peptide flip, with lower memory requests
# (Requesting 4GB, since that will give us lots of nodes with 2GB or 4GB)
echo "qsubbing no-flip residue jobs..."
qsub \
    -S /bin/bash \
    -l arch=linux-x64 \
    -l h_rt=96:0:0 \
    -l netapp=1G \
    -l scratch=1G \
    -l mem_free=4G \
    -cwd \
    -o residue_$first-$last.out -e residue_$first-$last.err -j n \
    -N r_nf_${basepdb}${chain}_${first}-${last}_${jobbase} \
    -t 1-$num_res_noflip \
    -V \
    $SCRIPTSDIR/residue.sh $pdb $mtz $chain $first $last "noflip" $res_numbers_noflip
echo "... done qsubbing no-flip residue jobs"

# We're now limiting the residue jobs to 4 days (96 hours), 
# since that should be pleeeenty of time for any given residue...  
# Beyond that we don't really care what qFit finds!

# Queue up finishing the chain once both the flip and no flip sets of residues are done
echo "qsubbing finish_chain job..."
qsub \
    -hold_jid r_f_${basepdb}${chain}_${first}-${last}_${jobbase},r_nf_${basepdb}${chain}_${first}-${last}_${jobbase} \
    -S /bin/bash \
    -l arch=linux-x64,netapp=1G,scratch=1G,mem_free=3G,h_rt=48:0:0 \
    -cwd -o finish_chain_${chain}.out -e finish_chain_${chain}.err -j n \
    -N fc_${basepdb}${chain}_${jobbase} \
    -V \
    $SCRIPTSDIR/finish_chain.sh $pdb $mtz $chain $first $last
echo "... done qsubbing finish_chain job"
